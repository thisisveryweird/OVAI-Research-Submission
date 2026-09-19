#!/usr/bin/env python3
"""Leakage-safe feature-selection study for the verified GSE94533 cohort.

This script deliberately stays inside the current GSE94533 cohort. It does not
load GSE106817, CPTAC, TCGA, or any other external dataset. Feature selection,
missing-value handling, scaling, and model fitting are repeated independently
inside every stratified training fold.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler


# ---------------------------------------------------------------------------
# Repository paths and immutable study settings
# ---------------------------------------------------------------------------

PROJECT = Path(__file__).resolve().parents[1]
MATRIX_PATH = PROJECT / "data" / "processed" / "gse94533_mirna_tpm.csv.gz"
METADATA_PATH = PROJECT / "data" / "processed" / "gse94533_metadata.csv"
REPORT_PATH = PROJECT / "reports" / "feature_selection_stability.md"
RESULTS_PATH = PROJECT / "reports" / "feature_selection_results.csv"
STABILITY_PATH = PROJECT / "reports" / "feature_selection_stability.csv"

SEED = 20260919
N_SPLITS = 5
SIGNATURE_SIZES = (5, 10, 14, 20, 30, 50)
SELECTION_METHODS = ("l1", "elastic_net")

SelectionMethod = Literal["l1", "elastic_net"]


@dataclass(frozen=True)
class Task:
    """A prespecified binary endpoint evaluated with the same fold protocol."""

    name: str
    frame_mask_description: str


@dataclass(frozen=True)
class FoldResult:
    """Metrics and selected features from one held-out fold."""

    task: str
    method: str
    signature_size: int
    fold: int
    roc_auc: float
    balanced_accuracy: float
    sensitivity: float
    specificity: float
    accuracy: float
    selected_features: tuple[str, ...]


@dataclass(frozen=True)
class StudyData:
    """Patient-level matrix and metadata used by this study."""

    frame: pd.DataFrame
    feature_names: tuple[str, ...]


def load_verified_data() -> StudyData:
    """Load and validate only the local, verified GSE94533 patient matrix."""

    metadata = pd.read_csv(METADATA_PATH)
    matrix = pd.read_csv(MATRIX_PATH)

    if len(metadata) != 179 or len(matrix) != 179:
        raise RuntimeError(
            "Phase 6 requires the verified 179-patient GSE94533 cohort. "
            f"Metadata={len(metadata)}, matrix={len(matrix)}."
        )

    if metadata["sample_id_normalized"].duplicated().any():
        raise RuntimeError("Duplicate GSE94533 patient IDs in metadata.")
    if matrix["sample_id_normalized"].duplicated().any():
        raise RuntimeError("Duplicate GSE94533 patient IDs in matrix.")

    frame = metadata[
        [
            "sample_id_normalized",
            "diagnosis",
            "binary_label",
            "stage",
        ]
    ].merge(
        matrix,
        on="sample_id_normalized",
        how="inner",
        validate="one_to_one",
    )

    feature_names = tuple(
        column for column in matrix.columns if column.startswith("hsa-")
    )
    if len(feature_names) < max(SIGNATURE_SIZES):
        raise RuntimeError(
            f"Expected at least {max(SIGNATURE_SIZES)} miRNA features; "
            f"found {len(feature_names)}."
        )

    # These are data-integrity checks, not feature selection. No feature is
    # discarded based on the complete dataset before cross-validation.
    values = frame[list(feature_names)].apply(pd.to_numeric, errors="coerce")
    if not np.isfinite(values.to_numpy(dtype=float)).all():
        raise RuntimeError("The verified GSE94533 matrix contains non-finite values.")
    if (values.to_numpy(dtype=float) < 0).any():
        raise RuntimeError("The verified GSE94533 TPM matrix contains negative values.")

    diagnosis_counts = frame["diagnosis"].value_counts().to_dict()
    expected_counts = {"Cancer": 96, "Benign": 46, "Borderline": 22, "Control": 15}
    if diagnosis_counts != expected_counts:
        raise RuntimeError(
            "Unexpected GSE94533 diagnosis composition: "
            f"{diagnosis_counts}; expected {expected_counts}."
        )

    return StudyData(frame=frame, feature_names=feature_names)


def task_subset(data: StudyData, task_name: str) -> tuple[pd.DataFrame, np.ndarray]:
    """Return one prespecified endpoint without relabeling borderline tumors."""

    frame = data.frame
    if task_name == "all_cancer_vs_non_cancer":
        subset = frame.copy()
        labels = subset["binary_label"].astype(int).to_numpy()
    elif task_name == "early_stage_vs_benign_control":
        # Borderline tumors are intentionally excluded, not silently treated as
        # benign or control. Positive cases are only Stage I/II cancer.
        keep = frame["diagnosis"].eq("Cancer") & frame["stage"].eq("Stage I/II")
        keep |= frame["diagnosis"].isin(["Benign", "Control"])
        subset = frame.loc[keep].copy()
        labels = (
            subset["diagnosis"].eq("Cancer")
            & subset["stage"].eq("Stage I/II")
        ).astype(int).to_numpy()
    else:
        raise ValueError(f"Unknown Phase 6 task: {task_name}")

    if labels.sum() == 0 or labels.sum() == len(labels):
        raise RuntimeError(f"Task {task_name} does not contain both classes.")
    return subset, labels


def selector_model(method: SelectionMethod) -> LogisticRegression:
    """Create the fold-local sparse selector with a fixed, prespecified setup."""

    if method == "l1":
        return LogisticRegression(
            penalty="l1",
            solver="saga",
            C=0.1,
            class_weight="balanced",
            max_iter=5000,
            random_state=SEED,
        )
    return LogisticRegression(
        penalty="elasticnet",
        l1_ratio=0.5,
        solver="saga",
        C=0.1,
        class_weight="balanced",
        max_iter=5000,
        random_state=SEED,
    )


def preprocessing_pipeline() -> Pipeline:
    """Create preprocessing that will be fitted only on a training fold."""

    return Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
    ])


def select_features(
    x_train: pd.DataFrame,
    y_train: np.ndarray,
    feature_names: tuple[str, ...],
    method: SelectionMethod,
    signature_size: int,
) -> tuple[Pipeline, tuple[str, ...], np.ndarray]:
    """Fit preprocessing and sparse selection on one training fold only.

    The selector is fit after the fold-local imputer and scaler. If the sparse
    model produces fewer non-zero coefficients than the requested signature
    size, remaining features are filled by descending absolute coefficient
    magnitude. This keeps every requested signature size defined while retaining
    the selector's ranking rather than selecting from held-out observations.
    """

    preprocessor = preprocessing_pipeline()
    x_train_processed = preprocessor.fit_transform(x_train)
    selector = selector_model(method)
    selector.fit(x_train_processed, y_train)
    coefficients = np.abs(selector.coef_[0])
    ranked_indices = np.argsort(-coefficients, kind="stable")
    selected_indices = ranked_indices[:signature_size]
    selected_names = tuple(feature_names[index] for index in selected_indices)
    return preprocessor, selected_names, selected_indices


def metrics(
    y_true: np.ndarray,
    probability: np.ndarray,
    threshold: float = 0.5,
) -> dict[str, float]:
    """Calculate threshold-free and thresholded metrics for one fold."""

    prediction = (probability >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, prediction, labels=[0, 1]).ravel()
    return {
        "roc_auc": float(roc_auc_score(y_true, probability)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, prediction)),
        "sensitivity": float(tp / (tp + fn)) if tp + fn else np.nan,
        "specificity": float(tn / (tn + fp)) if tn + fp else np.nan,
        "accuracy": float(accuracy_score(y_true, prediction)),
    }


def evaluate_task(
    data: StudyData,
    task: Task,
    x: pd.DataFrame,
    y: np.ndarray,
) -> tuple[list[FoldResult], dict[tuple[str, int], dict[str, list[int]]]]:
    """Run stratified OOF evaluation for all methods and signature sizes."""

    splitter = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED)
    fold_results: list[FoldResult] = []
    selections: dict[tuple[str, int], dict[str, list[int]]] = {
        (method, size): {feature: [] for feature in data.feature_names}
        for method in SELECTION_METHODS
        for size in SIGNATURE_SIZES
    }

    for fold, (train_indices, test_indices) in enumerate(splitter.split(x, y), start=1):
        x_train = x.iloc[train_indices]
        x_test = x.iloc[test_indices]
        y_train = y[train_indices]
        y_test = y[test_indices]

        for method in SELECTION_METHODS:
            for signature_size in SIGNATURE_SIZES:
                preprocessor, selected_names, selected_indices = select_features(
                    x_train,
                    y_train,
                    data.feature_names,
                    method,  # type: ignore[arg-type]
                    signature_size,
                )
                for feature_index in selected_indices:
                    selections[(method, signature_size)][data.feature_names[feature_index]].append(fold)

                train_processed = preprocessor.transform(x_train)[:, selected_indices]
                test_processed = preprocessor.transform(x_test)[:, selected_indices]
                classifier = LogisticRegression(
                    penalty="l2",
                    C=0.3,
                    solver="liblinear",
                    class_weight="balanced",
                    max_iter=5000,
                    random_state=SEED,
                )
                classifier.fit(train_processed, y_train)
                probability = classifier.predict_proba(test_processed)[:, 1]
                fold_metrics = metrics(y_test, probability)
                fold_results.append(FoldResult(
                    task=task.name,
                    method=method,
                    signature_size=signature_size,
                    fold=fold,
                    selected_features=selected_names,
                    **fold_metrics,
                ))

    return fold_results, selections


def result_frame(results: list[FoldResult]) -> pd.DataFrame:
    """Aggregate fold metrics without selecting a winning result post hoc."""

    rows = []
    for method in SELECTION_METHODS:
        for signature_size in SIGNATURE_SIZES:
            for task_name in sorted({result.task for result in results}):
                group = [
                    result for result in results
                    if result.method == method
                    and result.signature_size == signature_size
                    and result.task == task_name
                ]
                rows.append({
                    "task": task_name,
                    "selection_method": method,
                    "signature_size": signature_size,
                    "folds": len(group),
                    "roc_auc_mean": np.mean([item.roc_auc for item in group]),
                    "roc_auc_sd": np.std([item.roc_auc for item in group], ddof=1),
                    "balanced_accuracy_mean": np.mean([item.balanced_accuracy for item in group]),
                    "sensitivity_mean": np.mean([item.sensitivity for item in group]),
                    "specificity_mean": np.mean([item.specificity for item in group]),
                    "accuracy_mean": np.mean([item.accuracy for item in group]),
                })
    return pd.DataFrame(rows)


def stability_frame(
    data: StudyData,
    selections_by_task: dict[str, dict[tuple[str, int], dict[str, list[int]]]],
) -> pd.DataFrame:
    """Create fold-frequency rows for every feature, task, method, and size."""

    rows = []
    for task_name, selections in selections_by_task.items():
        for (method, signature_size), feature_folds in selections.items():
            for feature in data.feature_names:
                folds = feature_folds[feature]
                rows.append({
                    "task": task_name,
                    "selection_method": method,
                    "signature_size": signature_size,
                    "feature": feature,
                    "selected_folds": len(folds),
                    "selection_frequency": len(folds) / N_SPLITS,
                    "folds": ",".join(str(fold) for fold in folds),
                })
    return pd.DataFrame(rows)


def write_report(
    data: StudyData,
    results: pd.DataFrame,
    stability: pd.DataFrame,
) -> None:
    """Write a human-readable stability report with scientific guardrails."""

    lines = [
        "# Phase 6 Feature-Selection Stability",
        "",
        "This report was generated exclusively from the verified local GSE94533 cohort. No external dataset was loaded, and no feature was selected using held-out observations.",
        "",
        "## Cohort and Protocol",
        "",
        f"The study contains 179 patients and {len(data.feature_names)} measured miRNA features: 96 cancer, 46 benign, 22 borderline, and 15 controls.",
        "",
        "Each task uses shuffled stratified 5-fold cross-validation with a fixed seed. In every training fold, median imputation, standardization, L1 or Elastic-Net selection, and the final classifier are fit using training patients only. The requested signature sizes are 5, 10, 14, 20, 30, and 50.",
        "",
        "The early-stage endpoint is Stage I/II cancer versus benign/control. Borderline tumors are excluded from that endpoint and are not silently relabeled.",
        "",
        "## Cross-Validated Metrics",
        "",
        "Metrics below are means across held-out folds. They are descriptive OOF research estimates, not independent validation results.",
        "",
        "| Task | Selector | Features | ROC-AUC | Balanced accuracy | Sensitivity | Specificity | Accuracy |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in results.sort_values(["task", "selection_method", "signature_size"]).itertuples(index=False):
        lines.append(
            f"| {row.task} | {row.selection_method} | {row.signature_size} | "
            f"{row.roc_auc_mean:.3f} | {row.balanced_accuracy_mean:.3f} | "
            f"{row.sensitivity_mean:.3f} | {row.specificity_mean:.3f} | {row.accuracy_mean:.3f} |"
        )

    lines.extend([
        "",
        "## Feature Stability",
        "",
        "Selection frequency is the number of outer folds in which a feature was selected divided by five. High frequency indicates reproducibility under these folds; it is not proof of biological causality or external generalization.",
        "",
    ])
    for task_name in sorted(stability["task"].unique()):
        lines.append(f"### {task_name}")
        lines.append("")
        for method in SELECTION_METHODS:
            lines.append(f"#### {method}")
            lines.append("")
            for signature_size in SIGNATURE_SIZES:
                top = stability[
                    (stability["task"] == task_name)
                    & (stability["selection_method"] == method)
                    & (stability["signature_size"] == signature_size)
                ].sort_values(
                    ["selection_frequency", "feature"],
                    ascending=[False, True],
                ).head(15)
                compact = "; ".join(
                    f"{row.feature} ({int(row.selected_folds)}/5)"
                    for row in top.itertuples(index=False)
                    if row.selected_folds > 0
                )
                lines.append(f"- Top {signature_size}: {compact or 'No features selected.'}")
            lines.append("")

    lines.extend([
        "## Interpretation Guardrails",
        "",
        "- The 14-feature signature size is included as a sparsity point, not privileged as the historical signature.",
        "- Feature stability is estimated within GSE94533 only. It must not be interpreted as stability across GSE106817, CPTAC, TCGA, or other platforms.",
        "- Model and signature-size comparisons should consider ROC-AUC, balanced accuracy, sensitivity, specificity, and fold variability together; accuracy alone is insufficient.",
        "- Borderline tumors remain represented in the all-cancer-versus-non-cancer endpoint but are explicitly excluded from the early-stage-versus-benign/control endpoint.",
        "",
        "Machine-readable fold summaries: `reports/feature_selection_results.csv`.",
        "Machine-readable stability table: `reports/feature_selection_stability.csv`.",
    ])
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines) + "\n")


def main() -> None:
    """Run the complete leakage-safe Phase 6 study."""

    data = load_verified_data()
    all_results: list[FoldResult] = []
    all_selections: dict[str, dict[tuple[str, int], dict[str, list[int]]]] = {}

    for task in (
        Task("all_cancer_vs_non_cancer", "96 cancer versus all 83 non-cancer"),
        Task("early_stage_vs_benign_control", "52 Stage I/II cancer versus 61 benign/control"),
    ):
        subset, labels = task_subset(data, task.name)
        x = subset[list(data.feature_names)]
        task_results, task_selections = evaluate_task(data, task, x, labels)
        all_results.extend(task_results)
        all_selections[task.name] = task_selections

    result_table = result_frame(all_results)
    stability_table = stability_frame(data, all_selections)
    RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    result_table.to_csv(RESULTS_PATH, index=False)
    stability_table.to_csv(STABILITY_PATH, index=False)
    write_report(data, result_table, stability_table)

    print(f"Saved {REPORT_PATH}")
    print(f"Saved {RESULTS_PATH}")
    print(f"Saved {STABILITY_PATH}")


if __name__ == "__main__":
    main()
