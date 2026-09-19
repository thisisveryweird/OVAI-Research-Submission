#!/usr/bin/env python3
"""Leakage-safe Phase 7 clinical ablation study for GSE94533.

Only the verified local GSE94533 cohort is loaded. HE4 is intentionally absent
because it is not present in the verified clinical table. All feature selection,
clinical imputation, scaling, calibration, and model fitting occur inside each
training fold.
"""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.impute import KNNImputer, SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC

warnings.filterwarnings("ignore", category=FutureWarning)

PROJECT = Path(__file__).resolve().parents[1]
MATRIX_PATH = PROJECT / "data" / "processed" / "gse94533_mirna_tpm.csv.gz"
METADATA_PATH = PROJECT / "data" / "processed" / "gse94533_metadata.csv"
CLINICAL_PATH = PROJECT / "data" / "processed" / "gse94533_14mirna_clinical.csv.gz"
REPORT_PATH = PROJECT / "reports" / "clinical_ablation_results.md"
RESULTS_PATH = PROJECT / "reports" / "clinical_ablation_results.csv"

SEED = 20260919
N_SPLITS = 5
TOP_FEATURES = 5
MIRNA_PREFIX = "hsa-"
CLINICAL_FEATURES = ("age", "ca125")
IMPUTATION_MODES = ("median", "knn_ca125")
MODELS = ("regularized_logistic", "calibrated_svm")
PIPELINES = ("mirna_only", "clinical_only", "combined")

ModelName = Literal["regularized_logistic", "calibrated_svm"]
ImputationMode = Literal["median", "knn_ca125"]


@dataclass(frozen=True)
class Cohort:
    frame: pd.DataFrame
    mirna_features: tuple[str, ...]


@dataclass(frozen=True)
class Endpoint:
    name: str
    mask_description: str


def load_cohort() -> Cohort:
    """Load and validate only the local 179-patient GSE94533 cohort."""

    metadata = pd.read_csv(METADATA_PATH)
    matrix = pd.read_csv(MATRIX_PATH)
    clinical = pd.read_csv(CLINICAL_PATH)
    if len(metadata) != 179 or len(matrix) != 179:
        raise RuntimeError("Phase 7 requires the verified 179-patient GSE94533 cohort.")

    frame = metadata[
        ["sample_id_normalized", "diagnosis", "binary_label", "stage"]
    ].merge(
        matrix,
        on="sample_id_normalized",
        how="inner",
        validate="one_to_one",
    ).merge(
        clinical[["sample_id_normalized", "age", "ca125"]],
        on="sample_id_normalized",
        how="left",
        validate="one_to_one",
    )
    mirnas = tuple(column for column in matrix.columns if column.startswith(MIRNA_PREFIX))
    if len(mirnas) < TOP_FEATURES:
        raise RuntimeError("The local matrix has fewer than five miRNA features.")

    expected = {"Cancer": 96, "Benign": 46, "Borderline": 22, "Control": 15}
    observed = frame["diagnosis"].value_counts().to_dict()
    if observed != expected:
        raise RuntimeError(f"Unexpected diagnosis composition: {observed}")
    if "he4" in frame.columns:
        raise RuntimeError("HE4 must remain explicitly excluded from Phase 7.")
    return Cohort(frame=frame, mirna_features=mirnas)


def endpoint_data(cohort: Cohort, endpoint: Endpoint) -> tuple[pd.DataFrame, np.ndarray]:
    """Construct endpoints without relabeling borderline tumors."""

    frame = cohort.frame
    if endpoint.name == "all_cancer_vs_non_cancer":
        subset = frame.copy()
        labels = subset["binary_label"].astype(int).to_numpy()
    elif endpoint.name == "early_stage_vs_benign_control":
        keep = frame["diagnosis"].eq("Cancer") & frame["stage"].eq("Stage I/II")
        keep |= frame["diagnosis"].isin(["Benign", "Control"])
        subset = frame.loc[keep].copy()
        labels = (
            subset["diagnosis"].eq("Cancer")
            & subset["stage"].eq("Stage I/II")
        ).astype(int).to_numpy()
    else:
        raise ValueError(endpoint.name)
    return subset, labels


def impute_clinical(
    train: pd.DataFrame,
    validation: pd.DataFrame,
    mode: ImputationMode,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit clinical imputation on training rows only.

    Age always uses fold-local median imputation. The median mode uses a
    separate fold-local median for CA-125. The KNN mode uses KNNImputer over
    the two clinical columns; because age is complete for all but one row and
    is included as the observed covariate, CA-125 imputation is informed by
    the available age values while remaining training-fold-only.
    """

    age_imputer = SimpleImputer(strategy="median")
    train_age = age_imputer.fit_transform(train[["age"]])
    validation_age = age_imputer.transform(validation[["age"]])

    if mode == "median":
        ca125_imputer = SimpleImputer(strategy="median")
        train_ca125 = ca125_imputer.fit_transform(train[["ca125"]])
        validation_ca125 = ca125_imputer.transform(validation[["ca125"]])
    elif mode == "knn_ca125":
        knn = KNNImputer(n_neighbors=5, weights="distance")
        train_clinical = knn.fit_transform(train[["age", "ca125"]])
        validation_clinical = knn.transform(validation[["age", "ca125"]])
        train_ca125 = train_clinical[:, [1]]
        validation_ca125 = validation_clinical[:, [1]]
    else:
        raise ValueError(mode)

    return (
        np.column_stack([train_age[:, 0], train_ca125[:, 0]]),
        np.column_stack([validation_age[:, 0], validation_ca125[:, 0]]),
    )


def select_top_mirnas(
    train: pd.DataFrame,
    labels: np.ndarray,
    mirna_features: tuple[str, ...],
) -> tuple[str, ...]:
    """Select five miRNAs using only one training fold.

    The sparse selector is deliberately fit only to miRNA values. Clinical
    variables cannot influence which molecular features enter the ablation.
    """

    train_values = np.log1p(train[list(mirna_features)].clip(lower=0).to_numpy(dtype=float))
    scaler = StandardScaler()
    scaled = scaler.fit_transform(train_values)
    selector = LogisticRegression(
        penalty="l1",
        solver="liblinear",
        C=0.1,
        class_weight="balanced",
        max_iter=2000,
        random_state=SEED,
    )
    selector.fit(scaled, labels)
    ranking = np.argsort(-np.abs(selector.coef_[0]), kind="stable")
    return tuple(mirna_features[index] for index in ranking[:TOP_FEATURES])


def fit_model(name: ModelName) -> object:
    """Create a classifier whose calibration is fit inside the training fold."""

    if name == "regularized_logistic":
        return LogisticRegression(
            penalty="l2",
            C=0.3,
            solver="liblinear",
            class_weight="balanced",
            max_iter=5000,
            random_state=SEED,
        )
    if name == "calibrated_svm":
        return CalibratedClassifierCV(
            estimator=LinearSVC(
                C=0.3,
                class_weight="balanced",
                dual="auto",
                max_iter=10000,
                random_state=SEED,
            ),
            method="sigmoid",
            cv=3,
        )
    raise ValueError(name)


def metrics(labels: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    prediction = (probability >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(labels, prediction, labels=[0, 1]).ravel()
    return {
        "roc_auc": float(roc_auc_score(labels, probability)),
        "balanced_accuracy": float(balanced_accuracy_score(labels, prediction)),
        "sensitivity": float(tp / (tp + fn)) if tp + fn else np.nan,
        "specificity": float(tn / (tn + fp)) if tn + fp else np.nan,
        "brier_score": float(brier_score_loss(labels, probability)),
        "accuracy": float(accuracy_score(labels, prediction)),
    }


def evaluate_configuration(
    cohort: Cohort,
    endpoint: Endpoint,
    frame: pd.DataFrame,
    labels: np.ndarray,
    pipeline: str,
    model_name: ModelName,
    imputation_mode: ImputationMode,
    fold_selected_features: dict[int, tuple[str, ...]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Evaluate one ablation/model/imputation configuration."""

    splitter = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED)
    fold_rows: list[dict[str, object]] = []
    feature_rows: list[dict[str, object]] = []
    selected_counts: dict[str, int] = {feature: 0 for feature in cohort.mirna_features}

    for fold, (train_index, validation_index) in enumerate(splitter.split(frame, labels), start=1):
        train = frame.iloc[train_index]
        validation = frame.iloc[validation_index]
        train_labels = labels[train_index]
        validation_labels = labels[validation_index]

        selected = fold_selected_features[fold]
        for feature in selected:
            selected_counts[feature] += 1

        train_clinical, validation_clinical = impute_clinical(train, validation, imputation_mode)
        train_mirna = np.log1p(train[list(selected)].clip(lower=0).to_numpy(dtype=float))
        validation_mirna = np.log1p(validation[list(selected)].clip(lower=0).to_numpy(dtype=float))

        if pipeline == "mirna_only":
            train_values = train_mirna
            validation_values = validation_mirna
        elif pipeline == "clinical_only":
            train_values = train_clinical
            validation_values = validation_clinical
        elif pipeline == "combined":
            train_values = np.column_stack([train_mirna, train_clinical])
            validation_values = np.column_stack([validation_mirna, validation_clinical])
        else:
            raise ValueError(pipeline)

        scaler = StandardScaler()
        train_values = scaler.fit_transform(train_values)
        validation_values = scaler.transform(validation_values)
        model = fit_model(model_name)
        model.fit(train_values, train_labels)
        probability = model.predict_proba(validation_values)[:, 1]
        fold_metrics = metrics(validation_labels, probability)
        fold_rows.append({
            "endpoint": endpoint.name,
            "pipeline": pipeline,
            "model": model_name,
            "imputation": imputation_mode,
            "fold": fold,
            "n_train": len(train_index),
            "n_validation": len(validation_index),
            "selected_mirnas": ";".join(selected),
            **fold_metrics,
        })

    for feature, count in selected_counts.items():
        feature_rows.append({
            "endpoint": endpoint.name,
            "pipeline": pipeline,
            "model": model_name,
            "imputation": imputation_mode,
            "feature": feature,
            "selected_folds": count,
            "selection_frequency": count / N_SPLITS,
        })
    return fold_rows, feature_rows


def aggregate_results(fold_rows: list[dict[str, object]]) -> pd.DataFrame:
    """Aggregate fold metrics while retaining fold-level reproducibility."""

    frame = pd.DataFrame(fold_rows)
    grouping = ["endpoint", "pipeline", "model", "imputation"]
    rows = []
    for keys, group in frame.groupby(grouping, sort=True):
        row = dict(zip(grouping, keys))
        row.update({
            "folds": len(group),
            "roc_auc_mean": group["roc_auc"].mean(),
            "balanced_accuracy_mean": group["balanced_accuracy"].mean(),
            "sensitivity_mean": group["sensitivity"].mean(),
            "specificity_mean": group["specificity"].mean(),
            "brier_score_mean": group["brier_score"].mean(),
            "accuracy_mean": group["accuracy"].mean(),
        })
        rows.append(row)
    return pd.DataFrame(rows)


def write_report(results: pd.DataFrame, stability: pd.DataFrame) -> None:
    """Write the requested ablation comparison and stability summary."""

    lines = [
        "# Phase 7 Clinical Ablation Results",
        "",
        "This study uses only the verified local GSE94533 cohort. HE4 is not present and is excluded from every pipeline. No external dataset is loaded.",
        "",
        "## Protocol",
        "",
        "Each configuration uses shuffled stratified 5-fold cross-validation. Age is median-imputed inside each training fold. CA-125 is evaluated with either fold-local median imputation or fold-local KNN imputation. The Top-5 miRNA signature is selected by fold-local L1 logistic regression from the full local miRNA feature space; no global five-feature list is used.",
        "",
        "Borderline tumors are included in the all-cancer versus non-cancer endpoint and excluded from the Stage I/II versus benign/control endpoint.",
        "",
        "## Results",
        "",
        "| Endpoint | Pipeline | Model | Imputation | ROC-AUC | Balanced accuracy | Sensitivity | Specificity | Brier score |",
        "|---|---|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in results.sort_values(["endpoint", "pipeline", "model", "imputation"]).itertuples(index=False):
        lines.append(
            f"| {row.endpoint} | {row.pipeline} | {row.model} | {row.imputation} | "
            f"{row.roc_auc_mean:.3f} | {row.balanced_accuracy_mean:.3f} | "
            f"{row.sensitivity_mean:.3f} | {row.specificity_mean:.3f} | {row.brier_score_mean:.3f} |"
        )

    lines.extend([
        "",
        "## Top-5 Selection Stability",
        "",
        "Selection frequency is measured across the five folds for each configuration. It is a within-cohort stability measure, not evidence of biological causality or external generalization.",
        "",
    ])
    for endpoint in sorted(stability["endpoint"].unique()):
        lines.append(f"### {endpoint}")
        lines.append("")
        for pipeline in PIPELINES:
            for model in MODELS:
                for imputation in IMPUTATION_MODES:
                    subset = stability[
                        (stability["endpoint"] == endpoint)
                        & (stability["pipeline"] == pipeline)
                        & (stability["model"] == model)
                        & (stability["imputation"] == imputation)
                    ].sort_values(["selection_frequency", "feature"], ascending=[False, True]).head(10)
                    top = "; ".join(
                        f"{row.feature} ({int(row.selected_folds)}/5)"
                        for row in subset.itertuples(index=False)
                        if row.selected_folds > 0
                    )
                    lines.append(f"- {pipeline}/{model}/{imputation}: {top or 'No selected features.'}")
        lines.append("")

    lines.extend([
        "## Interpretation Guardrails",
        "",
        "- A clinical-only result is not a molecular diagnostic result; it is an ablation comparator.",
        "- Any combined-model improvement is an internal OOF estimate and requires independent validation before being considered robust.",
        "- CA-125 has 60 missing values of 179; KNN imputation is therefore a sensitivity analysis, not a substitute for measured values.",
        "- HE4 is unavailable and was excluded explicitly.",
        "",
        "Machine-readable aggregate results: `reports/clinical_ablation_results.csv`.",
        "Machine-readable Top-5 stability: `reports/clinical_ablation_stability.csv`.",
    ])
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines) + "\n")


def main() -> None:
    """Run all Phase 7 ablations for both endpoints."""

    cohort = load_cohort()
    fold_rows: list[dict[str, object]] = []
    stability_rows: list[dict[str, object]] = []
    endpoints = (
        Endpoint("all_cancer_vs_non_cancer", "96 cancer versus 83 non-cancer"),
        Endpoint("early_stage_vs_benign_control", "52 Stage I/II cancer versus 61 benign/control"),
    )
    for endpoint in endpoints:
        frame, labels = endpoint_data(cohort, endpoint)
        selection_splitter = StratifiedKFold(
            n_splits=N_SPLITS,
            shuffle=True,
            random_state=SEED,
        )
        fold_selected_features: dict[int, tuple[str, ...]] = {}
        for fold, (train_index, _validation_index) in enumerate(
            selection_splitter.split(frame, labels),
            start=1,
        ):
            fold_selected_features[fold] = select_top_mirnas(
                frame.iloc[train_index],
                labels[train_index],
                cohort.mirna_features,
            )
        for pipeline in PIPELINES:
            for model_name in MODELS:
                for imputation in IMPUTATION_MODES:
                    fold_result, feature_result = evaluate_configuration(
                        cohort,
                        endpoint,
                        frame,
                        labels,
                        pipeline,
                        model_name,  # type: ignore[arg-type]
                        imputation,  # type: ignore[arg-type]
                        fold_selected_features,
                    )
                    fold_rows.extend(fold_result)
                    stability_rows.extend(feature_result)

    fold_table = pd.DataFrame(fold_rows)
    result_table = aggregate_results(fold_rows)
    stability_table = pd.DataFrame(stability_rows)
    RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    result_table.to_csv(RESULTS_PATH, index=False)
    stability_table.to_csv(
        PROJECT / "reports" / "clinical_ablation_stability.csv",
        index=False,
    )
    write_report(result_table, stability_table)
    print(f"Saved {REPORT_PATH}")
    print(f"Saved {RESULTS_PATH}")
    print(f"Saved {len(fold_table)} fold-level rows")


if __name__ == "__main__":
    main()
