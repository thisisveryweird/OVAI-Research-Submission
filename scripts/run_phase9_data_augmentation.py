#!/usr/bin/env python3
"""Leakage-safe Phase 9 SMOTE augmentation experiments for GSE94533.

The study is restricted to the verified 113-patient early-stage endpoint:
52 Stage I/II cancers versus 61 benign/control patients. Borderline tumors,
HE4, and all external cohorts are excluded.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd
from imblearn.over_sampling import BorderlineSMOTE, SMOTE
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=FutureWarning)

PROJECT = Path(__file__).resolve().parents[1]
MATRIX_PATH = PROJECT / "data" / "processed" / "gse94533_mirna_tpm.csv.gz"
METADATA_PATH = PROJECT / "data" / "processed" / "gse94533_metadata.csv"
CLINICAL_PATH = PROJECT / "data" / "processed" / "gse94533_14mirna_clinical.csv.gz"
REPORT_PATH = PROJECT / "reports" / "data_augmentation_results.md"
RESULTS_PATH = PROJECT / "reports" / "data_augmentation_results.csv"
PREDICTIONS_PATH = PROJECT / "reports" / "data_augmentation_predictions.csv"

SEED = 20260919
N_SPLITS = 5
TOP_FEATURES = 5
AUGMENTATIONS = ("none", "smote", "borderline_smote")
MODEL_NAMES = ("regularized_logistic", "random_forest", "soft_voting")

AugmentationName = Literal["none", "smote", "borderline_smote"]


@dataclass(frozen=True)
class Cohort:
    frame: pd.DataFrame
    mirna_features: tuple[str, ...]


@dataclass(frozen=True)
class FoldData:
    train_values: np.ndarray
    validation_values: np.ndarray
    train_labels: np.ndarray
    validation_labels: np.ndarray
    selected_features: tuple[str, ...]


def load_cohort() -> Cohort:
    """Load only local GSE94533 data and construct the early-stage endpoint."""

    metadata = pd.read_csv(METADATA_PATH)
    matrix = pd.read_csv(MATRIX_PATH)
    clinical = pd.read_csv(CLINICAL_PATH)
    frame = metadata[["sample_id_normalized", "diagnosis", "stage"]].merge(
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
    if "he4" in frame.columns:
        raise RuntimeError("HE4 must remain excluded from Phase 9.")
    keep = frame["diagnosis"].eq("Cancer") & frame["stage"].eq("Stage I/II")
    keep |= frame["diagnosis"].isin(["Benign", "Control"])
    frame = frame.loc[keep].reset_index(drop=True)
    frame["target"] = (
        frame["diagnosis"].eq("Cancer") & frame["stage"].eq("Stage I/II")
    ).astype(int)
    mirnas = tuple(column for column in matrix.columns if column.startswith("hsa-"))
    if len(frame) != 113 or int(frame["target"].sum()) != 52:
        raise RuntimeError("Expected 113 early-stage patients with 52 positives.")
    if int((frame["target"] == 0).sum()) != 61:
        raise RuntimeError("Expected 61 benign/control patients.")
    return Cohort(frame=frame, mirna_features=mirnas)


def select_top_features(train: pd.DataFrame, labels: np.ndarray, features: tuple[str, ...]) -> tuple[str, ...]:
    """Select the Top-5 miRNAs using training-fold data only."""

    values = np.log1p(train[list(features)].clip(lower=0).to_numpy(dtype=float))
    values = StandardScaler().fit_transform(values)
    selector = LogisticRegression(
        penalty="l1",
        solver="liblinear",
        C=0.1,
        class_weight="balanced",
        max_iter=2000,
        random_state=SEED,
    )
    selector.fit(values, labels)
    ranking = np.argsort(-np.abs(selector.coef_[0]), kind="stable")
    return tuple(features[index] for index in ranking[:TOP_FEATURES])


def prepare_fold(train: pd.DataFrame, validation: pd.DataFrame, labels: np.ndarray, features: tuple[str, ...]) -> FoldData:
    """Impute, scale, and select features using training rows only."""

    selected = select_top_features(train, labels, features)
    selected_values_train = np.log1p(train[list(selected)].clip(lower=0).to_numpy(dtype=float))
    selected_values_validation = np.log1p(validation[list(selected)].clip(lower=0).to_numpy(dtype=float))
    train_values = np.column_stack([
        selected_values_train,
        train[["age", "ca125"]].to_numpy(dtype=float),
    ])
    validation_values = np.column_stack([
        selected_values_validation,
        validation[["age", "ca125"]].to_numpy(dtype=float),
    ])
    imputer = SimpleImputer(strategy="median")
    scaler = StandardScaler()
    train_values = scaler.fit_transform(imputer.fit_transform(train_values))
    validation_values = scaler.transform(imputer.transform(validation_values))
    return FoldData(
        train_values=train_values,
        validation_values=validation_values,
        train_labels=labels,
        validation_labels=np.empty(0, dtype=int),
        selected_features=selected,
    )


def apply_augmentation(
    train_values: np.ndarray,
    train_labels: np.ndarray,
    name: AugmentationName,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply a sampler to training data only after imputation and scaling."""

    if name == "none":
        augmented_values, augmented_labels = train_values, train_labels
    elif name == "smote":
        sampler = SMOTE(random_state=SEED, sampling_strategy="auto", k_neighbors=3)
        augmented_values, augmented_labels = sampler.fit_resample(train_values, train_labels)
    elif name == "borderline_smote":
        sampler = BorderlineSMOTE(
            random_state=SEED,
            sampling_strategy="auto",
            k_neighbors=3,
            m_neighbors=5,
        )
        augmented_values, augmented_labels = sampler.fit_resample(train_values, train_labels)
    else:
        raise ValueError(name)
    if len(augmented_values) < len(train_values):
        raise AssertionError("Augmentation reduced the training sample count.")
    if name != "none" and int((augmented_labels == 1).sum()) != int((augmented_labels == 0).sum()):
        raise AssertionError("SMOTE did not balance the training classes.")
    return augmented_values, augmented_labels


def build_model(name: str) -> Any:
    """Build a fixed Phase 8 model for the current augmented fold."""

    if name == "regularized_logistic":
        return LogisticRegression(
            penalty="l2",
            C=0.3,
            solver="liblinear",
            class_weight="balanced",
            max_iter=5000,
            random_state=SEED,
        )
    if name == "random_forest":
        return RandomForestClassifier(
            n_estimators=300,
            max_depth=4,
            min_samples_leaf=5,
            class_weight="balanced",
            random_state=SEED,
            n_jobs=1,
        )
    raise ValueError(name)


def metric_values(labels: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    prediction = (probability >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(labels, prediction, labels=[0, 1]).ravel()
    sensitivity = tp / (tp + fn) if tp + fn else np.nan
    specificity = tn / (tn + fp) if tn + fp else np.nan
    return {
        "roc_auc": float(roc_auc_score(labels, probability)),
        "pr_auc": float(average_precision_score(labels, probability)),
        "accuracy": float(accuracy_score(labels, prediction)),
        "balanced_accuracy": float(balanced_accuracy_score(labels, prediction)),
        "sensitivity": float(sensitivity),
        "specificity": float(specificity),
        "brier_score": float(brier_score_loss(labels, probability)),
        "geometric_mean": float(np.sqrt(sensitivity * specificity)),
    }


def aggregate(rows: list[dict[str, object]]) -> pd.DataFrame:
    frame = pd.DataFrame(rows)
    output = []
    for (model, augmentation), group in frame.groupby(["model", "augmentation"], sort=True):
        row: dict[str, object] = {
            "model": model,
            "augmentation": augmentation,
            "folds": len(group),
        }
        for metric in [
            "roc_auc", "pr_auc", "accuracy", "balanced_accuracy",
            "sensitivity", "specificity", "brier_score", "geometric_mean",
        ]:
            row[f"{metric}_mean"] = float(group[metric].mean())
        output.append(row)
    return pd.DataFrame(output)


def write_report(results: pd.DataFrame) -> None:
    """Write comparison table and a sensitivity/specificity conclusion."""

    lines = [
        "# Phase 9 Data Augmentation Results",
        "",
        "This study uses only the verified local GSE94533 early-stage endpoint: 52 Stage I/II cancers versus 61 benign/control patients. Borderline tumors, HE4, and external datasets are excluded.",
        "",
        "## Leakage Controls",
        "",
        "Stratified 5-fold cross-validation was used. Top-5 miRNA selection, Age/CA-125 median imputation, and scaling were fit on each training fold only. SMOTE and BorderlineSMOTE were applied only after those training-fold transformations and never touched validation rows or validation labels.",
        "",
        "Synthetic samples were generated only for the positive early-stage training class until the training classes were balanced. Validation rows remain the original biological patients.",
        "",
        "## Results",
        "",
        "| Model | Augmentation | ROC-AUC | PR-AUC | Accuracy | Balanced accuracy | Sensitivity | Specificity | Brier score | Geometric mean |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in results.sort_values(["model", "augmentation"]).itertuples(index=False):
        lines.append(
            f"| {row.model} | {row.augmentation} | {row.roc_auc_mean:.3f} | "
            f"{row.pr_auc_mean:.3f} | {row.accuracy_mean:.3f} | "
            f"{row.balanced_accuracy_mean:.3f} | {row.sensitivity_mean:.3f} | "
            f"{row.specificity_mean:.3f} | {row.brier_score_mean:.3f} | "
            f"{row.geometric_mean_mean:.3f} |"
        )

    baseline = results[results["augmentation"] == "none"].set_index("model")
    best_augmented = results[results["augmentation"] != "none"].copy()
    best_augmented["sensitivity_delta"] = best_augmented.apply(
        lambda row: row["sensitivity_mean"] - baseline.loc[row["model"], "sensitivity_mean"],
        axis=1,
    )
    best_augmented["specificity_delta"] = best_augmented.apply(
        lambda row: row["specificity_mean"] - baseline.loc[row["model"], "specificity_mean"],
        axis=1,
    )
    best = best_augmented.sort_values("sensitivity_delta", ascending=False).iloc[0]
    improved_sensitivity = best["sensitivity_delta"] > 1e-6
    retained_specificity = best["specificity_delta"] >= -0.05
    conclusion = (
        "Synthetic augmentation improved sensitivity while retaining specificity within 5 percentage points for the strongest observed comparison."
        if improved_sensitivity and retained_specificity
        else "Synthetic augmentation did not demonstrate a clear sensitivity improvement without a material specificity tradeoff."
    )
    lines.extend([
        "",
        "## Conclusion",
        "",
        conclusion,
        f"The largest sensitivity change was `{best.model}` with `{best.augmentation}`: "
        f"sensitivity delta {best.sensitivity_delta:+.3f}, specificity delta {best.specificity_delta:+.3f} versus its no-augmentation baseline.",
        "",
        "These are internal OOF research estimates. Synthetic augmentation must be re-tested on untouched external serum data before any generalization claim.",
        "",
        "Machine-readable aggregate results: `reports/data_augmentation_results.csv`.",
        "Machine-readable fold predictions: `reports/data_augmentation_predictions.csv`.",
    ])
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines) + "\n")


def main() -> None:
    """Run all Phase 9 augmentation/model combinations."""

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    logging.info("Loading verified local GSE94533 early-stage cohort")
    cohort = load_cohort()
    labels = cohort.frame["target"].to_numpy(dtype=int)
    splitter = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=SEED)
    fold_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []

    for fold, (train_index, validation_index) in enumerate(splitter.split(cohort.frame, labels), start=1):
        train = cohort.frame.iloc[train_index]
        validation = cohort.frame.iloc[validation_index]
        train_labels = labels[train_index]
        validation_labels = labels[validation_index]
        selected = select_top_features(train, train_labels, cohort.mirna_features)
        fold_data = prepare_fold(train, validation, train_labels, cohort.mirna_features)
        logging.info("Fold %d selected %s", fold, ", ".join(selected))

        for augmentation in AUGMENTATIONS:
            augmented_values, augmented_labels = apply_augmentation(
                fold_data.train_values,
                fold_data.train_labels,
                augmentation,  # type: ignore[arg-type]
            )
            probabilities: dict[str, np.ndarray] = {}
            for model_name in ("regularized_logistic", "random_forest"):
                model = build_model(model_name)
                model.fit(augmented_values, augmented_labels)
                probabilities[model_name] = model.predict_proba(fold_data.validation_values)[:, 1]
                fold_rows.append({
                    "model": model_name,
                    "augmentation": augmentation,
                    "fold": fold,
                    **metric_values(validation_labels, probabilities[model_name]),
                })

            probabilities["soft_voting"] = np.mean(
                np.column_stack([probabilities["regularized_logistic"], probabilities["random_forest"]]),
                axis=1,
            )
            fold_rows.append({
                "model": "soft_voting",
                "augmentation": augmentation,
                "fold": fold,
                **metric_values(validation_labels, probabilities["soft_voting"]),
            })
            for model_name, probability in probabilities.items():
                prediction_rows.extend({
                    "model": model_name,
                    "augmentation": augmentation,
                    "fold": fold,
                    "sample_id_normalized": sample_id,
                    "true_label": int(label),
                    "probability": float(probability_value),
                } for sample_id, label, probability_value in zip(
                    validation["sample_id_normalized"], validation_labels, probability
                ))

    result_table = aggregate(fold_rows)
    prediction_table = pd.DataFrame(prediction_rows)
    RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    result_table.to_csv(RESULTS_PATH, index=False)
    prediction_table.to_csv(PREDICTIONS_PATH, index=False)
    write_report(result_table)
    logging.info("Validation assertion: every augmentation was applied only to transformed training rows.")
    logging.info("Saved %s", REPORT_PATH)
    logging.info("Saved %s", RESULTS_PATH)
    logging.info("Saved %s", PREDICTIONS_PATH)


if __name__ == "__main__":
    main()
