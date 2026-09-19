#!/usr/bin/env python3
"""Leakage-safe Phase 8 model and ensemble experiments for GSE94533.

Only the verified local GSE94533 cohort is loaded. The endpoint is Stage I/II
cancer versus benign/control; borderline tumors are excluded. All preprocessing,
Top-5 selection, base-model tuning, calibration, and stacking meta-training are
contained within the outer training fold.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.calibration import CalibratedClassifierCV
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
from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from xgboost import XGBClassifier

warnings.filterwarnings("ignore", category=FutureWarning)

PROJECT = Path(__file__).resolve().parents[1]
MATRIX_PATH = PROJECT / "data" / "processed" / "gse94533_mirna_tpm.csv.gz"
METADATA_PATH = PROJECT / "data" / "processed" / "gse94533_metadata.csv"
CLINICAL_PATH = PROJECT / "data" / "processed" / "gse94533_14mirna_clinical.csv.gz"
REPORT_PATH = PROJECT / "reports" / "model_ensemble_comparison.md"
RESULTS_PATH = PROJECT / "reports" / "model_ensemble_results.csv"
PREDICTIONS_PATH = PROJECT / "reports" / "model_ensemble_predictions.csv"

SEED = 20260919
OUTER_SPLITS = 5
INNER_SPLITS = 3
TOP_FEATURES = 5
MODEL_NAMES = ("regularized_logistic", "calibrated_svm", "random_forest", "xgboost")
ENSEMBLE_NAMES = ("soft_voting", "stacked_logistic")
ALL_NAMES = MODEL_NAMES + ENSEMBLE_NAMES

ModelName = Literal[
    "regularized_logistic",
    "calibrated_svm",
    "random_forest",
    "xgboost",
]


@dataclass(frozen=True)
class Cohort:
    frame: pd.DataFrame
    mirna_features: tuple[str, ...]


def load_cohort() -> Cohort:
    """Load only local GSE94533 data and construct the early-stage endpoint."""

    metadata = pd.read_csv(METADATA_PATH)
    matrix = pd.read_csv(MATRIX_PATH)
    clinical = pd.read_csv(CLINICAL_PATH)
    frame = metadata[
        ["sample_id_normalized", "diagnosis", "stage"]
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
    if "he4" in frame.columns:
        raise RuntimeError("HE4 must remain excluded from Phase 8.")
    keep = frame["diagnosis"].eq("Cancer") & frame["stage"].eq("Stage I/II")
    keep |= frame["diagnosis"].isin(["Benign", "Control"])
    frame = frame.loc[keep].reset_index(drop=True)
    labels = (
        frame["diagnosis"].eq("Cancer")
        & frame["stage"].eq("Stage I/II")
    ).astype(int)
    frame["target"] = labels
    mirnas = tuple(column for column in matrix.columns if column.startswith("hsa-"))
    if labels.sum() != 52 or (labels == 0).sum() != 61:
        raise RuntimeError("Unexpected early-stage endpoint composition.")
    return Cohort(frame=frame, mirna_features=mirnas)


def select_top_features(train: pd.DataFrame, labels: np.ndarray, features: tuple[str, ...]) -> tuple[str, ...]:
    """Fit the Phase 7 L1 selector on one training split only."""

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


def transform_fold(
    train: pd.DataFrame,
    validation: pd.DataFrame,
    selected: tuple[str, ...],
) -> tuple[np.ndarray, np.ndarray]:
    """Impute age/CA-125 and scale the selected multimodal features fold-locally."""

    columns = list(selected) + ["age", "ca125"]
    imputer = SimpleImputer(strategy="median")
    scaler = StandardScaler()
    train_values = np.log1p(train[list(selected)].clip(lower=0).to_numpy(dtype=float))
    validation_values = np.log1p(validation[list(selected)].clip(lower=0).to_numpy(dtype=float))
    train_values = np.column_stack([train_values, train[["age", "ca125"]].to_numpy(dtype=float)])
    validation_values = np.column_stack([validation_values, validation[["age", "ca125"]].to_numpy(dtype=float)])
    train_values = imputer.fit_transform(train_values)
    validation_values = imputer.transform(validation_values)
    return scaler.fit_transform(train_values), scaler.transform(validation_values)


def model_spec(name: ModelName) -> tuple[Any, dict[str, list[Any]]]:
    """Return an estimator and small prespecified inner-CV grid."""

    if name == "regularized_logistic":
        return LogisticRegression(
            penalty="l2", solver="liblinear", class_weight="balanced", max_iter=5000,
            random_state=SEED,
        ), {"C": [0.03, 0.1, 0.3, 1.0]}
    if name == "calibrated_svm":
        return CalibratedClassifierCV(
            estimator=LinearSVC(
                class_weight="balanced", dual="auto", max_iter=10000, random_state=SEED,
            ),
            method="sigmoid",
            cv=3,
        ), {"estimator__C": [0.03, 0.1, 0.3, 1.0]}
    if name == "random_forest":
        return RandomForestClassifier(
            n_estimators=300, class_weight="balanced", random_state=SEED, n_jobs=1,
        ), {
            "max_depth": [2, 4, None],
            "min_samples_leaf": [2, 5, 10],
        }
    if name == "xgboost":
        return XGBClassifier(
            n_estimators=200,
            learning_rate=0.05,
            subsample=0.8,
            colsample_bytree=0.8,
            objective="binary:logistic",
            eval_metric="logloss",
            random_state=SEED,
            n_jobs=1,
        ), {
            "max_depth": [1, 2, 3],
            "min_child_weight": [1, 5],
        }
    raise ValueError(name)


def fit_tuned_model(name: ModelName, train_values: np.ndarray, labels: np.ndarray) -> Any:
    """Tune one base model inside the current outer training fold."""

    estimator, grid = model_spec(name)
    search = GridSearchCV(
        estimator,
        grid,
        scoring="roc_auc",
        cv=StratifiedKFold(n_splits=INNER_SPLITS, shuffle=True, random_state=SEED),
        n_jobs=1,
        refit=True,
        error_score="raise",
    )
    search.fit(train_values, labels)
    return search.best_estimator_


def metrics(labels: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    prediction = (probability >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(labels, prediction, labels=[0, 1]).ravel()
    return {
        "roc_auc": float(roc_auc_score(labels, probability)),
        "pr_auc": float(average_precision_score(labels, probability)),
        "accuracy": float(accuracy_score(labels, prediction)),
        "balanced_accuracy": float(balanced_accuracy_score(labels, prediction)),
        "sensitivity": float(tp / (tp + fn)) if tp + fn else np.nan,
        "specificity": float(tn / (tn + fp)) if tn + fp else np.nan,
        "brier_score": float(brier_score_loss(labels, probability)),
    }


def inner_oof_base_predictions(
    name: ModelName,
    train_values: np.ndarray,
    labels: np.ndarray,
) -> np.ndarray:
    """Generate strictly OOF base predictions for the stacking meta-model."""

    inner = StratifiedKFold(n_splits=INNER_SPLITS, shuffle=True, random_state=SEED + 1)
    out = np.full(len(labels), np.nan)
    for inner_train, inner_validation in inner.split(train_values, labels):
        model = fit_tuned_model(name, train_values[inner_train], labels[inner_train])
        out[inner_validation] = model.predict_proba(train_values[inner_validation])[:, 1]
    if not np.isfinite(out).all():
        raise RuntimeError("Stacking base predictions are not fully OOF.")
    return out


def aggregate(rows: list[dict[str, object]]) -> pd.DataFrame:
    frame = pd.DataFrame(rows)
    grouping = ["model"]
    output = []
    for model, group in frame.groupby(grouping, sort=True):
        model_name = model[0] if isinstance(model, tuple) else model
        row = {"model": model_name, "folds": len(group)}
        for metric_name in ["roc_auc", "pr_auc", "accuracy", "balanced_accuracy", "sensitivity", "specificity", "brier_score"]:
            row[f"{metric_name}_mean"] = group[metric_name].mean()
        output.append(row)
    return pd.DataFrame(output)


def write_report(results: pd.DataFrame) -> None:
    lines = [
        "# Phase 8 Model and Ensemble Comparison",
        "",
        "This study uses only verified local GSE94533 data. The endpoint is Stage I/II cancer versus benign/control; borderline tumors are excluded. HE4 and all external datasets are excluded.",
        "",
        "## Leakage Controls",
        "",
        "Stratified 5-fold outer evaluation was used. In every outer training fold, the Top-5 miRNA selector, median imputer, scaler, and inner-CV hyperparameter search were fit on training data only. The stacked logistic meta-model was trained only on inner out-of-fold base-model probabilities; no in-sample base predictions were used.",
        "",
        "## Results",
        "",
        "| Model | ROC-AUC | PR-AUC | Accuracy | Balanced accuracy | Sensitivity | Specificity | Brier score |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in results.sort_values("roc_auc_mean", ascending=False).itertuples(index=False):
        lines.append(
            f"| {row.model} | {row.roc_auc_mean:.3f} | {row.pr_auc_mean:.3f} | "
            f"{row.accuracy_mean:.3f} | {row.balanced_accuracy_mean:.3f} | "
            f"{row.sensitivity_mean:.3f} | {row.specificity_mean:.3f} | {row.brier_score_mean:.3f} |"
        )
    lines.extend([
        "",
        "## Validation Printout",
        "",
        "`stacked_logistic`: meta-training used inner OOF base predictions only; in-sample base predictions were not used.",
        "",
        "These are internal out-of-fold research estimates. No external validation claim is made.",
    ])
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines) + "\n")


def main() -> None:
    cohort = load_cohort()
    outer = StratifiedKFold(n_splits=OUTER_SPLITS, shuffle=True, random_state=SEED)
    labels = cohort.frame["target"].to_numpy(dtype=int)
    fold_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []

    for fold, (train_index, validation_index) in enumerate(outer.split(cohort.frame, labels), start=1):
        train = cohort.frame.iloc[train_index]
        validation = cohort.frame.iloc[validation_index]
        train_labels = labels[train_index]
        validation_labels = labels[validation_index]
        selected = select_top_features(train, train_labels, cohort.mirna_features)
        train_values, validation_values = transform_fold(train, validation, selected)
        base_probabilities: dict[str, np.ndarray] = {}
        for name in MODEL_NAMES:
            model = fit_tuned_model(name, train_values, train_labels)
            probability = model.predict_proba(validation_values)[:, 1]
            base_probabilities[name] = probability
            fold_rows.append({"model": name, "fold": fold, **metrics(validation_labels, probability)})
            prediction_rows.extend({
                "model": name,
                "fold": fold,
                "sample_id_normalized": sample_id,
                "true_label": int(label),
                "probability": float(probability_value),
            } for sample_id, label, probability_value in zip(
                validation["sample_id_normalized"], validation_labels, probability
            ))

        soft_probability = np.mean(np.column_stack(list(base_probabilities.values())), axis=1)
        fold_rows.append({"model": "soft_voting", "fold": fold, **metrics(validation_labels, soft_probability)})
        prediction_rows.extend({
            "model": "soft_voting", "fold": fold,
            "sample_id_normalized": sample_id,
            "true_label": int(label), "probability": float(probability_value),
        } for sample_id, label, probability_value in zip(
            validation["sample_id_normalized"], validation_labels, soft_probability
        ))

        meta_train_columns = []
        for name in MODEL_NAMES:
            meta_train_columns.append(inner_oof_base_predictions(name, train_values, train_labels))
        meta_train = np.column_stack(meta_train_columns)
        meta_model = LogisticRegression(C=0.3, class_weight="balanced", solver="liblinear", max_iter=5000, random_state=SEED)
        meta_model.fit(meta_train, train_labels)
        meta_probability = meta_model.predict_proba(np.column_stack(list(base_probabilities.values())))[:, 1]
        fold_rows.append({"model": "stacked_logistic", "fold": fold, **metrics(validation_labels, meta_probability)})
        prediction_rows.extend({
            "model": "stacked_logistic", "fold": fold,
            "sample_id_normalized": sample_id,
            "true_label": int(label), "probability": float(probability_value),
        } for sample_id, label, probability_value in zip(
            validation["sample_id_normalized"], validation_labels, meta_probability
        ))

    result_table = aggregate(fold_rows)
    result_table.to_csv(RESULTS_PATH, index=False)
    pd.DataFrame(prediction_rows).to_csv(PREDICTIONS_PATH, index=False)
    write_report(result_table)
    print("Validation: stacked_logistic meta-training used inner OOF base predictions only.")
    print(f"Saved {REPORT_PATH}")
    print(f"Saved {RESULTS_PATH}")
    print(f"Saved {PREDICTIONS_PATH}")


if __name__ == "__main__":
    main()
