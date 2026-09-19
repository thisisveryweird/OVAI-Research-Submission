#!/usr/bin/env python3
"""Run the reproducible patient-level GSE94533 baseline comparison.

The outer test folds are never used for preprocessing, feature selection,
hyperparameter selection, probability calibration, or threshold selection.
External cohorts are intentionally not loaded by this script.
"""

from __future__ import annotations

import json
import platform
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import FunctionTransformer, StandardScaler
from sklearn.svm import LinearSVC
from xgboost import XGBClassifier


PROJECT = Path(__file__).resolve().parents[1]
DATA = PROJECT / "data" / "processed" / "gse94533_14mirna_clinical.csv.gz"
REPORTS = PROJECT / "reports"
COMPARISON = REPORTS / "baseline_comparison.csv"
PREDICTIONS = REPORTS / "baseline_predictions.csv"
CALIBRATION = REPORTS / "baseline_calibration.csv"
METADATA = REPORTS / "baseline_run_metadata.json"

SEED = 20260919
OUTER_SPLITS = 5
INNER_SPLITS = 3
BOOTSTRAPS = 300
THRESHOLDS = np.linspace(0.10, 0.90, 81)

MIRNAS = [
    "hsa-miR-23b-3p", "hsa-miR-29a-3p", "hsa-miR-32-5p",
    "hsa-miR-92a-3p", "hsa-miR-150-5p", "hsa-miR-200a-3p",
    "hsa-miR-200c-3p", "hsa-miR-203a", "hsa-miR-320c",
    "hsa-miR-320d", "hsa-miR-335-5p", "hsa-miR-450b-5p",
    "hsa-miR-1246", "hsa-miR-1307-5p",
]
FEATURES = MIRNAS + ["age", "ca125"]
LOG_FEATURES = MIRNAS + ["ca125"]


@dataclass(frozen=True)
class ModelSpec:
    name: str
    estimator: Any
    grid: dict[str, list[Any]]
    scale: bool


def seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)


def safe_log1p(values: Any) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    return np.log1p(np.clip(array, a_min=0.0, a_max=None))


def make_pipeline(model: Any, feature_count: int, scale: bool) -> Pipeline:
    k_values = [k for k in [5, 10, 14, 16] if k <= feature_count]
    k_values.append("all")
    transformer = ColumnTransformer([
        (
            "log_positive_features",
            Pipeline([
                ("log1p", FunctionTransformer(safe_log1p, validate=False)),
                ("impute", SimpleImputer(strategy="median")),
            ]),
            LOG_FEATURES,
        ),
        ("clinical_age", SimpleImputer(strategy="median"), ["age"]),
    ])
    steps: list[tuple[str, Any]] = [("transform", transformer)]
    if scale:
        steps.append(("scale", StandardScaler()))
    steps.extend([
        ("select", SelectKBest(score_func=f_classif)),
        ("model", model),
    ])
    return Pipeline(steps)


def model_specs(feature_count: int) -> list[ModelSpec]:
    logistic = make_pipeline(
        LogisticRegression(max_iter=5000, solver="liblinear"),
        feature_count,
        scale=True,
    )
    regularized = make_pipeline(
        LogisticRegression(max_iter=5000, solver="liblinear", penalty="l1"),
        feature_count,
        scale=True,
    )
    linear_svm = make_pipeline(
        CalibratedClassifierCV(
            estimator=LinearSVC(class_weight="balanced", dual="auto", max_iter=10000),
            method="sigmoid",
            cv=3,
        ),
        feature_count,
        scale=True,
    )
    calibrated_svm = make_pipeline(
        CalibratedClassifierCV(
            estimator=LinearSVC(class_weight="balanced", dual="auto", max_iter=10000),
            method="isotonic",
            cv=3,
        ),
        feature_count,
        scale=True,
    )
    random_forest = make_pipeline(
        RandomForestClassifier(n_estimators=300, random_state=SEED, n_jobs=1),
        feature_count,
        scale=False,
    )
    xgboost = make_pipeline(
        XGBClassifier(
            n_estimators=200,
            learning_rate=0.05,
            subsample=0.8,
            colsample_bytree=0.8,
            objective="binary:logistic",
            eval_metric="logloss",
            random_state=SEED,
            n_jobs=1,
        ),
        feature_count,
        scale=False,
    )
    gradient_boosting = make_pipeline(
        GradientBoostingClassifier(random_state=SEED),
        feature_count,
        scale=False,
    )
    return [
        ModelSpec(
            "logistic_regression",
            logistic,
            {
                "select__k": [5, 14, "all"],
                "model__C": [0.1, 1.0],
                "model__class_weight": [None, "balanced"],
            },
            True,
        ),
        ModelSpec(
            "regularized_logistic",
            regularized,
            {
                "select__k": [5, 14, "all"],
                "model__C": [0.03, 0.3, 1.0],
                "model__class_weight": [None, "balanced"],
            },
            True,
        ),
        ModelSpec(
            "linear_svm",
            linear_svm,
            {
                "select__k": [5, 14, "all"],
                "model__estimator__C": [0.1, 1.0],
            },
            True,
        ),
        ModelSpec(
            "calibrated_svm",
            calibrated_svm,
            {
                "select__k": [5, 14, "all"],
                "model__estimator__C": [0.1, 1.0],
            },
            True,
        ),
        ModelSpec(
            "random_forest",
            random_forest,
            {
                "select__k": [5, 14, "all"],
                "model__max_depth": [2, 4, None],
                "model__min_samples_leaf": [2, 5, 10],
                "model__class_weight": [None, "balanced"],
            },
            False,
        ),
        ModelSpec(
            "xgboost",
            xgboost,
            {
                "select__k": [5, 14, "all"],
                "model__max_depth": [1, 2, 3],
                "model__min_child_weight": [1, 5],
            },
            False,
        ),
        ModelSpec(
            "gradient_boosting",
            gradient_boosting,
            {
                "select__k": [5, 14, "all"],
                "model__n_estimators": [50, 100],
                "model__max_depth": [1, 2],
                "model__learning_rate": [0.03, 0.1],
            },
            False,
        ),
    ]


def probabilities(estimator: Any, frame: pd.DataFrame) -> np.ndarray:
    return np.asarray(estimator.predict_proba(frame)[:, 1], dtype=float)


def choose_threshold(y_true: np.ndarray, probability: np.ndarray) -> float:
    scores = [
        balanced_accuracy_score(y_true, probability >= threshold)
        for threshold in THRESHOLDS
    ]
    best = max(scores)
    candidates = [
        threshold
        for threshold, score in zip(THRESHOLDS, scores)
        if np.isclose(score, best)
    ]
    return float(min(candidates, key=lambda value: abs(value - 0.5)))


def metric_values(y_true: np.ndarray, probability: np.ndarray, threshold: float) -> dict[str, float]:
    prediction = (probability >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, prediction, labels=[0, 1]).ravel()
    return {
        "roc_auc": float(roc_auc_score(y_true, probability)),
        "pr_auc": float(average_precision_score(y_true, probability)),
        "accuracy": float(accuracy_score(y_true, prediction)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, prediction)),
        "sensitivity": float(recall_score(y_true, prediction, zero_division=0)),
        "specificity": float(tn / (tn + fp)) if tn + fp else np.nan,
        "ppv": float(precision_score(y_true, prediction, zero_division=0)),
        "npv": float(tn / (tn + fn)) if tn + fn else np.nan,
        "f1": float(f1_score(y_true, prediction, zero_division=0)),
        "brier_score": float(brier_score_loss(y_true, probability)),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
    }


def calibration_values(y_true: np.ndarray, probability: np.ndarray, model: str) -> tuple[dict[str, float], pd.DataFrame]:
    fraction, mean_probability = calibration_curve(
        y_true,
        probability,
        n_bins=10,
        strategy="uniform",
    )
    bins = np.clip((probability * 10).astype(int), 0, 9)
    rows = []
    for bin_number in range(10):
        mask = bins == bin_number
        if not mask.any():
            continue
        rows.append({
            "model": model,
            "bin": bin_number,
            "count": int(mask.sum()),
            "mean_predicted_probability": float(probability[mask].mean()),
            "observed_fraction_positive": float(y_true[mask].mean()),
        })
    calibration_errors = []
    calibration_weights = []
    for bin_number in range(10):
        mask = bins == bin_number
        if mask.any():
            calibration_errors.append(
                abs(probability[mask].mean() - y_true[mask].mean())
            )
            calibration_weights.append(int(mask.sum()))
    calibration_error = float(np.average(calibration_errors, weights=calibration_weights))
    max_error = float(max(
        abs(row["mean_predicted_probability"] - row["observed_fraction_positive"])
        for row in rows
    ))
    return {
        "expected_calibration_error": calibration_error,
        "maximum_calibration_error": max_error,
    }, pd.DataFrame(rows)


def bootstrap_intervals(y_true: np.ndarray, probability: np.ndarray, threshold: float, seed: int) -> dict[str, tuple[float, float]]:
    rng = np.random.default_rng(seed)
    values: dict[str, list[float]] = {}
    indices = np.arange(len(y_true))
    for _ in range(BOOTSTRAPS):
        sample = rng.choice(indices, size=len(indices), replace=True)
        if len(np.unique(y_true[sample])) < 2:
            continue
        current = metric_values(y_true[sample], probability[sample], threshold)
        for name, value in current.items():
            if name in {"tn", "fp", "fn", "tp"} or not np.isfinite(value):
                continue
            values.setdefault(name, []).append(value)
    return {
        name: (float(np.percentile(items, 2.5)), float(np.percentile(items, 97.5)))
        for name, items in values.items()
    }


def run_model(spec: ModelSpec, X: pd.DataFrame, y: np.ndarray, outer: StratifiedKFold) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    oof_probability = np.full(len(X), np.nan)
    oof_threshold = np.full(len(X), np.nan)
    fold_rows = []
    best_params = []
    for fold, (train_idx, test_idx) in enumerate(outer.split(X, y), start=1):
        inner = StratifiedKFold(n_splits=INNER_SPLITS, shuffle=True, random_state=SEED + fold)
        search = GridSearchCV(
            spec.estimator,
            spec.grid,
            scoring="roc_auc",
            cv=inner,
            n_jobs=1,
            refit=True,
            error_score="raise",
        )
        search.fit(X.iloc[train_idx], y[train_idx])
        threshold_cv = StratifiedKFold(n_splits=INNER_SPLITS, shuffle=True, random_state=SEED + 100 + fold)
        threshold_probability = cross_val_predict(
            clone(search.best_estimator_),
            X.iloc[train_idx],
            y[train_idx],
            cv=threshold_cv,
            method="predict_proba",
            n_jobs=1,
        )[:, 1]
        threshold = choose_threshold(y[train_idx], threshold_probability)
        probability = probabilities(search.best_estimator_, X.iloc[test_idx])
        oof_probability[test_idx] = probability
        oof_threshold[test_idx] = threshold
        best_params.append(search.best_params_)
        fold_rows.append({
            "model": spec.name,
            "fold": fold,
            "train_n": len(train_idx),
            "test_n": len(test_idx),
            "threshold": threshold,
            **metric_values(y[test_idx], probability, threshold),
        })
    threshold = float(np.median(oof_threshold))
    metrics = metric_values(y, oof_probability, threshold)
    calibration_metrics, calibration = calibration_values(y, oof_probability, spec.name)
    metrics.update(calibration_metrics)
    intervals = bootstrap_intervals(y, oof_probability, threshold, SEED + len(spec.name))
    for name, (low, high) in intervals.items():
        metrics[f"{name}_ci95_low"] = low
        metrics[f"{name}_ci95_high"] = high
    metrics.update({
        "model": spec.name,
        "status": "completed",
        "patients": len(y),
        "threshold": threshold,
        "outer_folds": OUTER_SPLITS,
        "inner_folds": INNER_SPLITS,
        "best_params_by_fold": json.dumps(best_params, sort_keys=True),
    })
    prediction = pd.DataFrame({
        "model": spec.name,
        "patient_id": X.index.astype(str),
        "outer_fold": np.nan,
        "y_true": y,
        "probability": oof_probability,
        "threshold": oof_threshold,
        "prediction": (oof_probability >= oof_threshold).astype(int),
    })
    for fold, (_, test_idx) in enumerate(outer.split(X, y), start=1):
        prediction.loc[test_idx, "outer_fold"] = fold
    calibration["threshold"] = threshold
    return metrics, prediction, calibration


def unavailable_neural_row(reason: str, patients: int) -> dict[str, Any]:
    return {
        "model": "gated_fusion_neural",
        "status": "unavailable",
        "patients": patients,
        "reason": reason,
    }


def main() -> None:
    seed_everything(SEED)
    frame = pd.read_csv(DATA)
    missing = [column for column in FEATURES + ["binary_label"] if column not in frame]
    if missing:
        raise RuntimeError(f"Missing required baseline columns: {missing}")
    if frame["sample_id_normalized"].duplicated().any():
        raise RuntimeError("Duplicate patient IDs in baseline cohort")
    X = frame[FEATURES].copy()
    y = frame["binary_label"].astype(int).to_numpy()
    outer = StratifiedKFold(n_splits=OUTER_SPLITS, shuffle=True, random_state=SEED)
    rows = []
    predictions = []
    calibration = []
    for spec in model_specs(len(FEATURES)):
        print(f"Running {spec.name}...", flush=True)
        metrics, prediction, curve = run_model(spec, X, y, outer)
        rows.append(metrics)
        prediction["patient_id"] = frame["sample_id_normalized"].astype(str).to_numpy()
        predictions.append(prediction)
        calibration.append(curve)
    rows.append(unavailable_neural_row(
        "Existing src/models/gated_fusion.py requires optional torch, which is not installed; no neural result was substituted.",
        len(frame),
    ))
    REPORTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(COMPARISON, index=False)
    pd.concat(predictions, ignore_index=True).to_csv(PREDICTIONS, index=False)
    pd.concat(calibration, ignore_index=True).to_csv(CALIBRATION, index=False)
    METADATA.write_text(json.dumps({
        "script": str(Path(__file__).relative_to(PROJECT)),
        "data": str(DATA.relative_to(PROJECT)),
        "endpoint": "GSE94533 cancer versus non-cancer",
        "patients": len(frame),
        "class_counts": {"negative": int((y == 0).sum()), "positive": int((y == 1).sum())},
        "seed": SEED,
        "outer_folds": OUTER_SPLITS,
        "inner_folds": INNER_SPLITS,
        "bootstrap_replicates": BOOTSTRAPS,
        "threshold_policy": "maximize balanced accuracy on inner training OOF probabilities; ties closest to 0.5",
        "external_data_used": False,
        "external_validation_policy": "GSE106817 and all other external cohorts remain untouched by this baseline run",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
    }, indent=2) + "\n")
    print(f"Comparison: {COMPARISON}")
    print(f"Predictions: {PREDICTIONS}")
    print(f"Calibration: {CALIBRATION}")


if __name__ == "__main__":
    main()