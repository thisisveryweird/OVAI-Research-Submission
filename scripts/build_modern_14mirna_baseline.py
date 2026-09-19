#!/usr/bin/env python3

"""
MODERN 14-miRNA BASELINE
========================

Build a leakage-safe modern baseline using the exact 14 miRNAs from
the historical eLife PMML model.

This script DOES NOT use the historical neural-network weights.

It trains a new regularized logistic-regression model using:

    exact 14 historical miRNAs
        ->
    StandardScaler
        ->
    LogisticRegression

Evaluation:

    5-fold stratified patient-level cross-validation

All preprocessing occurs INSIDE each fold.

The output is an out-of-fold prediction for every patient.

This provides a fair modern baseline against which the fixed
historical PMML network can later be compared.

No neural network is trained here.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


# =====================================================================
# PATHS
# =====================================================================

ROOT = Path.home() / "ovarian-ai"

DATA_PATH = (
    ROOT
    / "data"
    / "processed"
    / "gse94533_14mirna_benchmark.csv.gz"
)

OUTPUT_DIR = (
    ROOT
    / "reports"
    / "tables"
)

PREDICTIONS_OUTPUT = (
    OUTPUT_DIR
    / "modern_14mirna_oof_predictions.csv"
)

RESULTS_OUTPUT = (
    OUTPUT_DIR
    / "modern_14mirna_baseline_results.json"
)


# =====================================================================
# EXACT HISTORICAL FEATURES
# =====================================================================

FEATURES = [
    "hsa-miR-23b-3p",
    "hsa-miR-29a-3p",
    "hsa-miR-32-5p",
    "hsa-miR-92a-3p",
    "hsa-miR-150-5p",
    "hsa-miR-200a-3p",
    "hsa-miR-200c-3p",
    "hsa-miR-203a",
    "hsa-miR-320c",
    "hsa-miR-320d",
    "hsa-miR-335-5p",
    "hsa-miR-450b-5p",
    "hsa-miR-1246",
    "hsa-miR-1307-5p",
]


# =====================================================================
# LOAD DATA
# =====================================================================

print("=" * 90)
print("MODERN 14-miRNA BASELINE")
print("=" * 90)

print()
print("Loading:")
print(DATA_PATH)

if not DATA_PATH.exists():
    raise RuntimeError(
        f"Dataset not found:\n{DATA_PATH}"
    )

df = pd.read_csv(DATA_PATH)

print()
print("Dataset shape:", df.shape)

print()
print("Columns:")
print(df.columns.tolist())


# =====================================================================
# IDENTIFY LABEL COLUMN
# =====================================================================

label_candidates = [
    "true_label",
    "binary_label",
    "y",
    "label",
]

label_column = None

for candidate in label_candidates:

    if candidate in df.columns:
        label_column = candidate
        break

if label_column is None:
    raise RuntimeError(
        "Could not find binary label column. "
        f"Available columns: {df.columns.tolist()}"
    )


# =====================================================================
# VALIDATE FEATURES
# =====================================================================

missing_features = [
    feature
    for feature in FEATURES
    if feature not in df.columns
]

if missing_features:
    raise RuntimeError(
        "Missing required 14-miRNA features:\n"
        + "\n".join(missing_features)
    )


# =====================================================================
# BUILD X / y
# =====================================================================

X = df[FEATURES].to_numpy(
    dtype=np.float64
)

y = df[label_column].to_numpy(
    dtype=int
)

if len(df) != 179:
    raise RuntimeError(
        f"Expected 179 patients, found {len(df)}."
    )

if len(np.unique(y)) != 2:
    raise RuntimeError(
        "Expected exactly two classes."
    )

if not np.all(np.isfinite(X)):
    raise RuntimeError(
        "Expression matrix contains non-finite values."
    )

print()
print("=" * 90)
print("DATA VALIDATION")
print("=" * 90)

print()
print(f"Patients: {len(df)}")
print(f"Features: {X.shape[1]}")
print(f"Cancer: {(y == 1).sum()}")
print(f"Non-cancer: {(y == 0).sum()}")


# =====================================================================
# CROSS-VALIDATION
# =====================================================================

cv = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=2026,
)

oof_probability = np.full(
    len(df),
    np.nan,
)

oof_prediction = np.full(
    len(df),
    -1,
    dtype=int,
)

fold_records = []


print()
print("=" * 90)
print("5-FOLD STRATIFIED PATIENT-LEVEL CROSS-VALIDATION")
print("=" * 90)


for fold, (train_idx, test_idx) in enumerate(
    cv.split(X, y),
    start=1,
):

    X_train = X[train_idx]
    X_test = X[test_idx]

    y_train = y[train_idx]
    y_test = y[test_idx]

    # ---------------------------------------------------------------
    # Pipeline prevents leakage:
    #
    # scaler is fitted ONLY on training data.
    # ---------------------------------------------------------------

    model = Pipeline(
        steps=[
            (
                "scaler",
                StandardScaler(),
            ),
            (
                "classifier",
                LogisticRegression(
                    penalty="l2",
                    C=1.0,
                    solver="liblinear",
                    max_iter=5000,
                    random_state=2026,
                ),
            ),
        ]
    )

    model.fit(
        X_train,
        y_train,
    )

    probability = model.predict_proba(
        X_test
    )[:, 1]

    prediction = (
        probability >= 0.5
    ).astype(int)

    oof_probability[test_idx] = probability
    oof_prediction[test_idx] = prediction

    fold_auc = roc_auc_score(
        y_test,
        probability,
    )

    fold_accuracy = accuracy_score(
        y_test,
        prediction,
    )

    fold_records.append(
        {
            "fold": fold,
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "cancer_train": int((y_train == 1).sum()),
            "noncancer_train": int((y_train == 0).sum()),
            "cancer_test": int((y_test == 1).sum()),
            "noncancer_test": int((y_test == 0).sum()),
            "roc_auc": float(fold_auc),
            "accuracy": float(fold_accuracy),
        }
    )

    print()
    print(
        f"Fold {fold}: "
        f"train={len(train_idx)}, "
        f"test={len(test_idx)}, "
        f"AUC={fold_auc:.6f}, "
        f"accuracy={fold_accuracy:.6f}"
    )


# =====================================================================
# OUT-OF-FOLD METRICS
# =====================================================================

if np.any(~np.isfinite(oof_probability)):
    raise RuntimeError(
        "Some patients did not receive an out-of-fold prediction."
    )

if np.any(oof_prediction < 0):
    raise RuntimeError(
        "Some patients did not receive an out-of-fold class prediction."
    )


oof_auc = roc_auc_score(
    y,
    oof_probability,
)

oof_pr_auc = average_precision_score(
    y,
    oof_probability,
)

oof_accuracy = accuracy_score(
    y,
    oof_prediction,
)

cm = confusion_matrix(
    y,
    oof_prediction,
)

tn, fp, fn, tp = cm.ravel()

sensitivity = (
    tp / (tp + fn)
    if (tp + fn) > 0
    else float("nan")
)

specificity = (
    tn / (tn + fp)
    if (tn + fp) > 0
    else float("nan")
)


# =====================================================================
# SAVE PREDICTIONS
# =====================================================================

predictions = df.copy()

predictions["oof_probability_cancer"] = (
    oof_probability
)

predictions["oof_prediction"] = (
    oof_prediction
)

predictions["fold"] = -1

for fold, (_, test_idx) in enumerate(
    cv.split(X, y),
    start=1,
):
    predictions.loc[
        test_idx,
        "fold",
    ] = fold


OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

predictions.to_csv(
    PREDICTIONS_OUTPUT,
    index=False,
)


# =====================================================================
# SAVE RESULTS
# =====================================================================

results = {
    "model": "standardized L2 logistic regression",
    "features": FEATURES,
    "n_patients": int(len(df)),
    "n_features": int(len(FEATURES)),
    "cancer_count": int((y == 1).sum()),
    "noncancer_count": int((y == 0).sum()),
    "cross_validation": {
        "method": "5-fold stratified patient-level CV",
        "shuffle": True,
        "random_state": 2026,
    },
    "metrics": {
        "roc_auc": float(oof_auc),
        "pr_auc": float(oof_pr_auc),
        "accuracy": float(oof_accuracy),
        "sensitivity_at_0.5": float(sensitivity),
        "specificity_at_0.5": float(specificity),
    },
    "confusion_matrix": {
        "TN": int(tn),
        "FP": int(fp),
        "FN": int(fn),
        "TP": int(tp),
    },
    "folds": fold_records,
}

with RESULTS_OUTPUT.open(
    "w"
) as f:

    json.dump(
        results,
        f,
        indent=2,
    )


# =====================================================================
# REPORT
# =====================================================================

print()
print("=" * 90)
print("MODERN 14-miRNA BASELINE RESULTS")
print("=" * 90)

print()
print(f"Patients:       {len(df)}")
print(f"Cancer:         {(y == 1).sum()}")
print(f"Non-cancer:     {(y == 0).sum()}")

print()
print(f"OOF ROC AUC:    {oof_auc:.6f}")
print(f"OOF PR AUC:     {oof_pr_auc:.6f}")
print(f"OOF Accuracy:   {oof_accuracy:.6f}")
print(f"Sensitivity:    {sensitivity:.6f}")
print(f"Specificity:    {specificity:.6f}")

print()
print("Confusion matrix:")
print(
    f"  TP={tp} FN={fn}"
)
print(
    f"  FP={fp} TN={tn}"
)

print()
print("=" * 90)
print("FILES")
print("=" * 90)

print()
print("OOF predictions:")
print(PREDICTIONS_OUTPUT)

print()
print("Results:")
print(RESULTS_OUTPUT)

print()
print("=" * 90)
print("IMPORTANT")
print("=" * 90)

print()
print(
    "This is a NEW modern model. It is NOT the historical "
    "eLife neural network."
)

print(
    "All scaling was performed separately inside each training fold."
)

print(
    "The out-of-fold AUC is therefore a patient-level cross-validated "
    "estimate rather than a resubstitution/full-cohort score."
)
