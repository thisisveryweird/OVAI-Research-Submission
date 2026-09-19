#!/usr/bin/env python3

from pathlib import Path
import json

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


ROOT = Path.home() / "ovarian-ai"

DATA = (
    ROOT
    / "data"
    / "processed"
    / "gse106817_external_validation.csv.gz"
)

OUT_PRED = (
    ROOT
    / "reports"
    / "tables"
    / "gse106817_170mirna_oof_predictions.csv"
)

OUT_RESULTS = (
    ROOT
    / "reports"
    / "tables"
    / "gse106817_170mirna_baseline_results.json"
)


print("=" * 80)
print("GSE106817 INDEPENDENT 170-miRNA BASELINE")
print("=" * 80)

df = pd.read_csv(DATA)

feature_cols = [
    c for c in df.columns
    if c.startswith("hsa-")
]

X = df[feature_cols].apply(
    pd.to_numeric,
    errors="coerce",
)

y = df["binary_label"].astype(int)

print()
print("Patients:", len(df))
print("Features:", len(feature_cols))
print("Cancer:", int(y.sum()))
print("Non-cancer:", int((y == 0).sum()))

# ---------------------------------------------------------------
# Remove rows with missing molecular data.
# ---------------------------------------------------------------

valid = (
    X.notna().all(axis=1)
    & y.notna()
)

X = X.loc[valid].reset_index(drop=True)
y = y.loc[valid].reset_index(drop=True)
meta = df.loc[valid].reset_index(drop=True)

print()
print("Patients after feature validation:", len(X))
print("Missing miRNA rows:", int((~valid).sum()))

# ---------------------------------------------------------------
# Model
#
# Class weighting is important here because GSE106817 has many
# non-cancer controls.
# ---------------------------------------------------------------

model = Pipeline([
    (
        "scaler",
        StandardScaler(),
    ),
    (
        "classifier",
        LogisticRegression(
            max_iter=5000,
            class_weight="balanced",
            solver="liblinear",
            random_state=42,
        ),
    ),
])

cv = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=42,
)

oof_prob = np.full(
    len(X),
    np.nan,
)

oof_pred = np.full(
    len(X),
    np.nan,
)

print()
print("=" * 80)
print("5-FOLD PATIENT-LEVEL OOF VALIDATION")
print("=" * 80)

for fold, (train_idx, test_idx) in enumerate(
    cv.split(X, y),
    start=1,
):

    X_train = X.iloc[train_idx]
    X_test = X.iloc[test_idx]

    y_train = y.iloc[train_idx]

    model.fit(
        X_train,
        y_train,
    )

    prob = model.predict_proba(
        X_test
    )[:, 1]

    pred = (
        prob >= 0.5
    ).astype(int)

    oof_prob[test_idx] = prob
    oof_pred[test_idx] = pred

    auc = roc_auc_score(
        y.iloc[test_idx],
        prob,
    )

    acc = accuracy_score(
        y.iloc[test_idx],
        pred,
    )

    print(
        f"Fold {fold}: "
        f"train={len(train_idx)} "
        f"test={len(test_idx)} "
        f"AUC={auc:.6f} "
        f"accuracy={acc:.6f}"
    )

# ---------------------------------------------------------------
# Overall OOF metrics
# ---------------------------------------------------------------

auc = roc_auc_score(
    y,
    oof_prob,
)

pr_auc = average_precision_score(
    y,
    oof_prob,
)

accuracy = accuracy_score(
    y,
    oof_pred,
)

cm = confusion_matrix(
    y,
    oof_pred,
)

tn, fp, fn, tp = cm.ravel()

sensitivity = (
    tp / (tp + fn)
    if (tp + fn)
    else np.nan
)

specificity = (
    tn / (tn + fp)
    if (tn + fp)
    else np.nan
)

print()
print("=" * 80)
print("RESULTS")
print("=" * 80)

print(
    f"OOF ROC AUC:    {auc:.6f}"
)

print(
    f"OOF PR AUC:     {pr_auc:.6f}"
)

print(
    f"OOF Accuracy:   {accuracy:.6f}"
)

print(
    f"Sensitivity:    {sensitivity:.6f}"
)

print(
    f"Specificity:    {specificity:.6f}"
)

print()
print("Confusion matrix:")
print(
    f"  TP={tp} FN={fn}"
)
print(
    f"  FP={fp} TN={tn}"
)

# ---------------------------------------------------------------
# Save OOF predictions
# ---------------------------------------------------------------

predictions = meta[
    [
        "sample_id",
        "title",
        "group",
        "age",
        "binary_label",
    ]
].copy()

predictions["prediction_probability"] = (
    oof_prob
)

predictions["prediction"] = (
    oof_pred.astype(int)
)

predictions["correct"] = (
    predictions["prediction"]
    == predictions["binary_label"]
)

OUT_PRED.parent.mkdir(
    parents=True,
    exist_ok=True,
)

predictions.to_csv(
    OUT_PRED,
    index=False,
)

results = {
    "dataset": "GSE106817",
    "patients": int(len(X)),
    "features": int(len(feature_cols)),
    "cancer": int(y.sum()),
    "non_cancer": int((y == 0).sum()),
    "roc_auc": float(auc),
    "pr_auc": float(pr_auc),
    "accuracy": float(accuracy),
    "sensitivity": float(sensitivity),
    "specificity": float(specificity),
    "tp": int(tp),
    "fn": int(fn),
    "fp": int(fp),
    "tn": int(tn),
    "model": "balanced logistic regression",
    "evaluation": "5-fold stratified patient-level OOF",
}

OUT_RESULTS.write_text(
    json.dumps(
        results,
        indent=2,
    )
)

print()
print("=" * 80)
print("FILES")
print("=" * 80)

print(OUT_PRED)
print(OUT_RESULTS)

print()
print("IMPORTANT:")
print(
    "This is an independent GSE106817 baseline."
)
print(
    "It is NOT the frozen GSE94533 model."
)
print(
    "Age is not used because most GSE106817 controls lack age."
)
