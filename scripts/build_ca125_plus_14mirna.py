#!/usr/bin/env python3

from __future__ import annotations

import json
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.model_selection import (
    GridSearchCV,
    StratifiedKFold,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path.home() / "ovarian-ai"

PREDICTIONS = (
    ROOT
    / "reports"
    / "tables"
    / "original_historical_clean_predictions.csv"
)

DATASET = (
    ROOT
    / "data"
    / "processed"
    / "gse94533_14mirna_benchmark.csv.gz"
)

CA125_PREDICTIONS = (
    ROOT
    / "reports"
    / "tables"
    / "ca125_baseline_oof_predictions.csv"
)

OUTPUT_DIR = (
    ROOT
    / "reports"
    / "tables"
)

OUT = (
    OUTPUT_DIR
    / "ca125_plus_14mirna_oof_predictions.csv"
)

RESULTS = (
    OUTPUT_DIR
    / "ca125_plus_14mirna_results.json"
)

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


def norm_id(x):
    s = str(x).strip().lower()
    if re.fullmatch(r"\d+\.0+", s):
        s = s.split(".", 1)[0]
    return re.sub(r"[^a-z0-9]", "", s)


print("=" * 90)
print("CA-125 + 14-miRNA MODERN BASELINE")
print("=" * 90)

pred = pd.read_csv(PREDICTIONS)

if len(pred) != 179:
    raise RuntimeError(f"Expected 179 patients, found {len(pred)}.")

pred["sample_key"] = pred["sample_id"].map(norm_id)

# ---------------------------------------------------------------
# Load exact 14-miRNA historical representation.
# ---------------------------------------------------------------

data = pd.read_csv(DATASET)

if len(data) != 179:
    raise RuntimeError(f"Expected 179 dataset rows, found {len(data)}.")

missing_features = [
    f for f in FEATURES
    if f not in data.columns
]

if missing_features:
    raise RuntimeError(
        f"Missing 14-miRNA features: {missing_features}"
    )

data["sample_key"] = data["sample_id"].map(norm_id)

data = data[
    ["sample_key"] + FEATURES
]

# ---------------------------------------------------------------
# Load CA-125.
# ---------------------------------------------------------------

ca = pd.read_csv(
    CA125_PREDICTIONS
)

if len(ca) != 179:
    raise RuntimeError(
        f"Expected 179 CA-125 rows, found {len(ca)}."
    )

ca["sample_key"] = ca["sample_id"].map(norm_id)

if "ca125" not in ca.columns:
    raise RuntimeError(
        "CA-125 prediction file does not contain 'ca125'."
    )

# ---------------------------------------------------------------
# Merge everything.
# ---------------------------------------------------------------

df = (
    pred[
        [
            "sample_id",
            "sample_key",
            "true_label",
        ]
    ]
    .merge(
        data,
        on="sample_key",
        how="left",
        validate="one_to_one",
    )
    .merge(
        ca[
            [
                "sample_key",
                "ca125",
            ]
        ],
        on="sample_key",
        how="left",
        validate="one_to_one",
    )
)

if len(df) != 179:
    raise RuntimeError(
        f"Expected 179 merged patients, got {len(df)}."
    )

# ---------------------------------------------------------------
# Build matrix.
# ---------------------------------------------------------------

X = df[
    ["ca125"] + FEATURES
].to_numpy(
    dtype=float
)

y = df[
    "true_label"
].to_numpy(
    dtype=int
)

if (y == 1).sum() != 97:
    raise RuntimeError("Expected 97 cancer patients.")

if (y == 0).sum() != 82:
    raise RuntimeError("Expected 82 non-cancer patients.")

print()
print("Patients:", len(df))
print("Features:", X.shape[1])
print("CA-125 + 14 miRNAs")
print("Cancer:", int((y == 1).sum()))
print("Non-cancer:", int((y == 0).sum()))

# ---------------------------------------------------------------
# Nested patient-level CV.
# ---------------------------------------------------------------

outer = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=2026,
)

oof_probability = np.full(
    len(y),
    np.nan,
)

oof_prediction = np.full(
    len(y),
    -1,
    dtype=int,
)

fold_records = []

for fold, (train_idx, test_idx) in enumerate(
    outer.split(X, y),
    start=1,
):

    print()
    print(f"FOLD {fold}/5")

    inner = StratifiedKFold(
        n_splits=3,
        shuffle=True,
        random_state=2026 + fold,
    )

    pipe = Pipeline(
        steps=[
            (
                "imputer",
                SimpleImputer(strategy="median"),
            ),
            (
                "scaler",
                StandardScaler(),
            ),
            (
                "classifier",
                LogisticRegression(
                    solver="saga",
                    max_iter=10000,
                    random_state=2026,
                ),
            ),
        ]
    )

    grid = {
        "classifier__C": [
            0.01,
            0.03,
            0.1,
            0.3,
            1.0,
            3.0,
        ],
        "classifier__l1_ratio": [
            0.0,
            0.25,
            0.5,
            0.75,
            1.0,
        ],
        "classifier__penalty": [
            "elasticnet",
        ],
    }

    search = GridSearchCV(
        pipe,
        grid,
        scoring="roc_auc",
        cv=inner,
        n_jobs=-1,
        refit=True,
    )

    search.fit(
        X[train_idx],
        y[train_idx],
    )

    probability = search.predict_proba(
        X[test_idx]
    )[:, 1]

    prediction = (
        probability >= 0.5
    ).astype(int)

    oof_probability[test_idx] = probability
    oof_prediction[test_idx] = prediction

    fold_auc = roc_auc_score(
        y[test_idx],
        probability,
    )

    fold_accuracy = accuracy_score(
        y[test_idx],
        prediction,
    )

    fold_records.append(
        {
            "fold": fold,
            "test_auc": float(fold_auc),
            "test_accuracy": float(fold_accuracy),
            "inner_cv_auc": float(search.best_score_),
            "best_params": search.best_params_,
        }
    )

    print(
        f"  Test AUC: {fold_auc:.6f}"
    )

    print(
        f"  Accuracy: {fold_accuracy:.6f}"
    )

    print(
        f"  Best params: {search.best_params_}"
    )

# ---------------------------------------------------------------
# Final OOF metrics.
# ---------------------------------------------------------------

auc = roc_auc_score(
    y,
    oof_probability,
)

pr_auc = average_precision_score(
    y,
    oof_probability,
)

accuracy = accuracy_score(
    y,
    oof_prediction,
)

tn, fp, fn, tp = confusion_matrix(
    y,
    oof_prediction,
).ravel()

sensitivity = tp / (tp + fn)
specificity = tn / (tn + fp)

print()
print("=" * 90)
print("CA-125 + 14-miRNA RESULTS")
print("=" * 90)

print()
print(f"OOF ROC AUC:   {auc:.6f}")
print(f"OOF PR AUC:    {pr_auc:.6f}")
print(f"OOF Accuracy:  {accuracy:.6f}")
print(f"Sensitivity:   {sensitivity:.6f}")
print(f"Specificity:   {specificity:.6f}")

print()
print("Confusion matrix:")
print(f"  TP={tp} FN={fn}")
print(f"  FP={fp} TN={tn}")

# ---------------------------------------------------------------
# Save.
# ---------------------------------------------------------------

df["oof_probability"] = oof_probability
df["oof_prediction"] = oof_prediction

df.to_csv(
    OUT,
    index=False,
)

results = {
    "model": "CA-125 + 14-miRNA elastic-net logistic regression",
    "patients": 179,
    "cancer": 97,
    "non_cancer": 82,
    "n_features": 15,
    "cross_validation": "5-fold nested patient-level CV",
    "metrics": {
        "roc_auc": float(auc),
        "pr_auc": float(pr_auc),
        "accuracy": float(accuracy),
        "sensitivity": float(sensitivity),
        "specificity": float(specificity),
    },
    "confusion_matrix": {
        "TN": int(tn),
        "FP": int(fp),
        "FN": int(fn),
        "TP": int(tp),
    },
    "folds": fold_records,
}

with RESULTS.open("w") as f:
    json.dump(
        results,
        f,
        indent=2,
    )

print()
print("Saved:")
print(OUT)
print(RESULTS)
