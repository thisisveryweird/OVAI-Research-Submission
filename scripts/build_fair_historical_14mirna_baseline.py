#!/usr/bin/env python3

"""
FAIR MODERN 14-miRNA BASELINE

Uses EXACTLY the same historical cohort and Dataset 2 representation
used by the reproduced eLife PMML model.

Historical cohort:
    179 patients
    97 cancer
    82 non-cancer

Representation:
    eLife Supplementary Dataset 2
    batch-adjusted
    log10-transformed

Features:
    exact 14 PMML miRNAs

Model:
    standardized L2 logistic regression

Evaluation:
    5-fold stratified patient-level OOF predictions

This is a NEW model.
It does NOT use the historical PMML weights.
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
from sklearn.model_selection import (
    GridSearchCV,
    StratifiedKFold,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path.home() / "ovarian-ai"

SUPPLEMENTARY_XLSX = (
    ROOT
    / "data"
    / "external"
    / "elife_28932"
    / "elife_poa_e28932_Supplementary_File_6.xlsx"
)

PREDICTIONS_PATH = (
    ROOT
    / "reports"
    / "tables"
    / "original_historical_clean_predictions.csv"
)

OUTPUT_DIR = (
    ROOT
    / "reports"
    / "tables"
)

PREDICTION_OUTPUT = (
    OUTPUT_DIR
    / "fair_historical_14mirna_oof_predictions.csv"
)

RESULTS_OUTPUT = (
    OUTPUT_DIR
    / "fair_historical_14mirna_baseline_results.json"
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


def norm_id(value):

    s = str(value).strip().lower()

    if s.endswith(".0") and s[:-2].isdigit():
        s = s[:-2]

    return s.replace("-", "")


# =====================================================================
# LOAD HISTORICAL SAMPLE ORDER
# =====================================================================

print("=" * 90)
print("FAIR HISTORICAL 14-miRNA MODERN BASELINE")
print("=" * 90)

print()
print("Loading historical PMML prediction file...")

pred = pd.read_csv(PREDICTIONS_PATH)

if len(pred) != 179:
    raise RuntimeError(
        f"Expected 179 historical patients, found {len(pred)}."
    )

pred["sample_key"] = pred["sample_id"].map(norm_id)

print(
    f"Historical cohort: {len(pred)} patients"
)

print(
    f"Cancer: {(pred['true_label'] == 1).sum()}"
)

print(
    f"Non-cancer: {(pred['true_label'] == 0).sum()}"
)


# =====================================================================
# READ DATASET 2
# =====================================================================

print()
print("=" * 90)
print("READING eLIFE SUPPLEMENTAL DATASET 2")
print("=" * 90)

raw = pd.read_excel(
    SUPPLEMENTARY_XLSX,
    sheet_name="Supplemental Dataset 2",
    header=None,
    engine="openpyxl",
)

print()
print("Raw shape:", raw.shape)

HEADER_ROW = 2

# Dataset 2:
#   column 0 = feature name
#   columns 1..179 = 179 retained historical samples

header = list(
    raw.iloc[
        HEADER_ROW,
        1:180,
    ]
)

if len(header) != 179:
    raise RuntimeError(
        f"Expected 179 Dataset 2 sample columns, found {len(header)}."
    )


# =====================================================================
# BUILD THE EXACT HISTORICAL SAMPLE ORDER FROM THE PMML REPRODUCTION
# =====================================================================

# The prediction file is already in the exact historical ordering
# established by the successful PMML reproduction.

historical_ids = pred["sample_key"].tolist()

if len(set(historical_ids)) != 179:
    raise RuntimeError(
        "Historical prediction file contains duplicate sample IDs."
    )


# =====================================================================
# EXTRACT THE 14 FEATURES
# =====================================================================

feature_series = (
    raw.iloc[:, 0]
    .astype(str)
    .str.strip()
)

rows = {}

for feature in FEATURES:

    matches = np.where(
        feature_series.values == feature
    )[0]

    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one row for {feature}, "
            f"found {len(matches)}."
        )

    rows[feature] = int(matches[0])


values = []

for feature in FEATURES:

    row_idx = rows[feature]

    vals = pd.to_numeric(
        raw.iloc[
            row_idx,
            1:180,
        ],
        errors="coerce",
    ).to_numpy(
        dtype=np.float64
    )

    if len(vals) != 179:
        raise RuntimeError(
            f"{feature} has {len(vals)} values."
        )

    if not np.all(np.isfinite(vals)):
        raise RuntimeError(
            f"{feature} contains non-finite values."
        )

    values.append(vals)


X = np.vstack(values).T

historical = pd.DataFrame(
    X,
    columns=FEATURES,
)

historical.insert(
    0,
    "sample_id",
    historical_ids,
)


# =====================================================================
# ALIGN LABELS
# =====================================================================

label_map = (
    pred
    .set_index("sample_key")["true_label"]
    .astype(int)
)

historical["sample_key"] = (
    historical["sample_id"]
    .map(norm_id)
)

historical["true_label"] = (
    historical["sample_key"]
    .map(label_map)
)

if historical["true_label"].isna().any():
    bad = historical.loc[
        historical["true_label"].isna(),
        "sample_id",
    ].tolist()

    raise RuntimeError(
        "Missing labels for historical samples:\n"
        + "\n".join(map(str, bad))
    )

historical["true_label"] = (
    historical["true_label"].astype(int)
)


# =====================================================================
# VALIDATION
# =====================================================================

y = historical["true_label"].to_numpy()

if len(historical) != 179:
    raise RuntimeError(
        "Historical dataset is not 179 patients."
    )

if X.shape != (179, 14):
    raise RuntimeError(
        f"Unexpected X shape: {X.shape}"
    )

if (y == 1).sum() != 97:
    raise RuntimeError(
        f"Expected 97 cancer patients, found {(y == 1).sum()}."
    )

if (y == 0).sum() != 82:
    raise RuntimeError(
        f"Expected 82 non-cancer patients, found {(y == 0).sum()}."
    )


print()
print("=" * 90)
print("HISTORICAL DATA VALIDATION")
print("=" * 90)

print()
print("Patients:", len(historical))
print("Features:", X.shape[1])
print("Cancer:", int((y == 1).sum()))
print("Non-cancer:", int((y == 0).sum()))


# =====================================================================
# 5-FOLD OOF
# =====================================================================

cv = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=2026,
)

oof_probability = np.full(
    len(historical),
    np.nan,
)

oof_prediction = np.full(
    len(historical),
    -1,
    dtype=int,
)

fold_records = []


print()
print("=" * 90)
print("5-FOLD PATIENT-LEVEL CROSS-VALIDATION")
print("=" * 90)


for fold, (train_idx, test_idx) in enumerate(
    cv.split(X, y),
    start=1,
):

    model = Pipeline(
        steps=[
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

    inner_cv = StratifiedKFold(
        n_splits=3,
        shuffle=True,
        random_state=2026 + fold,
    )

    param_grid = {
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
        model,
        param_grid=param_grid,
        scoring="roc_auc",
        cv=inner_cv,
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
            "n_train": int(len(train_idx)),
            "n_test": int(len(test_idx)),
            "roc_auc": float(fold_auc),
            "accuracy": float(fold_accuracy),
        }
    )

    print(
        f"Fold {fold}: "
        f"train={len(train_idx)} "
        f"test={len(test_idx)} "
        f"AUC={fold_auc:.6f} "
        f"accuracy={fold_accuracy:.6f}"
    )


# =====================================================================
# OOF METRICS
# =====================================================================

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


# =====================================================================
# SAVE
# =====================================================================

historical["oof_probability_cancer"] = (
    oof_probability
)

historical["oof_prediction"] = (
    oof_prediction
)

historical["fold"] = -1

for fold, (_, test_idx) in enumerate(
    cv.split(X, y),
    start=1,
):

    historical.loc[
        test_idx,
        "fold",
    ] = fold


OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

historical.to_csv(
    PREDICTION_OUTPUT,
    index=False,
)


results = {
    "model": "standardized L2 logistic regression",
    "cohort": {
        "patients": 179,
        "cancer": 97,
        "non_cancer": 82,
    },
    "representation": {
        "source": "eLife Supplementary Dataset 2",
        "batch_adjusted": True,
        "log10_transformed": True,
    },
    "features": FEATURES,
    "cross_validation": {
        "method": "5-fold stratified patient-level OOF",
        "random_state": 2026,
    },
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

with RESULTS_OUTPUT.open("w") as f:

    json.dump(
        results,
        f,
        indent=2,
    )


print()
print("=" * 90)
print("FAIR HISTORICAL 14-miRNA BASELINE RESULTS")
print("=" * 90)

print()
print(f"Patients:       {len(historical)}")
print(f"Cancer:         {(y == 1).sum()}")
print(f"Non-cancer:     {(y == 0).sum()}")

print()
print(f"OOF ROC AUC:    {auc:.6f}")
print(f"OOF PR AUC:     {pr_auc:.6f}")
print(f"OOF Accuracy:   {accuracy:.6f}")
print(f"Sensitivity:    {sensitivity:.6f}")
print(f"Specificity:    {specificity:.6f}")

print()
print("Confusion matrix:")
print(f"  TP={tp} FN={fn}")
print(f"  FP={fp} TN={tn}")

print()
print("Saved:")
print(PREDICTION_OUTPUT)
print(RESULTS_OUTPUT)
