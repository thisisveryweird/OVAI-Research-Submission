#!/usr/bin/env python3

"""
FULL HISTORICAL miRNA BASELINE

Uses the full miRNA feature set available in eLife Supplemental
Dataset 2 on the exact historical 179-patient cohort.

This is a NEW model.

It does NOT use:
    - historical PMML weights
    - historical neural-network architecture

Model:
    StandardScaler
        ->
    Elastic-net logistic regression

Evaluation:
    5-fold outer patient-level stratified CV

Hyperparameter tuning:
    3-fold inner CV within each outer training set

Therefore:

    test patients
        are never used for
        feature/model selection.

This prevents data leakage.

The purpose is to compare:

    historical 14-miRNA neural network
                vs
    modern 14-miRNA logistic regression
                vs
    modern full-miRNA regularized model
"""

from __future__ import annotations

import json
import re
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


# =====================================================================
# PATHS
# =====================================================================

ROOT = Path.home() / "ovarian-ai"

SUPPLEMENTARY_XLSX = (
    ROOT
    / "data"
    / "external"
    / "elife_28932"
    / "elife_poa_e28932_Supplementary_File_6.xlsx"
)

HISTORICAL_PREDICTIONS = (
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

PREDICTIONS_OUTPUT = (
    OUTPUT_DIR
    / "full_historical_mirna_oof_predictions.csv"
)

RESULTS_OUTPUT = (
    OUTPUT_DIR
    / "full_historical_mirna_baseline_results.json"
)


# =====================================================================
# SETTINGS
# =====================================================================

OUTER_SPLITS = 5
INNER_SPLITS = 3

RANDOM_STATE = 2026


# =====================================================================
# ID NORMALIZATION
# =====================================================================

def norm_id(value) -> str:

    s = str(value).strip().lower()

    if re.fullmatch(r"\d+\.0+", s):
        s = s.split(".", 1)[0]

    return s.replace("-", "")


# =====================================================================
# IDENTIFY miRNA FEATURE ROWS
# =====================================================================

def looks_like_mirna(value) -> bool:

    if pd.isna(value):
        return False

    s = str(value).strip()

    return (
        s.startswith("hsa-miR-")
        or s.startswith("hsa-let-")
    )


# =====================================================================
# LOAD HISTORICAL COHORT
# =====================================================================

print("=" * 90)
print("FULL HISTORICAL miRNA BASELINE")
print("=" * 90)

print()
print("Loading historical cohort...")

pred = pd.read_csv(
    HISTORICAL_PREDICTIONS
)

if len(pred) != 179:
    raise RuntimeError(
        f"Expected 179 historical patients, found {len(pred)}."
    )

pred["sample_key"] = (
    pred["sample_id"].map(norm_id)
)

if pred["sample_key"].duplicated().any():
    raise RuntimeError(
        "Duplicate historical sample IDs detected."
    )

print(
    f"Patients: {len(pred)}"
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

raw_header = list(
    raw.iloc[
        HEADER_ROW,
        1:180,
    ]
)

if len(raw_header) != 179:
    raise RuntimeError(
        f"Expected 179 sample columns, found {len(raw_header)}."
    )


# =====================================================================
# USE THE SAME HISTORICAL SAMPLE ORDER AS THE PMML REPRODUCTION
# =====================================================================

historical_ids = pred["sample_key"].tolist()

# Dataset 2's first 119 numeric columns correspond to the retained
# numeric historical samples; the remainder are the 60 e-samples.
#
# We already independently established the exact historical ordering
# in original_historical_clean_predictions.csv.

if len(historical_ids) != 179:
    raise RuntimeError(
        "Historical ID count is not 179."
    )


# =====================================================================
# FIND ALL miRNA ROWS
# =====================================================================

feature_names = []

feature_rows = []

for row_idx, value in enumerate(
    raw.iloc[:, 0]
):

    if not looks_like_mirna(value):
        continue

    feature = str(value).strip()

    # Avoid duplicate feature names.
    if feature in feature_names:
        continue

    feature_names.append(feature)
    feature_rows.append(row_idx)


if len(feature_names) < 50:
    raise RuntimeError(
        "Unexpectedly few miRNA features found: "
        f"{len(feature_names)}"
    )


print()
print("=" * 90)
print("FULL miRNA FEATURE SET")
print("=" * 90)

print()
print(
    f"miRNA features detected: {len(feature_names)}"
)

print()
print("First 20 features:")

for feature in feature_names[:20]:
    print(f"  {feature}")

print()
print("Last 10 features:")

for feature in feature_names[-10:]:
    print(f"  {feature}")


# =====================================================================
# EXTRACT MATRIX
# =====================================================================

values = []

for feature, row_idx in zip(
    feature_names,
    feature_rows,
):

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
        bad = np.where(
            ~np.isfinite(vals)
        )[0]

        raise RuntimeError(
            f"{feature} contains non-finite values "
            f"at columns {bad[:10].tolist()}."
        )

    values.append(vals)


X = np.vstack(values).T


# =====================================================================
# BUILD DATAFRAME
# =====================================================================

mirna = pd.DataFrame(
    X,
    columns=feature_names,
)

mirna.insert(
    0,
    "sample_id",
    historical_ids,
)

mirna["sample_key"] = (
    mirna["sample_id"]
    .map(norm_id)
)


# =====================================================================
# ALIGN LABELS
# =====================================================================

label_map = (
    pred
    .set_index("sample_key")["true_label"]
    .astype(int)
)

mirna["true_label"] = (
    mirna["sample_key"]
    .map(label_map)
)

if mirna["true_label"].isna().any():

    bad = mirna.loc[
        mirna["true_label"].isna(),
        "sample_id",
    ].tolist()

    raise RuntimeError(
        "Missing labels for:\n"
        + "\n".join(
            map(str, bad)
        )
    )

mirna["true_label"] = (
    mirna["true_label"].astype(int)
)


# =====================================================================
# VALIDATION
# =====================================================================

y = mirna["true_label"].to_numpy()

if X.shape[0] != 179:
    raise RuntimeError(
        f"Expected 179 patients, found {X.shape[0]}."
    )

if X.shape[1] != len(feature_names):
    raise RuntimeError(
        "Feature dimension mismatch."
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
print("DATA VALIDATION")
print("=" * 90)

print()
print("Patients:", X.shape[0])
print("miRNA features:", X.shape[1])
print("Cancer:", int((y == 1).sum()))
print("Non-cancer:", int((y == 0).sum()))


# =====================================================================
# OUTER CV
# =====================================================================

outer_cv = StratifiedKFold(
    n_splits=OUTER_SPLITS,
    shuffle=True,
    random_state=RANDOM_STATE,
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


print()
print("=" * 90)
print("NESTED CROSS-VALIDATION")
print("=" * 90)


for fold, (train_idx, test_idx) in enumerate(
    outer_cv.split(X, y),
    start=1,
):

    print()
    print(
        f"OUTER FOLD {fold}/{OUTER_SPLITS}"
    )

    X_train = X[train_idx]
    X_test = X[test_idx]

    y_train = y[train_idx]
    y_test = y[test_idx]


    # ---------------------------------------------------------------
    # Inner CV
    # ---------------------------------------------------------------

    inner_cv = StratifiedKFold(
        n_splits=INNER_SPLITS,
        shuffle=True,
        random_state=RANDOM_STATE + fold,
    )


    # ---------------------------------------------------------------
    # Pipeline
    #
    # Scaling occurs separately inside every inner/outer training
    # split.
    # ---------------------------------------------------------------

    pipeline = Pipeline(
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
                    random_state=RANDOM_STATE,
                ),
            ),
        ]
    )


    # ---------------------------------------------------------------
    # Hyperparameter grid
    # ---------------------------------------------------------------

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
        pipeline,
        param_grid=param_grid,
        scoring="roc_auc",
        cv=inner_cv,
        n_jobs=-1,
        refit=True,
    )


    search.fit(
        X_train,
        y_train,
    )


    probability = search.predict_proba(
        X_test
    )[:, 1]

    prediction = (
        probability >= 0.5
    ).astype(int)


    oof_probability[test_idx] = (
        probability
    )

    oof_prediction[test_idx] = (
        prediction
    )


    fold_auc = roc_auc_score(
        y_test,
        probability,
    )

    fold_accuracy = accuracy_score(
        y_test,
        prediction,
    )


    best_params = {
        key: value
        for key, value in search.best_params_.items()
    }


    fold_records.append(
        {
            "fold": fold,
            "n_train": int(len(train_idx)),
            "n_test": int(len(test_idx)),
            "roc_auc": float(fold_auc),
            "accuracy": float(fold_accuracy),
            "inner_cv_auc": float(
                search.best_score_
            ),
            "best_params": best_params,
        }
    )


    print(
        f"  test patients: {len(test_idx)}"
    )

    print(
        f"  inner-CV AUC: {search.best_score_:.6f}"
    )

    print(
        f"  test AUC:     {fold_auc:.6f}"
    )

    print(
        f"  accuracy:     {fold_accuracy:.6f}"
    )

    print(
        f"  best params:  {best_params}"
    )


# =====================================================================
# OOF METRICS
# =====================================================================

if np.any(
    ~np.isfinite(oof_probability)
):
    raise RuntimeError(
        "Some patients lack OOF probabilities."
    )

if np.any(
    oof_prediction < 0
):
    raise RuntimeError(
        "Some patients lack OOF predictions."
    )


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
# SAVE PREDICTIONS
# =====================================================================

mirna["oof_probability_cancer"] = (
    oof_probability
)

mirna["oof_prediction"] = (
    oof_prediction
)

mirna["fold"] = -1


for fold, (_, test_idx) in enumerate(
    outer_cv.split(X, y),
    start=1,
):

    mirna.loc[
        test_idx,
        "fold",
    ] = fold


OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

mirna.to_csv(
    PREDICTIONS_OUTPUT,
    index=False,
)


# =====================================================================
# SAVE RESULTS
# =====================================================================

results = {
    "model": "nested-CV elastic-net logistic regression",
    "representation": {
        "source": "eLife Supplemental Dataset 2",
        "batch_adjusted": True,
        "log10_transformed": True,
    },
    "cohort": {
        "patients": 179,
        "cancer": 97,
        "non_cancer": 82,
    },
    "n_mirna_features": int(
        len(feature_names)
    ),
    "cross_validation": {
        "outer": "5-fold stratified patient-level",
        "inner": "3-fold stratified",
        "random_state": RANDOM_STATE,
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
print("FULL HISTORICAL miRNA BASELINE RESULTS")
print("=" * 90)

print()
print(f"Patients:       {len(y)}")
print(f"miRNA features: {len(feature_names)}")

print()
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
print(PREDICTIONS_OUTPUT)
print(RESULTS_OUTPUT)

print()
print("=" * 90)
print("IMPORTANT")
print("=" * 90)

print()
print(
    "This is a NEW model and is not the historical eLife neural network."
)

print(
    "All scaling and hyperparameter selection were performed "
    "inside the cross-validation procedure."
)
