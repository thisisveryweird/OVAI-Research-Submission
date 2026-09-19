#!/usr/bin/env python3

from pathlib import Path
import json
import re

import numpy as np
import pandas as pd
import xenaPython as xena

from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path.home() / "ovarian-ai"

DATA = (
    ROOT
    / "data"
    / "processed"
    / "tcga_ov_223_multimodal_protein.csv.gz"
)

OUT_PRED = (
    ROOT
    / "reports"
    / "tables"
    / "tcga_multimodal_stage_oof_predictions.csv"
)

OUT_RESULTS = (
    ROOT
    / "reports"
    / "tables"
    / "tcga_multimodal_stage_results.json"
)

HUB = "https://tcga.xenahubs.net"

CLINICAL = (
    "TCGA.OV.sampleMap/"
    "OV_clinicalMatrix"
)

RNA = [
    "PAX8", "EPCAM", "KRT7", "KRT8", "KRT18",
    "MUC16", "WFDC2", "FOLR1", "WT1", "MKI67",
]

print("=" * 80)
print("TCGA-OV MULTIMODAL STAGE MODEL")
print("=" * 80)

# ------------------------------------------------------------------
# Load verified 223-patient multimodal dataset
# ------------------------------------------------------------------

df = pd.read_csv(DATA)
df["sample_id"] = df["sample_id"].astype(str)

print()
print("Input patients:", len(df))

# ------------------------------------------------------------------
# Pull clinical stage directly from TCGA
# ------------------------------------------------------------------

samples = df["sample_id"].tolist()

_, stage_values = xena.dataset_probe_values(
    HUB,
    CLINICAL,
    samples,
    ["clinical_stage"],
)

stage_raw = pd.Series(
    stage_values[0],
    index=df.index,
    dtype="object",
)

# ------------------------------------------------------------------
# TCGA clinical_stage is categorical and Xena returns numeric codes.
# Decode those codes using the field metadata supplied by Xena.
# ------------------------------------------------------------------

code_info = xena.field_codes(
    HUB,
    CLINICAL,
    ["clinical_stage"],
)

if not code_info or not code_info[0].get("code"):
    raise RuntimeError(
        "Could not retrieve categorical codes for clinical_stage."
    )

stage_codes = code_info[0]["code"].split("\t")

def decode_stage(value):

    if value is None:
        return np.nan

    if pd.isna(value):
        return np.nan

    try:
        idx = int(float(value))
        if 0 <= idx < len(stage_codes):
            return stage_codes[idx]
    except (ValueError, TypeError):
        pass

    return str(value)

stage_decoded = stage_raw.apply(
    decode_stage
)

print()
print("Decoded clinical_stage values:")

print(
    stage_decoded
    .value_counts(dropna=False)
    .to_string()
)

# ------------------------------------------------------------------
# Normalize stage
#
# Positive = Stage III/IV
# Negative = Stage I/II
# Unknown = excluded from this endpoint
# ------------------------------------------------------------------

def stage_class(value):

    if value is None:
        return np.nan

    text = str(value).strip().upper()

    if text in {
        "",
        "NAN",
        "NA",
        "N/A",
        "UNKNOWN",
        "NOT REPORTED",
    }:
        return np.nan

    text = text.replace(".", " ")

    # Stage III / IV, including substages such as IIIA, IIIB, IIIC, IVA.
    if re.search(r"\b(?:III|IV)[A-Z0-9]*\b", text):
        return 1

    # Stage I / II, including substages such as IA, IB, IIA, IIB.
    if re.search(r"\b(?:I|II)[A-Z0-9]*\b", text):
        return 0

    return np.nan


df["stage_raw"] = stage_raw
df["stage_label"] = stage_decoded.apply(stage_class)

print()
print("=" * 80)
print("STAGE AUDIT")
print("=" * 80)

print()
print(
    "Stage information observed:",
    int(df["stage_raw"].notna().sum()),
)

print(
    "Stage I/II:",
    int((df["stage_label"] == 0).sum()),
)

print(
    "Stage III/IV:",
    int((df["stage_label"] == 1).sum()),
)

print(
    "Unusable/unknown:",
    int(df["stage_label"].isna().sum()),
)

print()
print("Raw stage examples:")

print(
    df["stage_raw"]
    .dropna()
    .astype(str)
    .value_counts()
    .head(20)
    .to_string()
)

# ------------------------------------------------------------------
# Keep only patients with usable endpoint
# ------------------------------------------------------------------

df = df[
    df["stage_label"].notna()
].copy()

df["stage_label"] = df["stage_label"].astype(int)

print()
print(
    "Patients used for stage model:",
    len(df),
)

# ------------------------------------------------------------------
# Feature layers
# ------------------------------------------------------------------

mirnas = [
    c for c in df.columns
    if c.startswith("hsa-")
]

cnv = [
    c for c in df.columns
    if c.startswith("CNV_")
]

protein = [
    c for c in df.columns
    if c.startswith("PROT_")
]

required_rna = [
    c for c in RNA
    if c in df.columns
]

layer_sets = {
    "age_only": [
        "age",
    ],

    "age_plus_151mirna": [
        "age",
    ] + mirnas,

    "age_mirna_rna": [
        "age",
    ] + mirnas + required_rna,

    "age_mirna_rna_cnv": [
        "age",
    ] + mirnas + required_rna + cnv,

    "full_multimodal": [
        "age",
    ] + mirnas + required_rna + cnv + protein,
}

print()
print("=" * 80)
print("FEATURE LAYERS")
print("=" * 80)

for name, features in layer_sets.items():
    print(
        f"{name:28} "
        f"{len(features):4d} features"
    )

# ------------------------------------------------------------------
# Validate all feature data
# ------------------------------------------------------------------

for name, features in layer_sets.items():

    missing = [
        c for c in features
        if c not in df.columns
    ]

    if missing:
        raise RuntimeError(
            f"{name} missing columns: {missing}"
        )

# ------------------------------------------------------------------
# Nested CV
# ------------------------------------------------------------------

outer_cv = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=42,
)

inner_cv = StratifiedKFold(
    n_splits=4,
    shuffle=True,
    random_state=43,
)

all_predictions = pd.DataFrame({
    "sample_id": df["sample_id"].values,
    "stage_label": df["stage_label"].values,
})

results = {}

for model_name, features in layer_sets.items():

    print()
    print("=" * 80)
    print(model_name.upper())
    print("=" * 80)

    X = df[features].apply(
        pd.to_numeric,
        errors="coerce",
    )

    y = df["stage_label"].astype(int)

    if X.isna().any().any():
        raise RuntimeError(
            f"{model_name} contains missing feature values."
        )

    n_features = X.shape[1]

    # --------------------------------------------------------------
    # Age-only special case
    # --------------------------------------------------------------

    if n_features == 1:

        pipeline = Pipeline([
            (
                "scale",
                StandardScaler(),
            ),
            (
                "clf",
                LogisticRegression(
                    max_iter=5000,
                    class_weight="balanced",
                    solver="liblinear",
                ),
            ),
        ])

        param_grid = {
            "clf__C": [
                0.03,
                0.1,
                0.3,
                1.0,
                3.0,
            ],
        }

    else:

        candidate_k = [
            5,
            10,
            20,
            30,
            50,
            75,
            100,
            150,
            200,
        ]

        candidate_k = [
            k
            for k in candidate_k
            if k <= n_features
        ]

        candidate_k.append("all")

        pipeline = Pipeline([
            (
                "scale",
                StandardScaler(),
            ),
            (
                "select",
                SelectKBest(
                    score_func=f_classif,
                ),
            ),
            (
                "clf",
                LogisticRegression(
                    max_iter=5000,
                    class_weight="balanced",
                    solver="liblinear",
                ),
            ),
        ])

        param_grid = {
            "select__k": candidate_k,
            "clf__C": [
                0.01,
                0.03,
                0.1,
                0.3,
                1.0,
                3.0,
            ],
        }

    oof_prob = np.full(
        len(X),
        np.nan,
    )

    oof_pred = np.full(
        len(X),
        np.nan,
    )

    fold_results = []

    for fold, (train_idx, test_idx) in enumerate(
        outer_cv.split(X, y),
        start=1,
    ):

        search = GridSearchCV(
            estimator=pipeline,
            param_grid=param_grid,
            scoring="roc_auc",
            cv=inner_cv,
            n_jobs=-1,
            refit=True,
        )

        search.fit(
            X.iloc[train_idx],
            y.iloc[train_idx],
        )

        prob = search.predict_proba(
            X.iloc[test_idx]
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

        fold_results.append({
            "fold": fold,
            "test_n": int(len(test_idx)),
            "auc": float(auc),
            "accuracy": float(acc),
            "best_params": search.best_params_,
        })

        print(
            f"Fold {fold}: "
            f"test={len(test_idx)} "
            f"AUC={auc:.6f} "
            f"accuracy={acc:.6f}"
        )

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

    tn, fp, fn, tp = confusion_matrix(
        y,
        oof_pred,
    ).ravel()

    sensitivity = (
        tp / (tp + fn)
        if tp + fn
        else np.nan
    )

    specificity = (
        tn / (tn + fp)
        if tn + fp
        else np.nan
    )

    print()
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

    results[model_name] = {
        "patients": int(len(X)),
        "features": int(n_features),
        "stage_I_II": int((y == 0).sum()),
        "stage_III_IV": int((y == 1).sum()),
        "roc_auc": float(auc),
        "pr_auc": float(pr_auc),
        "accuracy": float(accuracy),
        "sensitivity": float(sensitivity),
        "specificity": float(specificity),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
        "folds": fold_results,
    }

    all_predictions[
        model_name + "_probability"
    ] = oof_prob

    all_predictions[
        model_name + "_prediction"
    ] = oof_pred.astype(int)

# ------------------------------------------------------------------
# Save
# ------------------------------------------------------------------

OUT_PRED.parent.mkdir(
    parents=True,
    exist_ok=True,
)

all_predictions.to_csv(
    OUT_PRED,
    index=False,
)

OUT_RESULTS.write_text(
    json.dumps(
        {
            "endpoint": (
                "TCGA-OV stage III/IV vs stage I/II"
            ),
            "dataset": (
                "TCGA-OV 223-patient multimodal cohort"
            ),
            "models": results,
            "important_exclusions": {
                "methylation": (
                    "excluded: only 9 overlapping patients"
                ),
                "mutations": (
                    "excluded: MC3 sparse extraction "
                    "not yet verified"
                ),
            },
        },
        indent=2,
    )
)

print()
print("=" * 80)
print("ABLATION SUMMARY")
print("=" * 80)

for name, r in results.items():

    print(
        f"{name:28} "
        f"AUC={r['roc_auc']:.6f}"
    )

print()
print("=" * 80)
print("FILES SAVED")
print("=" * 80)

print(OUT_PRED)
print(OUT_RESULTS)

print()
print("IMPORTANT:")
print(
    "This is a TCGA tumor-stage model."
)
print(
    "It is NOT a cancer-vs-control diagnostic model."
)
print(
    "All feature selection occurred inside the"
)
print(
    "cross-validation training folds."
)
