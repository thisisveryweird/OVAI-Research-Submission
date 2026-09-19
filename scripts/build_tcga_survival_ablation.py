#!/usr/bin/env python3

from pathlib import Path
import json
import warnings

import numpy as np
import pandas as pd
import xenaPython as xena

from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.model_selection import StratifiedKFold
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
    / "tcga_survival_ablation_oof_predictions.csv"
)

OUT_RESULTS = (
    ROOT
    / "reports"
    / "tables"
    / "tcga_survival_ablation_results.json"
)

HUB = "https://tcga.xenahubs.net"
SURV = "survival/OV_survival.txt"


print("=" * 80)
print("TCGA-OV MULTIMODAL OVERALL-SURVIVAL ABLATION")
print("=" * 80)

df = pd.read_csv(DATA)
df["sample_id"] = df["sample_id"].astype(str)

samples = df["sample_id"].tolist()

# ------------------------------------------------------------------
# Load survival
# ------------------------------------------------------------------

_, values = xena.dataset_probe_values(
    HUB,
    SURV,
    samples,
    ["OS", "OS.time"],
)

survival = pd.DataFrame({
    "sample_id": samples,
    "OS": pd.to_numeric(values[0], errors="coerce"),
    "OS.time": pd.to_numeric(values[1], errors="coerce"),
})

df = df.merge(
    survival,
    on="sample_id",
    how="inner",
    validate="one_to_one",
)

# Only patients with usable OS data.
df = df[
    df["OS"].notna()
    & df["OS.time"].notna()
    & (df["OS.time"] > 0)
].copy()

df["OS"] = df["OS"].astype(int)

print()
print("Usable patients:", len(df))
print("Events/deaths:", int(df["OS"].sum()))
print("Censored:", int((df["OS"] == 0).sum()))

# ------------------------------------------------------------------
# Feature layers
# ------------------------------------------------------------------

mirnas = [
    c for c in df.columns
    if c.startswith("hsa-")
]

rna = [
    c for c in [
        "PAX8", "EPCAM", "KRT7", "KRT8", "KRT18",
        "MUC16", "WFDC2", "FOLR1", "WT1", "MKI67"
    ]
    if c in df.columns
]

cnv = [
    c for c in df.columns
    if c.startswith("CNV_")
]

protein = [
    c for c in df.columns
    if c.startswith("PROT_")
]

layers = {
    "age_only": [
        "age"
    ],

    "age_plus_mirna": [
        "age"
    ] + mirnas,

    "age_mirna_rna": [
        "age"
    ] + mirnas + rna,

    "age_mirna_rna_cnv": [
        "age"
    ] + mirnas + rna + cnv,

    "full_multimodal": [
        "age"
    ] + mirnas + rna + cnv + protein,
}

print()
print("=" * 80)
print("FEATURE LAYERS")
print("=" * 80)

for name, features in layers.items():
    print(f"{name:28} {len(features):4d} features")

# ------------------------------------------------------------------
# Nested CV
# ------------------------------------------------------------------

outer = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=42,
)

inner = StratifiedKFold(
    n_splits=3,
    shuffle=True,
    random_state=43,
)

penalties = [
    0.001,
    0.003,
    0.01,
    0.03,
    0.1,
    0.3,
    1.0,
]

l1_ratios = [
    0.0,
    0.5,
    1.0,
]

predictions = pd.DataFrame({
    "sample_id": df["sample_id"].values,
    "OS": df["OS"].values,
    "OS.time": df["OS.time"].values,
})

results = {}

for model_name, features in layers.items():

    print()
    print("=" * 80)
    print(model_name.upper())
    print("=" * 80)

    X = df[features].apply(
        pd.to_numeric,
        errors="coerce",
    )

    if X.isna().any().any():
        raise RuntimeError(
            f"{model_name} contains missing feature values."
        )

    durations = df["OS"].astype(float).to_numpy()
    events = df["OS"].astype(int).to_numpy()

    oof_risk = np.full(len(df), np.nan)

    fold_results = []

    for fold, (train_idx, test_idx) in enumerate(
        outer.split(X, events),
        start=1,
    ):

        X_train = X.iloc[train_idx].copy()
        X_test = X.iloc[test_idx].copy()

        y_train = df.iloc[train_idx]
        y_test = df.iloc[test_idx]

        scaler = StandardScaler()

        X_train_scaled = pd.DataFrame(
            scaler.fit_transform(X_train),
            columns=features,
            index=X_train.index,
        )

        X_test_scaled = pd.DataFrame(
            scaler.transform(X_test),
            columns=features,
            index=X_test.index,
        )

        # ----------------------------------------------------------
        # Inner validation to choose penalization.
        # ----------------------------------------------------------

        best_score = -np.inf
        best_penalty = None
        best_l1 = None

        for penalty in penalties:

            for l1_ratio in l1_ratios:

                inner_scores = []

                for inner_train_idx, inner_test_idx in inner.split(
                    X_train_scaled,
                    y_train["OS"].to_numpy(),
                ):

                    inner_train = X_train_scaled.iloc[
                        inner_train_idx
                    ].copy()

                    inner_test = X_train_scaled.iloc[
                        inner_test_idx
                    ].copy()

                    inner_train["OS.time"] = y_train.iloc[
                        inner_train_idx
                    ]["OS.time"].to_numpy()

                    inner_train["OS"] = y_train.iloc[
                        inner_train_idx
                    ]["OS"].to_numpy()

                    model = CoxPHFitter(
                        penalizer=penalty,
                        l1_ratio=l1_ratio,
                    )

                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        model.fit(
                            inner_train,
                            duration_col="OS.time",
                            event_col="OS",
                        )

                    score = concordance_index(
                        y_train.iloc[inner_test_idx]["OS.time"],
                        -model.predict_partial_hazard(
                            inner_test
                        ).to_numpy().ravel(),
                        y_train.iloc[inner_test_idx]["OS"],
                    )

                    inner_scores.append(score)

                mean_score = float(np.mean(inner_scores))

                if mean_score > best_score:
                    best_score = mean_score
                    best_penalty = penalty
                    best_l1 = l1_ratio

        # ----------------------------------------------------------
        # Fit chosen model on entire outer training fold.
        # ----------------------------------------------------------

        train_final = X_train_scaled.copy()
        train_final["OS.time"] = y_train["OS.time"].to_numpy()
        train_final["OS"] = y_train["OS"].to_numpy()

        final_model = CoxPHFitter(
            penalizer=best_penalty,
            l1_ratio=best_l1,
        )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            final_model.fit(
                train_final,
                duration_col="OS.time",
                event_col="OS",
            )

        risk = final_model.predict_partial_hazard(
            X_test_scaled
        ).to_numpy().ravel()

        oof_risk[test_idx] = risk

        fold_cindex = concordance_index(
            y_test["OS.time"],
            -risk,
            y_test["OS"],
        )

        fold_results.append({
            "fold": fold,
            "test_n": int(len(test_idx)),
            "c_index": float(fold_cindex),
            "best_penalizer": float(best_penalty),
            "best_l1_ratio": float(best_l1),
            "inner_c_index": float(best_score),
        })

        print(
            f"Fold {fold}: "
            f"n={len(test_idx)} "
            f"C-index={fold_cindex:.6f} "
            f"penalizer={best_penalty} "
            f"l1={best_l1}"
        )

    overall_cindex = concordance_index(
        df["OS.time"],
        -oof_risk,
        df["OS"],
    )

    results[model_name] = {
        "patients": int(len(df)),
        "events": int(df["OS"].sum()),
        "features": int(len(features)),
        "oof_c_index": float(overall_cindex),
        "folds": fold_results,
    }

    predictions[
        model_name + "_risk"
    ] = oof_risk

    print()
    print(
        f"OOF C-index: {overall_cindex:.6f}"
    )

# ------------------------------------------------------------------
# Save
# ------------------------------------------------------------------

OUT_PRED.parent.mkdir(
    parents=True,
    exist_ok=True,
)

predictions.to_csv(
    OUT_PRED,
    index=False,
)

OUT_RESULTS.write_text(
    json.dumps(
        {
            "dataset": "TCGA-OV",
            "endpoint": "overall_survival",
            "usable_patients": int(len(df)),
            "events": int(df["OS"].sum()),
            "censored": int((df["OS"] == 0).sum()),
            "evaluation": (
                "nested 5-fold outer CV / "
                "3-fold inner CV"
            ),
            "results": results,
        },
        indent=2,
    )
)

print()
print("=" * 80)
print("ABLATION SUMMARY")
print("=" * 80)

for name, result in results.items():
    print(
        f"{name:28} "
        f"C-index={result['oof_c_index']:.6f}"
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
    "These are cross-validated exploratory survival results."
)
print(
    "They are not clinical validation."
)
