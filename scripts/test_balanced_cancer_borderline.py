#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd

from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
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

df = pd.read_csv(DATA)

df = df[
    df["group"].isin(["cancer", "borderline"])
].copy()

mirnas = [
    c for c in df.columns
    if c.startswith("hsa-")
]

features = mirnas + ["age"]

X = df[features].apply(pd.to_numeric, errors="coerce")
y = (df["group"] == "cancer").astype(int)

valid = X.notna().all(axis=1)

X = X.loc[valid].reset_index(drop=True)
y = y.loc[valid].reset_index(drop=True)

print("=" * 80)
print("BALANCED CANCER vs BORDERLINE ANALYSIS")
print("=" * 80)

print()
print("Original:")
print("  cancer:", int((y == 1).sum()))
print("  borderline:", int((y == 0).sum()))

# ---------------------------------------------------------------
# Randomly sample equal numbers from each class.
# Repeat many times to see stability.
# ---------------------------------------------------------------

n_borderline = int((y == 0).sum())
rng = np.random.default_rng(42)

N_REPEATS = 20

aucs = []
pr_aucs = []

for repeat in range(N_REPEATS):

    cancer_idx = np.where(
        y.to_numpy() == 1
    )[0]

    borderline_idx = np.where(
        y.to_numpy() == 0
    )[0]

    cancer_sample = rng.choice(
        cancer_idx,
        size=n_borderline,
        replace=False,
    )

    selected = np.concatenate([
        cancer_sample,
        borderline_idx,
    ])

    Xr = X.iloc[selected].reset_index(drop=True)
    yr = y.iloc[selected].reset_index(drop=True)

    cv = StratifiedKFold(
        n_splits=5,
        shuffle=True,
        random_state=100 + repeat,
    )

    oof = np.full(len(Xr), np.nan)

    for train_idx, test_idx in cv.split(Xr, yr):

        model = Pipeline([
            (
                "scale",
                StandardScaler(),
            ),
            (
                "select",
                SelectKBest(
                    f_classif,
                    k=30,
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

        model.fit(
            Xr.iloc[train_idx],
            yr.iloc[train_idx],
        )

        oof[test_idx] = model.predict_proba(
            Xr.iloc[test_idx]
        )[:, 1]

    auc = roc_auc_score(
        yr,
        oof,
    )

    pr = average_precision_score(
        yr,
        oof,
    )

    aucs.append(auc)
    pr_aucs.append(pr)

    print(
        f"Repeat {repeat + 1:02d}: "
        f"AUC={auc:.6f} "
        f"PR={pr:.6f}"
    )

print()
print("=" * 80)
print("SUMMARY")
print("=" * 80)

print(
    f"Mean AUC:   {np.mean(aucs):.6f}"
)

print(
    f"SD AUC:     {np.std(aucs, ddof=1):.6f}"
)

print(
    f"Min AUC:    {np.min(aucs):.6f}"
)

print(
    f"Max AUC:    {np.max(aucs):.6f}"
)

print(
    f"Mean PR AUC:{np.mean(pr_aucs):.6f}"
)

print()
print("This is a stability experiment.")
print("It is NOT external validation.")
