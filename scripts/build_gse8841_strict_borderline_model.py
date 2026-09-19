#!/usr/bin/env python3

from pathlib import Path
import json

import numpy as np
import pandas as pd

from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path.home() / "ovarian-ai"

DATA = ROOT / "data/processed/gse8841_borderline_expression.csv.gz"

OUT_PRED = ROOT / "reports/tables/gse8841_strict_borderline_predictions.csv"

OUT_RESULTS = ROOT / "reports/tables/gse8841_strict_borderline_results.json"


print("=" * 80)
print("GSE8841 STRICT BORDERLINE MODEL")
print("=" * 80)

df = pd.read_csv(DATA)

# ---------------------------------------------------------------
# Restrict to serous + mucinous.
# Histotype is used ONLY for stratification, never as a feature.
# ---------------------------------------------------------------

df = df[
    df["histotype"].str.lower().isin(["serous", "mucinous"])
].copy()

df["histotype"] = df["histotype"].str.lower()

y = (
    df["label"] == "stage_I_cancer"
).astype(int)

# Create strata so each fold preserves BOTH diagnosis and histotype.
strata = (
    df["histotype"]
    + "_"
    + df["label"]
)

metadata = [
    "sample_id",
    "title",
    "label",
    "figo_stage",
    "figo_substage",
    "histotype",
    "grade",
    "age",
    "overall_survival_days",
    "progression_free_survival_days",
    "relapsed",
    "status",
]

expression_cols = [
    c for c in df.columns
    if c not in metadata
]

X = df[expression_cols].apply(
    pd.to_numeric,
    errors="coerce",
)

print()
print("Patients:", len(df))
print("Stage-I cancer:", int((y == 1).sum()))
print("Borderline:", int((y == 0).sum()))
print("Expression features:", X.shape[1])

print()
print("Histotype × diagnosis:")
print(
    pd.crosstab(
        df["histotype"],
        df["label"],
    ).to_string()
)

# ---------------------------------------------------------------
# Remove expression features with missing values.
# ---------------------------------------------------------------

missing_fraction = X.isna().mean()

keep = missing_fraction == 0

X = X.loc[:, keep]

print()
print(
    "Features removed for missingness:",
    int((~keep).sum())
)

print(
    "Features retained:",
    X.shape[1]
)

if X.isna().any().any():
    raise RuntimeError("Missing values remain.")

# ---------------------------------------------------------------
# Pre-specified feature count.
#
# We are NOT searching for the best k after seeing the results.
# k=10 is fixed before evaluation.
# ---------------------------------------------------------------

K = min(10, X.shape[1])

N_REPEATS = 20

repeat_results = []
prediction_rows = []

for repeat in range(1, N_REPEATS + 1):

    cv = StratifiedKFold(
        n_splits=5,
        shuffle=True,
        random_state=1000 + repeat,
    )

    repeat_true = []
    repeat_prob = []
    repeat_pred = []

    for fold, (train_idx, test_idx) in enumerate(
        cv.split(X, strata),
        start=1,
    ):

        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]

        y_train = y.iloc[train_idx]
        y_test = y.iloc[test_idx]

        model = Pipeline([
            (
                "scale",
                StandardScaler(),
            ),
            (
                "select",
                SelectKBest(
                    score_func=f_classif,
                    k=K,
                ),
            ),
            (
                "clf",
                LogisticRegression(
                    max_iter=5000,
                    class_weight="balanced",
                    C=0.1,
                    solver="liblinear",
                ),
            ),
        ])

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

        repeat_true.extend(
            y_test.tolist()
        )

        repeat_prob.extend(
            prob.tolist()
        )

        repeat_pred.extend(
            pred.tolist()
        )

        for idx, p, pr in zip(
            test_idx,
            prob,
            pred,
        ):

            prediction_rows.append({
                "repeat": repeat,
                "fold": fold,
                "sample_id": df.iloc[idx]["sample_id"],
                "true_label": int(y.iloc[idx]),
                "cancer_probability": float(p),
                "predicted_label": int(pr),
            })

    repeat_true = np.asarray(
        repeat_true
    )

    repeat_prob = np.asarray(
        repeat_prob
    )

    repeat_pred = np.asarray(
        repeat_pred
    )

    auc = roc_auc_score(
        repeat_true,
        repeat_prob,
    )

    accuracy = accuracy_score(
        repeat_true,
        repeat_pred,
    )

    balanced_acc = balanced_accuracy_score(
        repeat_true,
        repeat_pred,
    )

    tn, fp, fn, tp = confusion_matrix(
        repeat_true,
        repeat_pred,
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

    repeat_results.append({
        "repeat": repeat,
        "auc": float(auc),
        "accuracy": float(accuracy),
        "balanced_accuracy": float(balanced_acc),
        "sensitivity": float(sensitivity),
        "specificity": float(specificity),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
    })

    print(
        f"Repeat {repeat:02d}/20: "
        f"AUC={auc:.6f} "
        f"balanced_acc={balanced_acc:.6f}"
    )

# ---------------------------------------------------------------
# Summary across complete repeats
# ---------------------------------------------------------------

aucs = np.array([
    r["auc"]
    for r in repeat_results
])

balanced_accs = np.array([
    r["balanced_accuracy"]
    for r in repeat_results
])

sensitivities = np.array([
    r["sensitivity"]
    for r in repeat_results
])

specificities = np.array([
    r["specificity"]
    for r in repeat_results
])

print()
print("=" * 80)
print("FINAL GSE8841 STRICT BORDERLINE RESULTS")
print("=" * 80)

print(
    f"Patients: {len(df)}"
)

print(
    f"Stage-I cancer: {int((y == 1).sum())}"
)

print(
    f"Borderline: {int((y == 0).sum())}"
)

print()
print(
    f"Mean 5-fold AUC:       {aucs.mean():.6f}"
)

print(
    f"SD of AUC:             {aucs.std(ddof=1):.6f}"
)

print(
    f"Median AUC:            {np.median(aucs):.6f}"
)

print(
    f"Min AUC:               {aucs.min():.6f}"
)

print(
    f"Max AUC:               {aucs.max():.6f}"
)

print()
print(
    f"Mean balanced accuracy:{balanced_accs.mean():.6f}"
)

print(
    f"Mean sensitivity:      {sensitivities.mean():.6f}"
)

print(
    f"Mean specificity:      {specificities.mean():.6f}"
)

print()
print("=" * 80)
print("INTERPRETATION")
print("=" * 80)

if aucs.mean() < 0.60:
    print(
        "Weak evidence of reproducible discrimination."
    )
elif aucs.mean() < 0.70:
    print(
        "Modest exploratory discrimination."
    )
elif aucs.mean() < 0.80:
    print(
        "Potentially meaningful exploratory discrimination."
    )
else:
    print(
        "Strong exploratory discrimination; independent "
        "replication is still required."
    )

print()
print("=" * 80)
print("SAVING")
print("=" * 80)

pred_df = pd.DataFrame(
    prediction_rows
)

pred_df.to_csv(
    OUT_PRED,
    index=False,
)

results = {
    "dataset": "GSE8841",
    "subset": "serous + mucinous",
    "patients": int(len(df)),
    "stage_I_cancer": int((y == 1).sum()),
    "borderline": int((y == 0).sum()),
    "expression_features_input": int(len(expression_cols)),
    "expression_features_retained": int(X.shape[1]),
    "features_selected_per_fold": int(K),
    "repeats": int(N_REPEATS),
    "folds_per_repeat": 5,
    "mean_auc": float(aucs.mean()),
    "sd_auc": float(aucs.std(ddof=1)),
    "median_auc": float(np.median(aucs)),
    "min_auc": float(aucs.min()),
    "max_auc": float(aucs.max()),
    "mean_balanced_accuracy": float(
        balanced_accs.mean()
    ),
    "mean_sensitivity": float(
        sensitivities.mean()
    ),
    "mean_specificity": float(
        specificities.mean()
    ),
    "repeat_results": repeat_results,
}

OUT_RESULTS.write_text(
    json.dumps(
        results,
        indent=2,
    )
)

print()
print(OUT_PRED)
print(OUT_RESULTS)

print()
print("NO CLINICAL CLAIM.")
print(
    "This is exploratory tissue-expression research."
)
print(
    "GSE8841 is not serum miRNA validation."
)
