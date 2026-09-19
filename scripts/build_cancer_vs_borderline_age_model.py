#!/usr/bin/env python3

from pathlib import Path
import json

import numpy as np
import pandas as pd

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
    / "gse106817_external_validation.csv.gz"
)

OUT_PRED = (
    ROOT
    / "reports"
    / "tables"
    / "gse106817_cancer_vs_borderline_age_predictions.csv"
)

OUT_RESULTS = (
    ROOT
    / "reports"
    / "tables"
    / "gse106817_cancer_vs_borderline_age_results.json"
)


print("=" * 80)
print("OVARIAN CANCER vs BORDERLINE — 170 miRNAs + AGE")
print("=" * 80)

df = pd.read_csv(DATA)

df = df[
    df["group"].isin([
        "cancer",
        "borderline",
    ])
].copy()

mirna_cols = [
    c for c in df.columns
    if c.startswith("hsa-")
]

feature_cols = mirna_cols + ["age"]

X = df[feature_cols].apply(
    pd.to_numeric,
    errors="coerce",
)

y = (
    df["group"] == "cancer"
).astype(int)

valid = X.notna().all(axis=1)

X = X.loc[valid].reset_index(drop=True)
y = y.loc[valid].reset_index(drop=True)
meta = df.loc[valid].reset_index(drop=True)

print()
print("Patients:", len(X))
print("Cancer:", int(y.sum()))
print("Borderline:", int((y == 0).sum()))
print("Features:", len(feature_cols))
print("Age missing:", int(X["age"].isna().sum()))

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

pipeline = Pipeline([
    (
        "scaler",
        StandardScaler(),
    ),
    (
        "select",
        SelectKBest(
            score_func=f_classif,
        ),
    ),
    (
        "classifier",
        LogisticRegression(
            max_iter=5000,
            class_weight="balanced",
            solver="liblinear",
        ),
    ),
])

param_grid = {
    "select__k": [
        5,
        10,
        20,
        30,
        50,
        75,
        100,
        "all",
    ],
    "classifier__C": [
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

print()
print("=" * 80)
print("NESTED CROSS-VALIDATION")
print("=" * 80)

for fold, (train_idx, test_idx) in enumerate(
    outer_cv.split(X, y),
    start=1,
):

    search = GridSearchCV(
        pipeline,
        param_grid,
        cv=inner_cv,
        scoring="roc_auc",
        n_jobs=-1,
        refit=True,
    )

    search.fit(
        X.iloc[train_idx],
        y.iloc[train_idx],
    )

    probability = search.predict_proba(
        X.iloc[test_idx]
    )[:, 1]

    prediction = (
        probability >= 0.5
    ).astype(int)

    oof_prob[test_idx] = probability
    oof_pred[test_idx] = prediction

    fold_auc = roc_auc_score(
        y.iloc[test_idx],
        probability,
    )

    fold_acc = accuracy_score(
        y.iloc[test_idx],
        prediction,
    )

    fold_results.append({
        "fold": fold,
        "auc": float(fold_auc),
        "accuracy": float(fold_acc),
        "best_params": search.best_params_,
    })

    print(
        f"Fold {fold}: "
        f"test={len(test_idx)} "
        f"AUC={fold_auc:.6f} "
        f"accuracy={fold_acc:.6f} "
        f"best={search.best_params_}"
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

sensitivity = tp / (tp + fn)
specificity = tn / (tn + fp)

print()
print("=" * 80)
print("FINAL RESULTS")
print("=" * 80)

print(f"ROC AUC:     {auc:.6f}")
print(f"PR AUC:      {pr_auc:.6f}")
print(f"Accuracy:    {accuracy:.6f}")
print(f"Sensitivity: {sensitivity:.6f}")
print(f"Specificity: {specificity:.6f}")

print()
print("Confusion matrix:")
print(f"  TN={tn} FP={fp}")
print(f"  FN={fn} TP={tp}")

pred = meta[
    [
        "sample_id",
        "title",
        "group",
        "age",
    ]
].copy()

pred["true_label"] = y.values
pred["cancer_probability"] = oof_prob
pred["predicted_label"] = oof_pred.astype(int)

pred.to_csv(
    OUT_PRED,
    index=False,
)

results = {
    "dataset": "GSE106817",
    "endpoint": "cancer_vs_borderline",
    "features": len(feature_cols),
    "mirnas": len(mirna_cols),
    "age_used": True,
    "patients": len(X),
    "cancer": int(y.sum()),
    "borderline": int((y == 0).sum()),
    "roc_auc": float(auc),
    "pr_auc": float(pr_auc),
    "accuracy": float(accuracy),
    "sensitivity": float(sensitivity),
    "specificity": float(specificity),
    "tn": int(tn),
    "fp": int(fp),
    "fn": int(fn),
    "tp": int(tp),
    "outer_cv": "5-fold stratified",
    "inner_cv": "4-fold stratified",
    "feature_selection_inside_training_folds": True,
    "fold_results": fold_results,
}

OUT_RESULTS.write_text(
    json.dumps(
        results,
        indent=2,
    )
)

print()
print("Saved:")
print(OUT_PRED)
print(OUT_RESULTS)
