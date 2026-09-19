#!/usr/bin/env python3

"""Evaluate GSE94533 serum-miRNA cancer detection without cohort leakage."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


PROJECT = Path(__file__).resolve().parents[1]
MATRIX = PROJECT / "data" / "processed" / "gse94533_mirna_tpm.csv.gz"
METADATA = PROJECT / "data" / "processed" / "gse94533_metadata.csv"
CLINICAL = PROJECT / "data" / "processed" / "gse94533_14mirna_clinical.csv.gz"
OUT_PRED = PROJECT / "reports" / "tables" / "gse94533_diagnostic_oof_predictions.csv"
OUT_RESULTS = PROJECT / "reports" / "tables" / "gse94533_diagnostic_results.json"

HISTORICAL_MIRNAS = [
    "hsa-miR-23b-3p", "hsa-miR-29a-3p", "hsa-miR-32-5p",
    "hsa-miR-92a-3p", "hsa-miR-150-5p", "hsa-miR-200a-3p",
    "hsa-miR-200c-3p", "hsa-miR-203a", "hsa-miR-320c",
    "hsa-miR-320d", "hsa-miR-335-5p", "hsa-miR-450b-5p",
    "hsa-miR-1246", "hsa-miR-1307-5p",
]


def build_pipeline(feature_count):
    ks = [k for k in [5, 10, 14, 20, 30, 50, 75, 100] if k <= feature_count]
    ks.append("all")
    return Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        ("select", SelectKBest(f_classif)),
        ("model", LogisticRegression(max_iter=5000, solver="liblinear")),
    ]), {
        "select__k": ks,
        "model__C": [0.01, 0.03, 0.1, 0.3, 1.0, 3.0],
        "model__class_weight": [None, "balanced"],
    }


def main():
    matrix = pd.read_csv(MATRIX)
    metadata = pd.read_csv(METADATA)
    frame = metadata[["sample_id_normalized", "binary_label"]].merge(
        matrix, on="sample_id_normalized", how="inner", validate="one_to_one"
    )
    feature_columns = [column for column in matrix if column.startswith("hsa-")]
    feature_columns = [
        column for column in feature_columns
        if frame[column].notna().any() and frame[column].nunique(dropna=True) > 1
    ]
    frame[feature_columns] = np.log1p(frame[feature_columns].astype(float))
    y = frame["binary_label"].astype(int)
    layers = {
        "historical_14mirna": [column for column in HISTORICAL_MIRNAS if column in frame],
        "full_mirna": feature_columns,
    }
    if CLINICAL.exists():
        clinical = pd.read_csv(CLINICAL)[["sample_id_normalized", "age", "ca125"]]
        frame = frame.merge(clinical, on="sample_id_normalized", how="left", validate="one_to_one")
        layers["historical_14mirna_age_ca125"] = [
            column for column in HISTORICAL_MIRNAS if column in frame
        ] + ["age", "ca125"]
    outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    inner = StratifiedKFold(n_splits=3, shuffle=True, random_state=43)
    predictions = frame[["sample_id_normalized", "binary_label"]].copy()
    results = {}

    for name, columns in layers.items():
        X = frame[columns]
        probabilities = np.full(len(frame), np.nan)
        predicted = np.full(len(frame), np.nan)
        folds = []
        for fold, (train_idx, test_idx) in enumerate(outer.split(X, y), start=1):
            pipeline, grid = build_pipeline(X.shape[1])
            search = GridSearchCV(
                pipeline, grid, scoring="balanced_accuracy", cv=inner,
                n_jobs=-1, refit=True,
            )
            search.fit(X.iloc[train_idx], y.iloc[train_idx])
            probability = search.predict_proba(X.iloc[test_idx])[:, 1]
            prediction = (probability >= 0.5).astype(int)
            probabilities[test_idx] = probability
            predicted[test_idx] = prediction
            folds.append({
                "fold": fold,
                "accuracy": float(accuracy_score(y.iloc[test_idx], prediction)),
                "balanced_accuracy": float(balanced_accuracy_score(y.iloc[test_idx], prediction)),
                "roc_auc": float(roc_auc_score(y.iloc[test_idx], probability)),
            })
        tn, fp, fn, tp = confusion_matrix(y, predicted.astype(int)).ravel()
        results[name] = {
            "patients": int(len(frame)),
            "features": int(len(columns)),
            "accuracy": float(accuracy_score(y, predicted)),
            "balanced_accuracy": float(balanced_accuracy_score(y, predicted)),
            "roc_auc": float(roc_auc_score(y, probabilities)),
            "pr_auc": float(average_precision_score(y, probabilities)),
            "sensitivity": float(tp / (tp + fn)),
            "specificity": float(tn / (tn + fp)),
            "confusion_matrix": {"tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp)},
            "folds": folds,
        }
        predictions[f"{name}_probability"] = probabilities
        predictions[f"{name}_prediction"] = predicted.astype(int)
        print(
            f"{name}: accuracy={results[name]['accuracy']:.4f} "
            f"balanced_accuracy={results[name]['balanced_accuracy']:.4f} "
            f"ROC-AUC={results[name]['roc_auc']:.4f}"
        )

    OUT_PRED.parent.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(OUT_PRED, index=False)
    OUT_RESULTS.write_text(json.dumps({
        "endpoint": "GSE94533 invasive ovarian cancer vs non-cancer",
        "patients": int(len(frame)),
        "class_counts": {"non_cancer": int((y == 0).sum()), "cancer": int((y == 1).sum())},
        "majority_class_accuracy": float(y.mean()),
        "evaluation": "nested 5-fold out-of-fold predictions",
        "models": results,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()