#!/usr/bin/env python3

"""Evaluate the public MethylBERT-EOC supplementary methylation matrix."""

import io
import json
import zipfile
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import RidgeClassifier
from sklearn.metrics import accuracy_score, balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


PROJECT = Path(__file__).resolve().parents[1]
ARCHIVE = PROJECT / "data" / "external" / "methylbert_supplementary.zip"
DATA = PROJECT / "data" / "external" / "MethylBERT_EOC_DataS1_mmc2.xlsx"
OUT = PROJECT / "reports" / "tables" / "methylbert_methylation_results.json"
URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11384945/supplementaryFiles"


def load_data():
    ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
    if not ARCHIVE.exists():
        ARCHIVE.write_bytes(urlopen(URL, timeout=180).read())
    if not DATA.exists():
        with zipfile.ZipFile(ARCHIVE) as archive:
            DATA.write_bytes(archive.read("mmc2.xlsx"))
    frame = pd.read_excel(DATA)
    feature_columns = [
        column for column in frame.columns
        if str(column).startswith(("chr", "chrom"))
    ]
    if not feature_columns:
        raise RuntimeError("Could not identify methylation feature columns.")
    frame[feature_columns] = frame[feature_columns].apply(pd.to_numeric, errors="coerce")
    feature_columns = [
        column for column in feature_columns
        if frame[column].notna().any() and frame[column].nunique(dropna=True) > 1
    ]
    return frame, feature_columns


def evaluate(frame, features, endpoint):
    if endpoint == "cancer_vs_healthy":
        subset = frame[frame["Sample.Type"].isin(["Ovarian cancer", "Healthy female"])].copy()
        y = (subset["Sample.Type"] == "Ovarian cancer").astype(int).to_numpy()
    else:
        subset = frame[frame["Stage"].isin(["Stage I/II", "Healthy"])].copy()
        y = (subset["Stage"] == "Stage I/II").astype(int).to_numpy()
    X = subset[features]
    outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    inner = StratifiedKFold(n_splits=3, shuffle=True, random_state=43)
    probabilities = np.full(len(subset), np.nan)
    folds = []
    k_values = [k for k in [25, 50, 100, 250, 500] if k <= len(features)] + ["all"]
    for fold, (train_idx, test_idx) in enumerate(outer.split(X, y), start=1):
        pipeline = Pipeline([
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            ("select", SelectKBest(f_classif)),
            ("model", RidgeClassifier(class_weight="balanced")),
        ])
        search = GridSearchCV(
            pipeline,
            {"select__k": k_values, "model__alpha": [0.01, 0.1, 1.0, 10.0, 100.0]},
            scoring="balanced_accuracy", cv=inner, n_jobs=-1, refit=True,
        )
        search.fit(X.iloc[train_idx], y[train_idx])
        score = search.decision_function(X.iloc[test_idx])
        probabilities[test_idx] = score
        prediction = (score >= 0).astype(int)
        folds.append({
            "fold": fold,
            "accuracy": float(accuracy_score(y[test_idx], prediction)),
            "balanced_accuracy": float(balanced_accuracy_score(y[test_idx], prediction)),
            "roc_auc": float(roc_auc_score(y[test_idx], score)),
        })
    prediction = (probabilities >= 0).astype(int)
    return {
        "patients": int(len(subset)),
        "positive": int(y.sum()),
        "negative": int((y == 0).sum()),
        "features_available": int(len(features)),
        "accuracy": float(accuracy_score(y, prediction)),
        "balanced_accuracy": float(balanced_accuracy_score(y, prediction)),
        "roc_auc": float(roc_auc_score(y, probabilities)),
        "folds": folds,
    }


def main():
    frame, features = load_data()
    results = {
        "source": "MethylBERT-EOC supplementary Dataset 1",
        "accession": "HRA011451 / PMC11384945 supplementary matrix",
        "raw_bam_used": False,
        "note": "Public supplement is a processed methylation matrix; no BAM/FASTQ was used.",
        "models": {
            "cancer_vs_healthy": evaluate(frame, features, "cancer_vs_healthy"),
            "early_stage_vs_healthy": evaluate(frame, features, "early_stage_vs_healthy"),
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2) + "\n")
    for name, result in results["models"].items():
        print(name, result["accuracy"], result["balanced_accuracy"], result["roc_auc"])


if __name__ == "__main__":
    main()