#!/usr/bin/env python3

"""Fit and persist the MethylBERT supplementary-matrix Ridge pipeline."""

import json
import zipfile
from pathlib import Path
from urllib.request import urlopen

import joblib
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import RidgeClassifier
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


PROJECT = Path(__file__).resolve().parents[1]
ARCHIVE = PROJECT / "data" / "external" / "methylbert_supplementary.zip"
DATA = PROJECT / "data" / "external" / "MethylBERT_EOC_DataS1_mmc2.xlsx"
MODEL = PROJECT / "models" / "methylbert_ridge_cancer_vs_healthy.joblib"
METADATA = PROJECT / "models" / "methylbert_ridge_cancer_vs_healthy.json"
URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11384945/supplementaryFiles"


def load_training_data():
    if not ARCHIVE.exists():
        ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        ARCHIVE.write_bytes(urlopen(URL, timeout=180).read())
    if not DATA.exists():
        with zipfile.ZipFile(ARCHIVE) as archive:
            DATA.write_bytes(archive.read("mmc2.xlsx"))
    frame = pd.read_excel(DATA)
    subset = frame[frame["Sample.Type"].isin(["Ovarian cancer", "Healthy female"])].copy()
    features = [
        column for column in frame.columns
        if str(column).startswith("chr")
        and subset[column].notna().any()
        and subset[column].nunique(dropna=True) > 1
    ]
    return subset, features


def main():
    frame, features = load_training_data()
    X = frame[features]
    y = (frame["Sample.Type"] == "Ovarian cancer").astype(int)
    pipeline = Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        ("select", SelectKBest(f_classif)),
        ("model", RidgeClassifier(class_weight="balanced")),
    ])
    search = GridSearchCV(
        pipeline,
        {
            "select__k": [25, 50, 100, 250, 493, "all"],
            "model__alpha": [0.01, 0.1, 1.0, 10.0, 100.0],
        },
        scoring="balanced_accuracy",
        cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=42),
        n_jobs=-1,
        refit=True,
    )
    search.fit(X, y)
    MODEL.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(search.best_estimator_, MODEL)
    METADATA.write_text(json.dumps({
        "endpoint": "MethylBERT-EOC ovarian cancer vs healthy",
        "patients": len(frame),
        "features": features,
        "best_params": search.best_params_,
        "training_only_fit": True,
        "external_refit": False,
    }, indent=2) + "\n")
    print(f"Saved frozen pipeline: {MODEL}")
    print(f"Saved metadata: {METADATA}")
    print(f"Best parameters: {search.best_params_}")


if __name__ == "__main__":
    main()