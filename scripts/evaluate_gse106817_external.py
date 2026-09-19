#!/usr/bin/env python3

"""Evaluate the frozen GSE94533 miRNA model on independent serum GSE106817."""

import csv
import gzip
import io
import json
import re
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import QuantileTransformer, StandardScaler


PROJECT = Path(__file__).resolve().parents[1]
TRAIN = PROJECT / "data" / "processed" / "gse94533_14mirna_clinical.csv.gz"
RAW = PROJECT / "data" / "external" / "gse106817_series_matrix.txt.gz"
PLATFORM = PROJECT / "data" / "external" / "GPL21263.txt"
OUT = PROJECT / "reports" / "tables" / "gse106817_external_results.json"

MIRNAS = [
    "hsa-miR-23b-3p", "hsa-miR-29a-3p", "hsa-miR-32-5p", "hsa-miR-92a-3p",
    "hsa-miR-150-5p", "hsa-miR-200a-3p", "hsa-miR-200c-3p", "hsa-miR-203a",
    "hsa-miR-320c", "hsa-miR-320d", "hsa-miR-335-5p", "hsa-miR-450b-5p",
    "hsa-miR-1246", "hsa-miR-1307-5p",
]
MIMAT_BY_NAME = {
    "hsa-miR-23b-3p": "MIMAT0000418", "hsa-miR-29a-3p": "MIMAT0000086",
    "hsa-miR-32-5p": "MIMAT0000090", "hsa-miR-92a-3p": "MIMAT0000093",
    "hsa-miR-150-5p": "MIMAT0000451", "hsa-miR-200a-3p": "MIMAT0000682",
    "hsa-miR-200c-3p": "MIMAT0000617", "hsa-miR-203a": "MIMAT0000264",
    "hsa-miR-320c": "MIMAT0005791", "hsa-miR-320d": "MIMAT0006764",
    "hsa-miR-335-5p": "MIMAT0000765", "hsa-miR-450b-5p": "MIMAT0004784",
    "hsa-miR-1246": "MIMAT0005898", "hsa-miR-1307-5p": "MIMAT0005953",
}


def download(url, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_bytes(urlopen(url, timeout=180).read())


def parse_series(path):
    metadata = {}
    lines = gzip.open(path, "rt", encoding="utf-8", errors="replace").readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_begin"))
    for line in lines[:start]:
        if line.startswith("!Sample_"):
            parts = next(csv.reader([line.rstrip("\n")], delimiter="\t", quotechar='"'))
            metadata[parts[0]] = parts[1:]
    table = pd.read_csv(io.StringIO("".join(lines[start + 1:])), sep="\t", comment="!", quotechar='"')
    table = table.rename(columns={table.columns[0]: "feature_id"})
    table = table[table.feature_id.notna() & (table.feature_id != "!series_matrix_table_end")]
    return metadata, table


def platform_names(path):
    names = {}
    started = False
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("!platform_table_begin"):
            started = True
            continue
        if started and line.startswith("!platform_table_end"):
            break
        if started:
            parts = re.split(r"\s+", line.strip())
            if len(parts) >= 3:
                names[parts[0]] = parts[2]
    return names


def main():
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106817/matrix/GSE106817_series_matrix.txt.gz",
        RAW,
    )
    download(
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc=GPL21263&form=text&view=full",
        PLATFORM,
    )
    metadata, table = parse_series(RAW)
    names = platform_names(PLATFORM)
    feature_by_name = {name: feature for feature, name in names.items()}
    if "hsa-miR-203a-3p" in feature_by_name:
        feature_by_name["hsa-miR-203a"] = feature_by_name["hsa-miR-203a-3p"]
    missing = [name for name in MIRNAS if name not in feature_by_name]
    if missing:
        raise RuntimeError(f"Missing external platform features: {missing}")
    sample_ids = metadata["!Sample_geo_accession"]
    titles = metadata["!Sample_title"]
    values = table.set_index("feature_id").loc[
        [feature_by_name[name] for name in MIRNAS]
    ].T
    values.index = sample_ids
    values.columns = MIRNAS
    labels = []
    keep = []
    for sample_id, title in zip(sample_ids, titles):
        if re.match(r"^Ovarian Cancer\b", title, re.I):
            keep.append(sample_id)
            labels.append(1)
        elif re.match(r"^(?:Non-Cancer|Benign Ovarian Disease)\b", title, re.I):
            keep.append(sample_id)
            labels.append(0)
    external = values.loc[keep].apply(pd.to_numeric, errors="coerce")
    y_external = np.asarray(labels, dtype=int)

    train = pd.read_csv(TRAIN)
    X_train = np.log1p(train[MIRNAS].clip(lower=0))
    y_train = train["binary_label"].astype(int)
    grids = {
        "standard_scale": (StandardScaler(), "standard_scale"),
        "train_fitted_quantile": (
            QuantileTransformer(
                output_distribution="normal",
                n_quantiles=min(100, len(X_train)),
                random_state=42,
            ),
            "train_fitted_quantile",
        ),
    }
    results = {
        "training_cohort": "GSE94533",
        "external_cohort": "GSE106817",
        "endpoint": "ovarian cancer vs non-cancer or benign ovarian disease",
        "patients": int(len(external)),
        "cancer": int(y_external.sum()),
        "non_cancer": int((y_external == 0).sum()),
        "features": MIRNAS,
        "platform": "GPL21263 3D-Gene Human miRNA V21",
        "combat": "not applied: external test must remain untouched",
    }
    for name, (transformer, normalization) in grids.items():
        pipeline = Pipeline([
            ("impute", SimpleImputer(strategy="median")),
            ("normalize", transformer),
            ("select", SelectKBest(f_classif)),
            ("model", LogisticRegression(max_iter=5000, solver="liblinear")),
        ])
        search = GridSearchCV(
            pipeline,
            {
                "select__k": [5, 10, 14, "all"],
                "model__C": [0.03, 0.1, 0.3, 1.0],
                "model__class_weight": [None, "balanced"],
            },
            scoring="balanced_accuracy",
            cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=42),
            refit=True, n_jobs=-1,
        )
        search.fit(X_train, y_train)
        probability = search.predict_proba(np.log1p(external.clip(lower=0)))[:, 1]
        train_probability = cross_val_predict(
            search.best_estimator_, X_train, y_train,
            cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=43),
            method="predict_proba", n_jobs=-1,
        )[:, 1]
        thresholds = np.linspace(0.10, 0.90, 81)
        threshold = max(
            thresholds,
            key=lambda candidate: balanced_accuracy_score(
                y_train, (train_probability >= candidate).astype(int)
            ),
        )
        prediction = (probability >= threshold).astype(int)
        tn, fp, fn, tp = confusion_matrix(y_external, prediction).ravel()
        results[name] = {
            "normalization": normalization,
            "best_params": search.best_params_,
            "frozen_threshold": float(threshold),
            "accuracy": float(accuracy_score(y_external, prediction)),
            "balanced_accuracy": float(balanced_accuracy_score(y_external, prediction)),
            "roc_auc": float(roc_auc_score(y_external, probability)),
            "sensitivity": float(tp / (tp + fn)),
            "specificity": float(tn / (tn + fp)),
            "confusion_matrix": {"tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp)},
        }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()