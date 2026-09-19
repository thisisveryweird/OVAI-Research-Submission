#!/usr/bin/env python3

"""Evaluate patient-linked TCGA-OV stage models with nested cross-validation.

This is a tumor-stage research endpoint, not a cancer-detection model. The
features come from the public TCGA-OV Xena mirrors used by the existing
project, while the cohort is restricted to samples with RNA, miRNA, and
clinical records. All imputation, feature selection, and hyperparameter
selection happen inside training folds.
"""

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import xenaPython as xena
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
HUB = "https://tcga.xenahubs.net"
RNA_DATASET = "TCGA.OV.sampleMap/HiSeqV2"
MIRNA_DATASET = "TCGA.OV.sampleMap/miRNA_HiSeq_gene"
CLINICAL_DATASET = "TCGA.OV.sampleMap/OV_clinicalMatrix"
OUT_DATA = PROJECT / "data" / "processed" / "tcga_ov_stage_features.csv.gz"
OUT_PRED = PROJECT / "reports" / "tables" / "tcga_ov_stage_oof_predictions.csv"
OUT_RESULTS = PROJECT / "reports" / "tables" / "tcga_ov_stage_results.json"

RNA_GENES = [
    "PAX8", "EPCAM", "KRT7", "KRT8", "KRT18",
    "MUC16", "WFDC2", "FOLR1", "WT1", "MKI67",
]


def patient_id(sample_id):
    parts = str(sample_id).split("-")
    return "-".join(parts[:3]) if len(parts) >= 3 else str(sample_id)


def decode_stage(value, codes):
    if value is None or pd.isna(value):
        return np.nan
    try:
        index = int(float(value))
        if 0 <= index < len(codes):
            value = codes[index]
    except (TypeError, ValueError):
        pass
    text = str(value).strip().upper().replace(".", " ")
    if not text or text in {"NAN", "NA", "UNKNOWN", "NOT REPORTED"}:
        return np.nan
    if re.search(r"\b(?:III|IV)[A-Z0-9]*\b", text):
        return 1
    if re.search(r"\b(?:I|II)[A-Z0-9]*\b", text):
        return 0
    return np.nan


def query_matrix(dataset, samples, fields):
    _, values = xena.dataset_probe_values(HUB, dataset, samples, fields)
    return np.asarray(values).T


def load_features():
    rna_samples = [str(value) for value in xena.dataset_samples(HUB, RNA_DATASET, None)]
    mirna_samples = [str(value) for value in xena.dataset_samples(HUB, MIRNA_DATASET, None)]
    clinical_samples = {str(value) for value in xena.dataset_samples(HUB, CLINICAL_DATASET, None)}
    samples = [
        sample for sample in rna_samples
        if sample in set(mirna_samples) and sample in clinical_samples
    ]

    # Prevent repeated aliquots from creating duplicate patient rows.
    samples = list(dict.fromkeys(samples))
    by_patient = {}
    for sample in samples:
        by_patient.setdefault(patient_id(sample), sample)
    samples = list(by_patient.values())

    rna_fields = set(str(value) for value in xena.dataset_field(HUB, RNA_DATASET))
    missing_genes = [gene for gene in RNA_GENES if gene not in rna_fields]
    if missing_genes:
        raise RuntimeError(f"Missing RNA fields: {missing_genes}")

    mirna_fields = [str(value) for value in xena.dataset_field(HUB, MIRNA_DATASET)]
    rna_values = query_matrix(RNA_DATASET, samples, RNA_GENES)
    mirna_values = query_matrix(MIRNA_DATASET, samples, mirna_fields)
    frame = pd.DataFrame(rna_values, columns=RNA_GENES)
    frame = pd.concat([
        frame,
        pd.DataFrame(mirna_values, columns=[f"mirna_{field}" for field in mirna_fields]),
    ], axis=1)
    frame.insert(0, "sample_id", samples)
    frame.insert(1, "patient_id", [patient_id(sample) for sample in samples])

    clinical_fields = set(str(value) for value in xena.dataset_field(HUB, CLINICAL_DATASET))
    if "clinical_stage" not in clinical_fields:
        raise RuntimeError("clinical_stage is not available in the TCGA clinical matrix")
    _, stage_values = xena.dataset_probe_values(
        HUB, CLINICAL_DATASET, samples, ["clinical_stage"]
    )
    code_info = xena.field_codes(HUB, CLINICAL_DATASET, ["clinical_stage"])
    codes = code_info[0].get("code", "").split("\t") if code_info else []
    frame["stage_label"] = [decode_stage(value, codes) for value in stage_values[0]]

    age_field = "age_at_initial_pathologic_diagnosis"
    if age_field in clinical_fields:
        _, age_values = xena.dataset_probe_values(HUB, CLINICAL_DATASET, samples, [age_field])
        frame["age"] = pd.to_numeric(age_values[0], errors="coerce")
    frame = frame[frame["stage_label"].notna()].copy()
    frame["stage_label"] = frame["stage_label"].astype(int)
    metadata = {"sample_id", "patient_id", "stage_label"}
    usable = [
        column for column in frame.columns
        if column in metadata
        or (
            frame[column].notna().mean() >= 0.80
            and frame[column].nunique(dropna=True) > 1
        )
    ]
    frame = frame[usable]
    return frame


def model_pipeline(feature_count):
    candidates = [k for k in [5, 10, 20, 30, 50, 75, 100, 150, 200] if k <= feature_count]
    candidates.append("all")
    return Pipeline([
        ("impute", SimpleImputer(strategy="median", add_indicator=True)),
        ("scale", StandardScaler()),
        ("select", SelectKBest(score_func=f_classif)),
        ("classifier", LogisticRegression(
            max_iter=5000, class_weight="balanced", solver="liblinear"
        )),
    ]), {
        "select__k": candidates,
        "classifier__C": [0.01, 0.03, 0.1, 0.3, 1.0, 3.0],
    }


def main():
    print("Loading patient-linked TCGA-OV features...")
    frame = load_features()
    OUT_DATA.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT_DATA, index=False, compression="gzip")
    print(f"Patients with usable stage: {len(frame)}")
    print(frame["stage_label"].value_counts().sort_index().to_string())

    mirna = [column for column in frame if column.startswith("mirna_")]
    layers = {
        "age_only": ["age"],
        "rna_only": [column for column in RNA_GENES if column in frame],
        "mirna_only": mirna,
        "age_rna_mirna": ["age"] + RNA_GENES + mirna,
    }
    outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    inner = StratifiedKFold(n_splits=3, shuffle=True, random_state=43)
    predictions = frame[["patient_id", "sample_id", "stage_label"]].copy()
    results = {}

    for name, columns in layers.items():
        if any(column not in frame for column in columns):
            raise RuntimeError(f"{name} has missing columns")
        X = frame[columns].apply(pd.to_numeric, errors="coerce")
        y = frame["stage_label"]
        probabilities = np.full(len(frame), np.nan)
        predicted = np.full(len(frame), np.nan)
        folds = []
        for fold, (train_idx, test_idx) in enumerate(outer.split(X, y), start=1):
            pipeline, grid = model_pipeline(X.shape[1])
            search = GridSearchCV(
                pipeline,
                grid,
                scoring="balanced_accuracy",
                cv=inner,
                n_jobs=-1,
                refit=True,
            )
            search.fit(X.iloc[train_idx], y.iloc[train_idx])
            probability = search.predict_proba(X.iloc[test_idx])[:, 1]
            prediction = (probability >= 0.5).astype(int)
            probabilities[test_idx] = probability
            predicted[test_idx] = prediction
            folds.append({
                "fold": fold,
                "accuracy": float(accuracy_score(y.iloc[test_idx], prediction)),
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
            "sensitivity": float(tp / (tp + fn)) if tp + fn else None,
            "specificity": float(tn / (tn + fp)) if tn + fp else None,
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
        "endpoint": "TCGA-OV stage III/IV vs stage I/II",
        "source": "TCGA-OV Xena public mirrors",
        "patient_linking": "TCGA barcode first three components",
        "evaluation": "nested 5-fold out-of-fold predictions",
        "class_counts": {
            "stage_I_II": int((frame["stage_label"] == 0).sum()),
            "stage_III_IV": int((frame["stage_label"] == 1).sum()),
        },
        "majority_class_accuracy": float(frame["stage_label"].mean()),
        "models": results,
    }, indent=2) + "\n")
    print(f"Results: {OUT_RESULTS}")


if __name__ == "__main__":
    main()