#!/usr/bin/env python3
"""Compare high-capacity models against the frozen local logistic baseline."""

from __future__ import annotations

import argparse
import logging
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

warnings.filterwarnings("ignore", category=FutureWarning)

PROJECT = Path(__file__).resolve().parents[1]
CLINICAL_DATA = PROJECT / "data" / "processed" / "gse94533_14mirna_clinical.csv.gz"
BASE_MATRIX = PROJECT / "data" / "processed" / "gse94533_mirna_tpm.csv.gz"
FALLBACK_REPORT = PROJECT / "reports" / "error_analysis_predictions.csv"
REPORT = PROJECT / "reports" / "high_capacity_model_results.md"
SEED = 20260919
SIGNATURE = (
    "hsa-miR-200b-3p",
    "hsa-miR-296-3p",
    "hsa-miR-6826-5p",
    "hsa-miR-32-3p",
    "hsa-miR-1255a",
)
BASE_FEATURES = (*SIGNATURE, "age", "ca125")
LOGGER = logging.getLogger("ovai.high_capacity_models")


def validate_endpoint(frame: pd.DataFrame) -> pd.DataFrame:
    required = set(BASE_FEATURES) | {"sample_id_normalized", "diagnosis", "stage", "binary_label"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"Training source missing columns: {missing}")
    keep = frame["diagnosis"].eq("Cancer") & frame["stage"].eq("Stage I/II")
    keep |= frame["diagnosis"].isin(["Benign", "Control"])
    frame = frame.loc[keep].drop_duplicates("sample_id_normalized").reset_index(drop=True)
    labels = frame["binary_label"].astype(int)
    if len(frame) != 113 or labels.sum() != 52 or (labels == 0).sum() != 61:
        raise ValueError(f"Unexpected endpoint composition: n={len(frame)}, positives={int(labels.sum())}, negatives={int((labels == 0).sum())}")
    return frame


def load_training_data() -> tuple[pd.DataFrame, np.ndarray]:
    if CLINICAL_DATA.exists():
        frame = pd.read_csv(CLINICAL_DATA)
        source = CLINICAL_DATA
    else:
        if not FALLBACK_REPORT.exists() or not BASE_MATRIX.exists():
            raise FileNotFoundError(f"Missing clinical fusion and fallback inputs: {CLINICAL_DATA}, {FALLBACK_REPORT}, {BASE_MATRIX}")
        clinical = pd.read_csv(FALLBACK_REPORT).drop_duplicates("sample_id_normalized")
        matrix = pd.read_csv(BASE_MATRIX, usecols=["sample_id_normalized", *SIGNATURE])
        frame = clinical.merge(matrix, on="sample_id_normalized", how="inner", validate="one_to_one")
        source = FALLBACK_REPORT
        LOGGER.warning("Canonical clinical fusion absent; using explicit derived source %s", source)
    frame = validate_endpoint(frame)
    values = frame[list(BASE_FEATURES)].copy()
    values[list(SIGNATURE)] = np.log1p(values[list(SIGNATURE)].clip(lower=0))
    return values, frame["binary_label"].astype(int).to_numpy()


def metrics(labels: np.ndarray, probabilities: np.ndarray) -> dict[str, float]:
    predictions = (probabilities >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(labels, predictions, labels=[0, 1]).ravel()
    return {
        "roc_auc": float(roc_auc_score(labels, probabilities)),
        "balanced_accuracy": float(balanced_accuracy_score(labels, predictions)),
        "sensitivity": float(tp / (tp + fn)) if tp + fn else np.nan,
        "specificity": float(tn / (tn + fp)) if tn + fp else np.nan,
    }


def base_pipeline(estimator: Any) -> Pipeline:
    return Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        ("model", estimator),
    ])


def engineered_frame(values: pd.DataFrame) -> pd.DataFrame:
    engineered = values.copy()
    engineered["mirna_200b_ca125"] = engineered["hsa-miR-200b-3p"] * engineered["ca125"]
    engineered["mirna_296_age"] = engineered["hsa-miR-296-3p"] * engineered["age"]
    return engineered


def model_factories() -> dict[str, Any]:
    return {
        "frozen_logistic": lambda: base_pipeline(LogisticRegression(C=0.3, class_weight="balanced", solver="liblinear", max_iter=5000, random_state=SEED)),
        "rbf_svm": lambda: GridSearchCV(
            base_pipeline(SVC(kernel="rbf", class_weight="balanced", probability=True, random_state=SEED)),
            {"model__C": [0.03, 0.1, 0.3, 1.0, 3.0], "model__gamma": ["scale", 0.01, 0.1]},
            scoring="roc_auc",
            cv=StratifiedKFold(3, shuffle=True, random_state=SEED + 1),
            n_jobs=-1,
            refit=True,
            error_score="raise",
        ),
        "shallow_random_forest": lambda: base_pipeline(RandomForestClassifier(n_estimators=300, max_depth=3, min_samples_leaf=2, class_weight="balanced", random_state=SEED, n_jobs=-1)),
        "engineered_logistic": lambda: base_pipeline(LogisticRegression(C=0.3, class_weight="balanced", solver="liblinear", max_iter=5000, random_state=SEED)),
    }


def run(values: pd.DataFrame, labels: np.ndarray) -> pd.DataFrame:
    factories = model_factories()
    outer = StratifiedKFold(5, shuffle=True, random_state=SEED)
    predictions = {name: np.full(len(labels), np.nan) for name in factories}
    for fold, (train_index, validation_index) in enumerate(outer.split(values, labels), start=1):
        LOGGER.info("Outer fold %d/5", fold)
        train = values.iloc[train_index]
        validation = values.iloc[validation_index]
        engineered_train = engineered_frame(train)
        engineered_validation = engineered_frame(validation)
        for name, factory in factories.items():
            if name == "engineered_logistic":
                estimator = factory()
                estimator.fit(engineered_train, labels[train_index])
                probability = estimator.predict_proba(engineered_validation)[:, 1]
            else:
                estimator = factory()
                estimator.fit(train, labels[train_index])
                probability = estimator.predict_proba(validation)[:, 1]
            predictions[name][validation_index] = probability
    rows = []
    for name, probability in predictions.items():
        result = metrics(labels, probability)
        result["model"] = name
        result["exceeds_90_auc"] = bool(result["roc_auc"] > 0.90)
        rows.append(result)
    return pd.DataFrame(rows)[["model", "roc_auc", "balanced_accuracy", "sensitivity", "specificity", "exceeds_90_auc"]]


def write_report(results: pd.DataFrame) -> None:
    lines = [
        "# High-Capacity Local Model Comparison",
        "",
        "All preprocessing, interaction construction, and SVM hyperparameter tuning occur inside each outer training fold. Endpoint: GSE94533 Stage I/II cancer versus benign/control.",
        "",
        "| Model | Mean OOF ROC-AUC | Balanced Accuracy | Sensitivity | Specificity | >0.90 ROC-AUC |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for row in results.itertuples(index=False):
        lines.append(f"| {row.model} | {row.roc_auc:.4f} | {row.balanced_accuracy:.4f} | {row.sensitivity:.4f} | {row.specificity:.4f} | {'Yes' if row.exceeds_90_auc else 'No'} |")
    REPORT.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    values, labels = load_training_data()
    results = run(values, labels)
    print(results.to_string(index=False, float_format=lambda value: f"{value:.4f}"))
    write_report(results)
    LOGGER.info("Wrote %s", REPORT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
