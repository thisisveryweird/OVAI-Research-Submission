#!/usr/bin/env python3
"""Freeze the GSE94533 Top-5 miRNA plus Age/CA-125 inference pipeline."""

from __future__ import annotations

import argparse
import json
import logging
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=FutureWarning)

PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = PROJECT / "data" / "processed" / "gse94533_14mirna_clinical.csv.gz"
MATRIX_PATH = PROJECT / "data" / "processed" / "gse94533_mirna_tpm.csv.gz"
DEFAULT_SOURCE_REPORT = PROJECT / "reports" / "error_analysis_predictions.csv"
ARTIFACT_DIR = PROJECT / "production_artifacts"
FEATURES = (
    "hsa-miR-200b-3p",
    "hsa-miR-296-3p",
    "hsa-miR-6826-5p",
    "hsa-miR-32-3p",
    "hsa-miR-1255a",
)
CLINICAL_FEATURES = ("age", "ca125")
INTERACTION_FEATURES = ("mirna_200b_ca125", "mirna_296_age")
MODEL_FEATURES = FEATURES + CLINICAL_FEATURES + INTERACTION_FEATURES
SEED = 20260919
LOGGER = logging.getLogger("ovai.production_export")


def validate_data(frame: pd.DataFrame) -> pd.DataFrame:
    required = set(FEATURES) | set(CLINICAL_FEATURES) | {"binary_label", "diagnosis", "stage"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"Production source missing columns: {missing}")
    keep = frame["diagnosis"].eq("Cancer") & frame["stage"].eq("Stage I/II")
    keep |= frame["diagnosis"].isin(["Benign", "Control"])
    frame = frame.loc[keep].copy()
    labels = frame["binary_label"].astype(int)
    if len(frame) != 113 or labels.sum() != 52 or (labels == 0).sum() != 61:
        raise ValueError(
            "Unexpected Phase 8 endpoint composition: "
            f"n={len(frame)}, positives={int(labels.sum())}, negatives={int((labels == 0).sum())}"
        )
    if frame["sample_id_normalized"].duplicated().any():
        raise ValueError("Duplicate patient IDs in production source")
    return frame.reset_index(drop=True)


def load_data(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Curated clinical fusion not found: {path}")
    return validate_data(pd.read_csv(path))


def load_explicit_fallback(path: Path) -> pd.DataFrame:
    """Load a derived report only when explicitly requested by the operator."""

    frame = pd.read_csv(path)
    frame = frame.drop_duplicates("sample_id_normalized", keep="first")
    if not set(FEATURES).issubset(frame.columns):
        if not MATRIX_PATH.exists():
            raise FileNotFoundError(f"Fallback requires the processed miRNA matrix: {MATRIX_PATH}")
        matrix = pd.read_csv(MATRIX_PATH, usecols=["sample_id_normalized", *FEATURES])
        frame = frame.drop(columns=[feature for feature in FEATURES if feature in frame], errors="ignore").merge(
            matrix, on="sample_id_normalized", how="inner", validate="one_to_one"
        )
    return validate_data(frame)


def transformed_values(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mirna = np.log1p(frame[list(FEATURES)].clip(lower=0).to_numpy(dtype=float))
    clinical = frame[list(CLINICAL_FEATURES)].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    imputer = SimpleImputer(strategy="median")
    imputed = imputer.fit_transform(clinical)
    engineered = np.column_stack([
        mirna[:, 0] * imputed[:, 1],
        mirna[:, 1] * imputed[:, 0],
    ])
    return np.column_stack([mirna, imputed, engineered]), imputer.statistics_, clinical


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA, help="Canonical curated clinical fusion input.")
    parser.add_argument("--source-report", type=Path, help="Explicit fallback derived report, e.g. error_analysis_predictions.csv.")
    parser.add_argument("--output-dir", type=Path, default=ARTIFACT_DIR)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    if args.source_report:
        frame = load_explicit_fallback(args.source_report)
        source = str(args.source_report)
        LOGGER.warning("Using explicitly requested derived source report: %s", source)
    else:
        frame = load_data(args.data)
        source = str(args.data)
    values, medians, _ = transformed_values(frame)
    labels = frame["binary_label"].astype(int).to_numpy()
    scaler = StandardScaler().fit(values)
    scaled = scaler.transform(values)
    model = LogisticRegression(C=0.3, class_weight="balanced", solver="liblinear", max_iter=5000, random_state=SEED)
    model.fit(scaled, labels)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(model, args.output_dir / "frozen_model.joblib")
    joblib.dump(scaler, args.output_dir / "feature_scaler.joblib")
    (args.output_dir / "imputation_defaults.json").write_text(json.dumps({
        "age": float(medians[0]),
        "ca125": float(medians[1]),
        "clinical_features": list(CLINICAL_FEATURES),
        "policy": "median values fitted on the curated GSE94533 Stage I/II cancer versus benign/control training cohort",
    }, indent=2) + "\n")
    manifest = {
        "source": source,
        "endpoint": "GSE94533 Stage I/II invasive cancer versus benign/control",
        "n_training_samples": int(len(frame)),
        "positive_samples": int(labels.sum()),
        "negative_samples": int((labels == 0).sum()),
        "features": list(MODEL_FEATURES),
        "miRNA_transform": "log1p(max(value, 0))",
        "interaction_features": {
            "mirna_200b_ca125": "hsa-miR-200b-3p_log1p * ca125_imputed",
            "mirna_296_age": "hsa-miR-296-3p_log1p * age_imputed",
        },
        "clinical_transform": "median imputation using imputation_defaults.json",
        "scaler": "StandardScaler fitted after imputation",
        "model": "LogisticRegression(C=0.3, class_weight=balanced, solver=liblinear)",
        "random_state": SEED,
    }
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    LOGGER.info("Exported production artifacts to %s", args.output_dir)
    LOGGER.info("Training samples=%d; positives=%d; negatives=%d", len(frame), int(labels.sum()), int((labels == 0).sum()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
