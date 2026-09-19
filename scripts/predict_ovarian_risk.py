#!/usr/bin/env python3
"""Run production ovarian-cancer risk inference from seven raw inputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_ARTIFACT_DIR = PROJECT / "production_artifacts"
FEATURES = (
    "hsa-miR-200b-3p",
    "hsa-miR-296-3p",
    "hsa-miR-6826-5p",
    "hsa-miR-32-3p",
    "hsa-miR-1255a",
)
CLINICAL_FEATURES = ("age", "ca125")


def parse_number(value: str) -> float:
    if value.strip().lower() in {"", "none", "nan", "na", "null"}:
        return np.nan
    number = float(value)
    if not np.isfinite(number):
        raise ValueError(f"Non-finite input: {value}")
    return number


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for feature in FEATURES:
        parser.add_argument(f"--{feature}", required=True, type=parse_number)
    parser.add_argument("--age", required=True, type=parse_number)
    parser.add_argument("--ca125", required=True, type=parse_number)
    parser.add_argument("--artifact-dir", type=Path, default=DEFAULT_ARTIFACT_DIR)
    parser.add_argument("--high-confidence", type=float, default=0.80, help="Confidence threshold for the high-confidence flag.")
    args = parser.parse_args()
    if not 0.5 <= args.high_confidence < 1.0:
        raise ValueError("--high-confidence must be at least 0.5 and below 1.0")
    artifact_dir = args.artifact_dir
    model = joblib.load(artifact_dir / "frozen_model.joblib")
    scaler = joblib.load(artifact_dir / "feature_scaler.joblib")
    defaults = json.loads((artifact_dir / "imputation_defaults.json").read_text())
    row = {feature: getattr(args, feature.replace("-", "_")) for feature in FEATURES}
    row.update({"age": args.age, "ca125": args.ca125})
    values = pd.DataFrame([row], columns=list(FEATURES) + list(CLINICAL_FEATURES))
    values.loc[:, FEATURES] = np.log1p(values.loc[:, FEATURES].clip(lower=0))
    for feature in CLINICAL_FEATURES:
        values[feature] = values[feature].fillna(float(defaults[feature]))
    probability = float(model.predict_proba(scaler.transform(values.to_numpy(dtype=float)))[:, 1][0])
    classification = "Cancer" if probability >= 0.5 else "Healthy/Benign"
    confidence = max(probability, 1.0 - probability)
    result = {
        "cancer_probability": probability,
        "classification": classification,
        "confidence": confidence,
        "high_confidence": bool(confidence >= args.high_confidence),
        "threshold": 0.5,
        "high_confidence_threshold": args.high_confidence,
    }
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
