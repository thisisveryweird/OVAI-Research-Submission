#!/usr/bin/env python3
"""Run research-only, schema-gated model inference from seven raw inputs."""

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


EXPECTED_FEATURES = FEATURES + CLINICAL_FEATURES


def validate_artifact_schema(manifest, scaler=None, model=None) -> None:
    expected = list(EXPECTED_FEATURES)
    if not isinstance(manifest, dict) or manifest.get("features") != expected:
        raise ValueError("Artifact manifest must declare exact ordered native features")
    for artifact in (scaler, model):
        if artifact is not None:
            dimension = getattr(artifact, "n_features_in_", None)
            if dimension != len(expected):
                raise ValueError("Artifact feature dimension missing or incompatible")


    if model is not None:
        classes = getattr(model, "classes_", None)
        if classes is None or np.asarray(classes).shape != (2,) or not np.array_equal(classes, [0, 1]):
            raise ValueError("Model classes_ must be exactly [0, 1] in that order")


def validate_probabilities(probabilities):
    values = np.asarray(probabilities, dtype=float)
    if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
        raise ValueError("Model probabilities must be finite and within [0, 1]")
    return values


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
    parser.add_argument("--high-confidence", type=float, default=0.80, help="Model-score flag cutoff; not calibrated confidence.")
    args = parser.parse_args()
    if not 0.5 <= args.high_confidence < 1.0:
        raise ValueError("--high-confidence must be at least 0.5 and below 1.0")
    artifact_dir = args.artifact_dir
    manifest = json.loads((artifact_dir / "manifest.json").read_text())
    validate_artifact_schema(manifest)
    model = joblib.load(artifact_dir / "frozen_model.joblib")
    scaler = joblib.load(artifact_dir / "feature_scaler.joblib")
    validate_artifact_schema(manifest, scaler, model)
    defaults = json.loads((artifact_dir / "imputation_defaults.json").read_text())
    row = {feature: getattr(args, feature.replace("-", "_")) for feature in FEATURES}
    row.update({"age": args.age, "ca125": args.ca125})
    for feature, value in row.items():
        if feature in FEATURES and (not np.isfinite(value) or value < 0):
            raise ValueError(f"{feature} must be finite and non-negative")
        if feature in CLINICAL_FEATURES:
            default = float(defaults[feature])
            if not np.isfinite(default) or default < 0:
                raise ValueError(f"Invalid default for {feature}")
            if np.isinf(value) or (not np.isnan(value) and value < 0):
                raise ValueError(f"Invalid measured {feature}")
    values = pd.DataFrame([row], columns=list(FEATURES) + list(CLINICAL_FEATURES))
    values.loc[:, FEATURES] = np.log1p(values.loc[:, FEATURES].clip(lower=0))
    for feature in CLINICAL_FEATURES:
        values[feature] = values[feature].fillna(float(defaults[feature]))
    probability = float(model.predict_proba(scaler.transform(values.to_numpy(dtype=float)))[:, 1][0])
    validate_probabilities(probability)
    classification = "Above legacy research cutoff" if probability >= 0.5 else "Below legacy research cutoff"
    confidence = max(probability, 1.0 - probability)
    result = {
        "research_model_probability": probability,
        "research_classification": classification,
        "research_only": True,
        "score_note": "max(p, 1-p) is not calibrated confidence",
        "model_class_score": confidence,
        "high_model_score": bool(confidence >= args.high_confidence),
        "threshold": 0.5,
        "threshold_selection_provenance": "Unknown",
        "high_confidence_threshold": args.high_confidence,
    }
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
