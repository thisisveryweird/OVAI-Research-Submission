#!/usr/bin/env python3

"""
Create the exact 14-miRNA benchmark matrix used by the historical
eLife ovarian-cancer neural-network model.

This script DOES NOT train a model.

It:

1. Loads the frozen 179-patient GSE94533 dataset.
2. Uses the exact 14 miRNAs recovered from the original PMML.
3. Verifies every feature exists exactly once.
4. Joins diagnosis labels.
5. Writes a compact patient-level benchmark matrix.
6. Produces detailed provenance/audit information.

The goal is to create a reproducible starting point for:
    historical reproduction
        ->
    modern baseline
        ->
    multimodal model
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


# ================================================================
# PATHS
# ================================================================

ROOT = Path.home() / "ovarian-ai"

MATRIX_PATH = (
    ROOT
    / "data"
    / "processed"
    / "gse94533_mirna_tpm.csv.gz"
)

METADATA_PATH = (
    ROOT
    / "data"
    / "processed"
    / "gse94533_metadata.csv"
)

FEATURE_JSON = (
    ROOT
    / "reports"
    / "tables"
    / "original_14_mirna_features.json"
)

OUTPUT_DIR = (
    ROOT
    / "data"
    / "processed"
)

REPORT_DIR = (
    ROOT
    / "reports"
    / "tables"
)

OUTPUT_MATRIX = (
    OUTPUT_DIR
    / "gse94533_14mirna_benchmark.csv.gz"
)

OUTPUT_AUDIT = (
    REPORT_DIR
    / "gse94533_14mirna_benchmark_audit.json"
)


# ================================================================
# FALLBACK EXACT FEATURE LIST
# ================================================================
#
# This is intentionally written explicitly as a second safety
# mechanism. The PMML-derived JSON remains the provenance source,
# but these are the exact 14 inputs we already verified.
#
# The script will compare the JSON against this list.
#

EXPECTED_14 = [
    "hsa-miR-23b-3p",
    "hsa-miR-29a-3p",
    "hsa-miR-32-5p",
    "hsa-miR-92a-3p",
    "hsa-miR-150-5p",
    "hsa-miR-200a-3p",
    "hsa-miR-200c-3p",
    "hsa-miR-203a",
    "hsa-miR-320c",
    "hsa-miR-320d",
    "hsa-miR-335-5p",
    "hsa-miR-450b-5p",
    "hsa-miR-1246",
    "hsa-miR-1307-5p",
]


# ================================================================
# LOAD FEATURE PROVENANCE
# ================================================================

def load_pmml_features() -> list[str]:

    if not FEATURE_JSON.exists():

        raise RuntimeError(
            "Could not find:\n"
            f"{FEATURE_JSON}\n\n"
            "Run extract_original_14mirna_model.py first."
        )

    data = json.loads(
        FEATURE_JSON.read_text()
    )

    fields = data.get(
        "mirna_named_fields",
        [],
    )

    fields = list(
        dict.fromkeys(
            str(x).strip()
            for x in fields
        )
    )

    return fields


# ================================================================
# LOAD DATA
# ================================================================

def load_matrix() -> pd.DataFrame:

    if not MATRIX_PATH.exists():

        raise RuntimeError(
            f"Missing expression matrix:\n{MATRIX_PATH}"
        )

    matrix = pd.read_csv(
        MATRIX_PATH,
        index_col=0,
    )

    return matrix


def load_metadata() -> pd.DataFrame:

    if not METADATA_PATH.exists():

        raise RuntimeError(
            f"Missing metadata:\n{METADATA_PATH}"
        )

    metadata = pd.read_csv(
        METADATA_PATH,
    )

    required = {
        "sample_id_normalized",
        "diagnosis",
        "binary_label",
    }

    missing = required - set(
        metadata.columns
    )

    if missing:

        raise RuntimeError(
            f"Metadata missing columns: {missing}"
        )

    return metadata


# ================================================================
# MAIN BUILD
# ================================================================

def main() -> None:

    print("=" * 90)
    print("BUILDING EXACT 14-miRNA BENCHMARK DATASET")
    print("=" * 90)

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    REPORT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ------------------------------------------------------------
    # Load PMML-derived features.
    # ------------------------------------------------------------

    pmml_features = load_pmml_features()

    print()
    print("Features recovered from original PMML:")

    for feature in pmml_features:
        print(
            f"  {feature}"
        )

    # ------------------------------------------------------------
    # Require exactly the known 14 inputs.
    # ------------------------------------------------------------

    if set(pmml_features) != set(EXPECTED_14):

        print()
        print("EXPECTED:")
        for x in EXPECTED_14:
            print(" ", x)

        print()
        print("PMML:")
        for x in pmml_features:
            print(" ", x)

        raise RuntimeError(
            "PMML-derived feature set does not match the expected "
            "14-feature historical model."
        )

    if len(pmml_features) != 14:

        raise RuntimeError(
            f"Expected 14 features, found {len(pmml_features)}."
        )

    # ------------------------------------------------------------
    # Load patient expression data.
    # ------------------------------------------------------------

    matrix = load_matrix()

    print()
    print(
        f"Loaded expression matrix: "
        f"{matrix.shape[0]} patients x "
        f"{matrix.shape[1]} features"
    )

    # ------------------------------------------------------------
    # Confirm all 14 features exist.
    # ------------------------------------------------------------

    missing = [
        feature
        for feature in pmml_features
        if feature not in matrix.columns
    ]

    if missing:

        raise RuntimeError(
            "Historical model features missing from cleaned dataset:\n"
            + "\n".join(
                f"  {x}"
                for x in missing
            )
        )

    duplicate_features = [
        feature
        for feature in pmml_features
        if list(matrix.columns).count(feature) != 1
    ]

    if duplicate_features:

        raise RuntimeError(
            "Historical features are duplicated in expression matrix:\n"
            + "\n".join(
                f"  {x}"
                for x in duplicate_features
            )
        )

    # ------------------------------------------------------------
    # Extract exact 14 columns IN PMML ORDER.
    # ------------------------------------------------------------

    benchmark = matrix[
        pmml_features
    ].copy()

    # ------------------------------------------------------------
    # Load metadata.
    # ------------------------------------------------------------

    metadata = load_metadata()

    if len(metadata) != len(matrix):

        raise RuntimeError(
            "Metadata patient count does not match expression matrix."
        )

    # ------------------------------------------------------------
    # Verify patient alignment.
    # ------------------------------------------------------------

    metadata_ids = (
        metadata[
            "sample_id_normalized"
        ]
        .astype(str)
        .tolist()
    )

    matrix_ids = [
        str(x)
        for x in matrix.index
    ]

    if metadata_ids != matrix_ids:

        raise RuntimeError(
            "Patient ordering differs between expression matrix "
            "and metadata."
        )

    # ------------------------------------------------------------
    # Add metadata to benchmark.
    # ------------------------------------------------------------

    benchmark.insert(
        0,
        "binary_label",
        metadata["binary_label"].astype(int).values,
    )

    benchmark.insert(
        0,
        "diagnosis",
        metadata["diagnosis"].astype(str).values,
    )

    benchmark.insert(
        0,
        "sample_id",
        metadata[
            "sample_id_normalized"
        ].astype(str).values,
    )

    # Add stage/histology where available.
    if "stage" in metadata.columns:

        benchmark["stage"] = (
            metadata["stage"]
            .astype("string")
            .values
        )

    if "histology" in metadata.columns:

        benchmark["histology"] = (
            metadata["histology"]
            .astype("string")
            .values
        )

    # ------------------------------------------------------------
    # Numeric checks.
    # ------------------------------------------------------------

    feature_values = benchmark[
        pmml_features
    ].to_numpy(
        dtype=float
    )

    if not np.isfinite(
        feature_values
    ).all():

        raise RuntimeError(
            "14-miRNA benchmark contains NaN or infinite values."
        )

    if (
        feature_values < 0
    ).any():

        raise RuntimeError(
            "14-miRNA TPM matrix contains negative values."
        )

    # ------------------------------------------------------------
    # Confirm patient uniqueness.
    # ------------------------------------------------------------

    if benchmark[
        "sample_id"
    ].duplicated().any():

        raise RuntimeError(
            "Duplicate patient IDs found."
        )

    if len(benchmark) != 179:

        raise RuntimeError(
            f"Expected 179 patients, found {len(benchmark)}."
        )

    # ------------------------------------------------------------
    # Summary.
    # ------------------------------------------------------------

    diagnosis_counts = (
        benchmark[
            "diagnosis"
        ]
        .value_counts()
        .to_dict()
    )

    binary_counts = (
        benchmark[
            "binary_label"
        ]
        .value_counts()
        .sort_index()
        .to_dict()
    )

    audit = {

        "dataset": "GSE94533",

        "purpose": (
            "Exact 14-miRNA feature matrix corresponding to "
            "the historical eLife ovarian-cancer neural-network model."
        ),

        "n_patients": int(
            len(benchmark)
        ),

        "n_features": 14,

        "features": pmml_features,

        "original_model_architecture": {
            "input_features": 14,
            "hidden_neurons": 7,
            "output_classes": 2,
            "model_name": (
                "Dataset complete with classificat_MLP 14-7-2"
            ),
        },

        "feature_provenance": str(
            FEATURE_JSON
        ),

        "source_expression_matrix": str(
            MATRIX_PATH
        ),

        "source_metadata": str(
            METADATA_PATH
        ),

        "diagnosis_counts": diagnosis_counts,

        "binary_counts": binary_counts,

        "binary_endpoint": {
            "1": "Cancer",
            "0": "Benign, Borderline, or Control",
        },

        "validation_checks": {
            "patient_count_is_179": (
                len(benchmark) == 179
            ),
            "exactly_14_features": (
                len(pmml_features) == 14
            ),
            "all_features_present": (
                len(missing) == 0
            ),
            "all_values_finite": bool(
                np.isfinite(
                    feature_values
                ).all()
            ),
            "all_values_nonnegative": bool(
                (
                    feature_values >= 0
                ).all()
            ),
            "unique_patient_ids": (
                not benchmark[
                    "sample_id"
                ].duplicated().any()
            ),
            "patient_order_aligned": True,
        },

        "important_note": (
            "The labels in this benchmark come from the current "
            "GSE94533 GEO metadata used to construct the cleaned dataset. "
            "They are intentionally not altered to force agreement with "
            "historical class counts."
        ),
    }

    # ------------------------------------------------------------
    # Save compact matrix.
    # ------------------------------------------------------------

    benchmark.to_csv(
        OUTPUT_MATRIX,
        index=False,
        compression="gzip",
    )

    OUTPUT_AUDIT.write_text(
        json.dumps(
            audit,
            indent=2,
        )
    )

    # ------------------------------------------------------------
    # Terminal report.
    # ------------------------------------------------------------

    print()
    print("=" * 90)
    print("14-miRNA BENCHMARK DATASET")
    print("=" * 90)

    print()
    print(
        f"Patients: {len(benchmark)}"
    )

    print(
        f"Features: {len(pmml_features)}"
    )

    print()
    print("Feature order:")

    for i, feature in enumerate(
        pmml_features,
        1,
    ):

        print(
            f"  {i:2d}. {feature}"
        )

    print()
    print("Diagnosis counts:")

    for diagnosis, count in diagnosis_counts.items():

        print(
            f"  {diagnosis:12s} {count}"
        )

    print()
    print("Binary counts:")

    print(
        f"  0 = non-cancer   "
        f"{binary_counts.get(0, 0)}"
    )

    print(
        f"  1 = cancer       "
        f"{binary_counts.get(1, 0)}"
    )

    print()
    print("Validation:")

    for name, result in audit[
        "validation_checks"
    ].items():

        print(
            f"  {name:30s} {result}"
        )

    print()
    print("=" * 90)
    print("FILES CREATED")
    print("=" * 90)

    print()
    print(
        f"Benchmark matrix:\n"
        f"  {OUTPUT_MATRIX}"
    )

    print()
    print(
        f"Audit:\n"
        f"  {OUTPUT_AUDIT}"
    )

    print()
    print("=" * 90)
    print("NO MODEL TRAINING PERFORMED")
    print("=" * 90)

    print(
        """
We now have a frozen 14-feature patient-level dataset.

The next scientific step is to reproduce/benchmark the historical
classifier and establish a modern baseline using proper patient-level
cross-validation.

We should NOT immediately train a large neural network on all 2,582
features because 179 patients versus thousands of features creates a
serious overfitting risk.
"""
    )


if __name__ == "__main__":
    main()
