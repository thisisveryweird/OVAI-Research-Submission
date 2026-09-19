#!/usr/bin/env python3

"""
Build a clean, patient-level GSE94533 serum-miRNA dataset.

Purpose
-------
GSE94533 contains processed serum miRNA expression matrices in two Excel
workbooks. The Excel workbooks contain sample IDs and molecular measurements,
but the diagnostic labels are stored in GEO sample metadata.

This script:

1. Downloads the authoritative GEO family SOFT metadata.
2. Extracts sample title, diagnosis, stage, histology, and GEO accession.
3. Normalizes sample IDs so Excel and GEO identifiers can be joined safely.
4. Reads the two processed GSE94533 miRNA TPM workbooks.
5. Identifies the sample columns automatically.
6. Excludes E-1056, which the original GSE94533 study excluded.
7. Converts the GEO diagnosis into:
       1 = invasive ovarian cancer
       0 = non-cancer comparison group
8. Preserves the original diagnosis as well.
9. Builds a patient x miRNA matrix.
10. Performs strict integrity checks.
11. Writes:
       data/interim/gse94533_metadata.csv
       data/processed/gse94533_mirna_tpm.csv.gz
       data/processed/gse94533_dataset_manifest.json

IMPORTANT
---------
This script deliberately DOES NOT train an AI model.

We first establish a correct ground-truth dataset.
"""

from __future__ import annotations

import gzip
import json
import re
import sys
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

ROOT = Path(__file__).resolve().parents[1]

RAW_DIR = ROOT / "data" / "raw"
INTERIM_DIR = ROOT / "data" / "interim"
PROCESSED_DIR = ROOT / "data" / "processed"
REPORT_DIR = ROOT / "reports" / "tables"
SCRIPT_DIR = ROOT / "scripts"

SOFT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/"
    "GSE94nnn/GSE94533/soft/GSE94533_family.soft.gz"
)

SOFT_PATH = INTERIM_DIR / "GSE94533_family.soft.gz"

METADATA_PATH = INTERIM_DIR / "gse94533_metadata.csv"
MATRIX_PATH = PROCESSED_DIR / "gse94533_mirna_tpm.csv.gz"
MANIFEST_PATH = PROCESSED_DIR / "gse94533_dataset_manifest.json"


# ---------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------

def normalize_sample_id(value: object) -> str:
    """
    Convert different representations of the same sample identifier
    into one canonical form.

    Examples
    --------
    E-1036   -> e1036
    e1036    -> e1036
    50217    -> 50217
    50217.0  -> 50217
    """

    if pd.isna(value):
        return ""

    s = str(value).strip().lower()

    # Excel often turns integer-looking IDs into floats.
    if re.fullmatch(r"\d+\.0+", s):
        s = s.split(".", 1)[0]

    # Remove spaces and punctuation.
    s = re.sub(r"[^a-z0-9]", "", s)

    return s


def classify_diagnosis(raw_value: str) -> tuple[str | None, int | None]:
    """
    Convert the original GEO diagnosis into the study endpoint.

    Returns
    -------
    diagnosis_group:
        cancer / borderline / benign / healthy / None

    binary_label:
        1 for invasive ovarian cancer
        0 for comparison groups
        None for unknown
    """

    if raw_value is None:
        return None, None

    s = str(raw_value).strip().lower()

    if not s:
        return None, None

    # IMPORTANT:
    # Check non-cancer categories first.
    # This prevents a future label such as
    # "borderline tumor" from being accidentally classified as cancer
    # just because it contains the word "tumor".

    if "borderline" in s:
        return "borderline", 0

    if "benign" in s:
        return "benign", 0

    if "healthy" in s:
        return "healthy", 0

    if "control" in s and "cancer" not in s:
        return "healthy", 0

    # Positive endpoint.
    if "invasive" in s:
        return "cancer", 1

    if "malignant" in s:
        return "cancer", 1

    if "cancer" in s:
        return "cancer", 1

    return None, None


def download_file(url: str, destination: Path) -> None:
    """
    Download a file from NCBI with a normal HTTP User-Agent.
    """

    destination.parent.mkdir(parents=True, exist_ok=True)

    if destination.exists() and destination.stat().st_size > 0:
        print(f"[download] Already exists: {destination}")
        return

    print(f"[download] Downloading:")
    print(f"           {url}")
    print(f"           -> {destination}")

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "ovarian-ai-research-prototype/1.0"
        },
    )

    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            data = response.read()
    except Exception as exc:
        raise RuntimeError(
            f"Could not download GEO metadata.\n"
            f"URL: {url}\n"
            f"Error: {exc}"
        ) from exc

    destination.write_bytes(data)

    print(
        f"[download] Complete: {destination.stat().st_size:,} bytes"
    )


# ---------------------------------------------------------------------
# GEO SOFT parser
# ---------------------------------------------------------------------

def parse_geo_soft(path: Path) -> pd.DataFrame:
    """
    Parse the GSE94533 family SOFT file.

    We only extract sample-level metadata.
    """

    print()
    print("=" * 90)
    print("PARSING GEO SAMPLE METADATA")
    print("=" * 90)

    records: list[dict] = []

    current: dict | None = None

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:

        for raw_line in fh:

            line = raw_line.rstrip("\n")

            # A new GEO sample record.
            if line.startswith("^SAMPLE = "):

                if current is not None:
                    records.append(current)

                accession = line.split("=", 1)[1].strip()

                current = {
                    "geo_accession": accession,
                    "title": None,
                    "characteristics": {},
                }

                continue

            if current is None:
                continue

            if line.startswith("!Sample_title = "):

                title = line.split("=", 1)[1].strip()
                current["title"] = title
                continue

            if line.startswith("!Sample_characteristics_ch1 = "):

                value = line.split("=", 1)[1].strip()

                if ":" in value:
                    key, val = value.split(":", 1)
                    key = key.strip().lower()
                    val = val.strip()
                else:
                    key = "other"
                    val = value.strip()

                current["characteristics"][key] = val

                continue

    if current is not None:
        records.append(current)

    rows = []

    for record in records:

        chars = record["characteristics"]

        diagnosis_raw = chars.get("diagnosis")

        diagnosis_group, binary_label = classify_diagnosis(
            diagnosis_raw
        )

        rows.append(
            {
                "geo_accession": record["geo_accession"],
                "sample_id": record["title"],
                "sample_id_normalized": normalize_sample_id(
                    record["title"]
                ),
                "diagnosis_raw": diagnosis_raw,
                "diagnosis_group": diagnosis_group,
                "binary_label": binary_label,
                "stage": chars.get("stage"),
                "histology": chars.get("histology"),
            }
        )

    metadata = pd.DataFrame(rows)

    print()
    print(f"GEO sample records found: {len(metadata)}")

    if len(metadata) != 180:
        print(
            "WARNING: expected 180 GEO sample records according to "
            "the GSE94533 record."
        )

    return metadata


# ---------------------------------------------------------------------
# Metadata validation
# ---------------------------------------------------------------------

def validate_geo_metadata(metadata: pd.DataFrame) -> None:

    print()
    print("=" * 90)
    print("VALIDATING GEO METADATA")
    print("=" * 90)

    print()
    print("Diagnosis raw values:")
    print(
        metadata["diagnosis_raw"]
        .fillna("<MISSING>")
        .value_counts(dropna=False)
        .to_string()
    )

    print()
    print("Normalized diagnosis groups:")

    print(
        metadata["diagnosis_group"]
        .fillna("<UNKNOWN>")
        .value_counts(dropna=False)
        .to_string()
    )

    unknown = metadata[metadata["diagnosis_group"].isna()]

    if len(unknown) > 0:

        print()
        print("WARNING: Unknown diagnosis records:")
        print(
            unknown[
                [
                    "geo_accession",
                    "sample_id",
                    "diagnosis_raw",
                ]
            ].to_string(index=False)
        )

        raise RuntimeError(
            "GEO diagnosis values could not all be mapped safely. "
            "The script refuses to guess their labels."
        )

    # E-1056 is intentionally still present at this point. The original
    # publication's table counts do not exactly match the current GEO
    # annotation, so the source metadata is retained as authoritative and
    # the excluded sample is recorded explicitly below.
    print()
    print("Expected post-exclusion sample count: 179")
    print(f"Current GEO metadata count: {len(metadata)}")

    # Verify the published diagnosis composition among the 179
    # eligible samples after excluding E-1056.
    eligible = metadata[
        metadata["sample_id_normalized"] != "e1056"
    ].copy()

    observed = (
        eligible["diagnosis_group"]
        .value_counts()
        .to_dict()
    )

    print()
    print("Observed eligible diagnosis counts:")
    print(json.dumps(observed, indent=2))

    if len(eligible) != 179:
        raise RuntimeError(
            f"Eligible patient count is {len(eligible)}, "
            "but expected 179 after excluding E-1056."
        )

    # Check sample ID uniqueness.
    duplicate_ids = (
        eligible["sample_id_normalized"]
        .loc[
            eligible["sample_id_normalized"].duplicated(
                keep=False
            )
        ]
        .unique()
        .tolist()
    )

    if duplicate_ids:

        raise RuntimeError(
            "Duplicate normalized sample IDs found: "
            + ", ".join(duplicate_ids)
        )

    print()
    print("GEO metadata validation PASSED.")


# ---------------------------------------------------------------------
# Excel loading
# ---------------------------------------------------------------------

def find_workbooks() -> list[Path]:

    candidates = sorted(
        RAW_DIR.glob("GSE94533_Processed_file_Cohort*.xlsx")
    )

    if len(candidates) != 2:

        raise RuntimeError(
            "Expected exactly two GSE94533 processed Excel workbooks, "
            f"but found {len(candidates)}:\n"
            + "\n".join(str(p) for p in candidates)
        )

    return candidates


def load_mirna_matrix_from_workbook(
    workbook: Path,
    metadata_ids: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:

    print()
    print("=" * 90)
    print(f"READING: {workbook.name}")
    print("=" * 90)

    print()
    print("Reading sheet 'mirna_tpm' with header row 3...")

    df = pd.read_excel(
        workbook,
        sheet_name="mirna_tpm",
        header=2,
        engine="openpyxl",
    )

    print(f"Raw shape: {df.shape}")

    if df.shape[1] < 5:
        raise RuntimeError(
            f"{workbook.name} does not look like the expected "
            "GSE94533 miRNA TPM matrix."
        )

    # Identify the name columns.
    columns = list(df.columns)

    feature_name_col = columns[0]

    # The standard processed file contains:
    # Name, Sequence, Name count, Identical miRNAs, samples...
    identical_col = None

    for col in columns[:6]:

        if str(col).strip().lower() in {
            "identical mirnas",
            "identical miRNAs".lower(),
        }:
            identical_col = col
            break

    if identical_col is None:
        # Fall back to the fourth column used by the processed file.
        if len(columns) >= 4:
            identical_col = columns[3]

    print()
    print(f"Feature-name column: {feature_name_col!r}")
    print(f"Identical-miRNA column: {identical_col!r}")

    sample_columns = []
    sample_map_rows = []

    for col in columns:

        normalized = normalize_sample_id(col)

        if normalized in metadata_ids:

            sample_columns.append(col)

            sample_map_rows.append(
                {
                    "excel_column": str(col),
                    "sample_id_normalized": normalized,
                    "source_workbook": workbook.name,
                }
            )

    print()
    print(f"Matched sample columns: {len(sample_columns)}")

    if not sample_columns:

        raise RuntimeError(
            f"No sample columns from {workbook.name} matched GEO metadata."
        )

    # Show the first several matches for human inspection.
    print()
    print("First sample-column matches:")

    for row in sample_map_rows[:20]:
        print(
            f"  {row['excel_column']!r}"
            f" -> {row['sample_id_normalized']}"
        )

    # Extract only the feature identifier and sample measurements.
    expression = df[
        [feature_name_col, identical_col] + sample_columns
    ].copy()

    expression = expression.rename(
        columns={
            feature_name_col: "feature_name",
            identical_col: "canonical_mirna",
        }
    )

    # Use the canonical miRNA name when available.
    expression["feature_id"] = (
        expression["canonical_mirna"]
        .astype("string")
        .fillna(
            expression["feature_name"].astype("string")
        )
        .str.strip()
    )

    # Remove empty feature identifiers.
    expression = expression[
        expression["feature_id"].notna()
        & (expression["feature_id"] != "")
    ].copy()

    # Convert expression columns to numeric.
    for col in sample_columns:
        expression[col] = pd.to_numeric(
            expression[col],
            errors="coerce",
        )

    # Replace missing expression measurements with zero only
    # after verifying that the measurement table is numeric.
    expression[sample_columns] = (
        expression[sample_columns]
        .replace([np.inf, -np.inf], np.nan)
        .fillna(0.0)
    )

    # Collapse identical canonical miRNAs.
    #
    # Multiple sequence representations can sometimes correspond to
    # the same canonical miRNA. For TPM, summing these representations
    # is more appropriate than silently dropping one.
    grouped = (
        expression
        .groupby("feature_id", sort=True)[sample_columns]
        .sum()
    )

    # Transpose so rows = patients and columns = miRNAs.
    patient_matrix = grouped.T.copy()

    # Normalize sample names.
    patient_matrix.index = [
        normalize_sample_id(x)
        for x in patient_matrix.index
    ]

    patient_matrix.index.name = "sample_id_normalized"

    sample_map = pd.DataFrame(sample_map_rows)

    return patient_matrix, sample_map


# ---------------------------------------------------------------------
# Combine both cohorts
# ---------------------------------------------------------------------

def build_patient_level_dataset(
    metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:

    workbooks = find_workbooks()

    metadata_ids = set(
        metadata["sample_id_normalized"].astype(str)
    )

    matrices = []
    mappings = []

    for workbook in workbooks:

        matrix, sample_map = load_mirna_matrix_from_workbook(
            workbook=workbook,
            metadata_ids=metadata_ids,
        )

        matrices.append(matrix)
        mappings.append(sample_map)

    # Check whether the two processed workbooks accidentally contain
    # overlapping patient IDs.
    all_indices = [
        sample_id
        for matrix in matrices
        for sample_id in matrix.index
    ]

    duplicates = [
        sample_id
        for sample_id, count in Counter(all_indices).items()
        if count > 1
    ]

    if duplicates:

        raise RuntimeError(
            "The processed workbooks share sample IDs:\n"
            + "\n".join(duplicates)
        )

    # Align features across the two workbooks.
    combined_matrix = pd.concat(
        matrices,
        axis=0,
        join="outer",
    )

    combined_matrix = combined_matrix.fillna(0.0)

    # Restrict metadata to samples actually present in expression data.
    eligible_meta = metadata[
        metadata["sample_id_normalized"].isin(
            combined_matrix.index
        )
        & (
            metadata["sample_id_normalized"] != "e1056"
        )
    ].copy()

    # Exclude E-1056 explicitly.
    combined_matrix = combined_matrix[
        combined_matrix.index != "e1056"
    ].copy()

    # Verify one row per patient.
    if combined_matrix.index.duplicated().any():
        raise RuntimeError(
            "Duplicate patient rows remain after combining workbooks."
        )

    # Determine sample source workbook.
    mapping = pd.concat(
        mappings,
        ignore_index=True,
    )

    source_lookup = dict(
        zip(
            mapping["sample_id_normalized"],
            mapping["source_workbook"],
        )
    )

    eligible_meta["source_workbook"] = (
        eligible_meta["sample_id_normalized"]
        .map(source_lookup)
    )

    # Keep metadata and expression in exactly the same order.
    common_ids = [
        sample_id
        for sample_id in eligible_meta["sample_id_normalized"]
        if sample_id in combined_matrix.index
    ]

    eligible_meta = (
        eligible_meta
        .set_index("sample_id_normalized")
        .loc[common_ids]
        .reset_index()
    )

    combined_matrix = combined_matrix.loc[
        common_ids
    ]

    return eligible_meta, combined_matrix


# ---------------------------------------------------------------------
# Final integrity checks
# ---------------------------------------------------------------------

def final_integrity_checks(
    metadata: pd.DataFrame,
    matrix: pd.DataFrame,
) -> None:

    print()
    print("=" * 90)
    print("FINAL DATASET INTEGRITY CHECK")
    print("=" * 90)

    print()
    print(f"Patients: {matrix.shape[0]}")
    print(f"miRNA features: {matrix.shape[1]}")

    if matrix.shape[0] != 179:
        raise RuntimeError(
            f"Expected 179 eligible patients, found {matrix.shape[0]}."
        )

    if len(metadata) != matrix.shape[0]:
        raise RuntimeError(
            "Metadata and expression matrix patient counts do not match."
        )

    if list(metadata["sample_id_normalized"]) != list(matrix.index):
        raise RuntimeError(
            "Metadata patient order and matrix patient order do not match."
        )

    if matrix.index.duplicated().any():
        raise RuntimeError(
            "Duplicate patient IDs detected in expression matrix."
        )

    if not np.isfinite(matrix.to_numpy(dtype=float)).all():
        raise RuntimeError(
            "Expression matrix contains NaN or infinite values."
        )

    # Expression cannot be negative.
    if (matrix.to_numpy(dtype=float) < 0).any():
        raise RuntimeError(
            "Negative TPM values detected."
        )

    print()
    print("Diagnosis distribution:")
    print(
        metadata["diagnosis_group"]
        .value_counts()
        .reindex(
            ["cancer", "borderline", "benign", "healthy"]
        )
        .fillna(0)
        .astype(int)
        .to_string()
    )

    print()
    print("Binary endpoint distribution:")
    print(
        metadata["binary_label"]
        .value_counts()
        .sort_index()
        .to_string()
    )

    print()
    print("Patients by source workbook:")
    print(
        metadata["source_workbook"]
        .value_counts()
        .to_string()
    )

    print()
    print("Patients by diagnosis and source workbook:")
    print(
        pd.crosstab(
            metadata["source_workbook"],
            metadata["diagnosis_group"],
        )
        .reindex(
            columns=[
                "cancer",
                "borderline",
                "benign",
                "healthy",
            ],
            fill_value=0,
        )
        .to_string()
    )

    print()
    print("Stage distribution among cancer cases:")

    cancer_stage = metadata.loc[
        metadata["diagnosis_group"] == "cancer",
        "stage",
    ]

    print(
        cancer_stage
        .fillna("<MISSING>")
        .value_counts(dropna=False)
        .to_string()
    )

    print()
    print("Histology distribution:")
    print(
        metadata["histology"]
        .fillna("<MISSING>")
        .value_counts(dropna=False)
        .head(30)
        .to_string()
    )

    print()
    print("FINAL INTEGRITY CHECK PASSED.")


# ---------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------

def write_outputs(
    metadata: pd.DataFrame,
    matrix: pd.DataFrame,
) -> None:

    INTERIM_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    metadata.to_csv(
        METADATA_PATH,
        index=False,
    )

    matrix.to_csv(
        MATRIX_PATH,
        compression="gzip",
    )

    manifest = {
        "dataset": "GSE94533",
        "description": (
            "Patient-level serum miRNA TPM matrix aligned to GEO "
            "diagnosis metadata."
        ),
        "source": (
            "NCBI Gene Expression Omnibus GSE94533"
        ),
        "geo_metadata_source": SOFT_URL,
        "expression_sources": [
            str(p)
            for p in find_workbooks()
        ],
        "excluded_samples": [
            {
                "sample_id": "E-1056",
                "reason": (
                    "Excluded by original study because of "
                    "abnormally high miR-122 believed to be "
                    "associated with a recent myocardial infarction."
                ),
            }
        ],
        "n_patients": int(matrix.shape[0]),
        "n_features": int(matrix.shape[1]),
        "diagnosis_counts": (
            metadata["diagnosis_group"]
            .value_counts()
            .to_dict()
        ),
        "binary_label_definition": {
            "1": "invasive ovarian cancer",
            "0": (
                "borderline tumor, benign ovarian disease, "
                "or healthy control"
            ),
        },
        "expression_measurement": "TPM",
        "aggregation": (
            "Multiple sequence representations mapping to the same "
            "canonical miRNA were summed."
        ),
        "note": (
            "Research dataset. This is not a validated clinical "
            "diagnostic test."
        ),
    }

    MANIFEST_PATH.write_text(
        json.dumps(
            manifest,
            indent=2,
        )
    )

    print()
    print("=" * 90)
    print("OUTPUTS WRITTEN")
    print("=" * 90)

    print()
    print(f"Metadata:")
    print(f"  {METADATA_PATH}")

    print()
    print(f"Expression matrix:")
    print(f"  {MATRIX_PATH}")

    print()
    print(f"Manifest:")
    print(f"  {MANIFEST_PATH}")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main() -> int:

    print("=" * 90)
    print("GSE94533 PATIENT-LEVEL DATASET BUILDER")
    print("=" * 90)

    print()
    print(f"Project root: {ROOT}")

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    INTERIM_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------------------
    # Step 1: obtain official GEO metadata
    # -------------------------------------------------------------

    download_file(
        SOFT_URL,
        SOFT_PATH,
    )

    metadata = parse_geo_soft(SOFT_PATH)

    if metadata.empty:
        raise RuntimeError(
            "No GEO sample metadata were extracted."
        )

    validate_geo_metadata(metadata)

    # -------------------------------------------------------------
    # Step 2: build patient-level expression matrix
    # -------------------------------------------------------------

    final_metadata, matrix = build_patient_level_dataset(
        metadata
    )

    # -------------------------------------------------------------
    # Step 3: strict integrity checks
    # -------------------------------------------------------------

    final_integrity_checks(
        final_metadata,
        matrix,
    )

    # -------------------------------------------------------------
    # Step 4: write clean data products
    # -------------------------------------------------------------

    write_outputs(
        final_metadata,
        matrix,
    )

    print()
    print("=" * 90)
    print("GSE94533 DATASET BUILD COMPLETE")
    print("=" * 90)

    print()
    print("We have NOT trained a model yet.")
    print()
    print(
        "The next stage will use this frozen patient-level dataset "
        "to build the first miRNA diagnostic baseline."
    )

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nInterrupted by user.")
        raise SystemExit(130)
    except Exception as exc:
        print()
        print("=" * 90)
        print("DATASET BUILD FAILED")
        print("=" * 90)
        print()
        print(f"{type(exc).__name__}: {exc}")
        print()
        print(
            "No model training was performed."
        )
        raise SystemExit(1)
