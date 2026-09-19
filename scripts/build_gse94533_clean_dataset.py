#!/usr/bin/env python3

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


# ================================================================
# PROJECT PATHS
# ================================================================

ROOT = Path(__file__).resolve().parents[1]

RAW_DIR = ROOT / "data" / "raw"
INTERIM_DIR = ROOT / "data" / "interim"
PROCESSED_DIR = ROOT / "data" / "processed"
REPORT_DIR = ROOT / "reports" / "tables"

METADATA_INPUT = INTERIM_DIR / "gse94533_geo_metadata_full.csv"

METADATA_OUTPUT = PROCESSED_DIR / "gse94533_metadata.csv"
MATRIX_OUTPUT = PROCESSED_DIR / "gse94533_mirna_tpm.csv.gz"
AUDIT_OUTPUT = REPORT_DIR / "gse94533_clean_dataset_audit.json"


# ================================================================
# SAMPLE-ID NORMALIZATION
# ================================================================

def normalize_sample_id(value: object) -> str:

    if pd.isna(value):
        return ""

    s = str(value).strip().lower()

    # Correct Excel-style numeric IDs.
    #
    # Examples:
    #   50057.0 -> 50057
    #   50217.0 -> 50217
    #
    if re.fullmatch(r"\d+\.0+", s):
        s = s.split(".", 1)[0]

    # Remove punctuation.
    #
    # E-1001 -> e1001
    # 50057  -> 50057
    #
    s = re.sub(r"[^a-z0-9]", "", s)

    return s


# ================================================================
# FIND INPUT WORKBOOKS
# ================================================================

def find_workbooks() -> list[Path]:

    files = sorted(
        RAW_DIR.glob("GSE94533_Processed_file_Cohort*.xlsx")
    )

    if len(files) != 2:

        raise RuntimeError(
            "Expected exactly two processed GSE94533 workbooks.\n"
            f"Found {len(files)}:\n"
            + "\n".join(str(x) for x in files)
        )

    return files


# ================================================================
# LOAD GEO METADATA
# ================================================================

def load_metadata() -> pd.DataFrame:

    if not METADATA_INPUT.exists():

        raise RuntimeError(
            f"Could not find GEO metadata at:\n{METADATA_INPUT}\n\n"
            "Run scripts/investigate_gse94533_labels.py first."
        )

    meta = pd.read_csv(
        METADATA_INPUT,
        dtype=str,
    )

    required = {
        "geo_accession",
        "title",
        "diagnosis",
    }

    missing = required - set(meta.columns)

    if missing:

        raise RuntimeError(
            f"GEO metadata is missing required columns: {missing}"
        )

    meta = meta.copy()

    meta["sample_id_normalized"] = meta["title"].map(
        normalize_sample_id
    )

    # ------------------------------------------------------------
    # Normalize diagnosis
    # ------------------------------------------------------------

    def normalize_diagnosis(value: object) -> str:

        if pd.isna(value):
            return ""

        s = str(value).strip().lower()

        if s == "cancer":
            return "Cancer"

        if s == "benign":
            return "Benign"

        if s == "borderline":
            return "Borderline"

        if s == "control":
            return "Control"

        return str(value).strip()

    meta["diagnosis"] = meta["diagnosis"].map(
        normalize_diagnosis
    )

    # ------------------------------------------------------------
    # Binary endpoint
    # ------------------------------------------------------------

    valid = {
        "Cancer",
        "Benign",
        "Borderline",
        "Control",
    }

    bad = sorted(
        set(meta["diagnosis"].dropna())
        - valid
    )

    if bad:

        raise RuntimeError(
            "Unexpected GEO diagnosis values:\n"
            + "\n".join(map(str, bad))
        )

    meta["binary_label"] = (
        meta["diagnosis"] == "Cancer"
    ).astype(int)

    # ------------------------------------------------------------
    # Explicit exclusion
    # ------------------------------------------------------------

    meta["excluded_from_analysis"] = (
        meta["sample_id_normalized"] == "e1056"
    )

    return meta


# ================================================================
# FIND SAMPLE COLUMNS IN EXPRESSION WORKBOOK
# ================================================================


def find_sample_columns(
    workbook: Path,
    valid_ids: set[str],
) -> list[tuple[object, str]]:

    """
    Identify expression-matrix sample columns without changing the
    original pandas column type.

    This matters because Cohort 1 uses identifiers such as:

        e1036

    while Cohort 2 uses numeric Excel headers such as:

        50217.0

    We must preserve the original column object for dataframe indexing
    and only use the normalized string for matching to GEO metadata.

    Returns
    -------
    list of:
        (original_excel_column_object, normalized_sample_id)
    """

    header = pd.read_excel(
        workbook,
        sheet_name="mirna_tpm",
        header=2,
        nrows=0,
        engine="openpyxl",
    )

    results = []

    for column in header.columns:

        normalized = normalize_sample_id(column)

        if normalized in valid_ids:

            # IMPORTANT:
            # Keep `column` exactly as pandas read it.
            #
            # Do NOT use str(column) here.
            #
            results.append(
                (column, normalized)
            )

    return results


def read_expression_workbook(
    workbook: Path,
    valid_ids: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:

    print()
    print("=" * 90)
    print(f"READING {workbook.name}")
    print("=" * 90)

    sample_columns = find_sample_columns(
        workbook,
        valid_ids,
    )

    if not sample_columns:

        raise RuntimeError(
            f"No GEO sample IDs matched {workbook.name}."
        )

    print()
    print(f"Matched sample columns: {len(sample_columns)}")

    # ------------------------------------------------------------
    # Read the actual matrix.
    # ------------------------------------------------------------

    df = pd.read_excel(
        workbook,
        sheet_name="mirna_tpm",
        header=2,
        engine="openpyxl",
    )

    columns = list(df.columns)

    feature_column = columns[0]

    print(
        f"Feature column: {feature_column!r}"
    )

    # ------------------------------------------------------------
    # Select sample columns using the ORIGINAL pandas objects.
    # ------------------------------------------------------------
    #
    # Cohort 2 contains headers such as 50217.0.
    #
    # We must select:
    #
    #     50217.0
    #
    # rather than:
    #
    #     "50217"
    #
    # because those are different pandas column labels.

    selected_columns = [
        original
        for original, _ in sample_columns
    ]

    expression = df[
        [feature_column] + selected_columns
    ].copy()

    expression = expression.rename(
        columns={
            feature_column: "feature_name"
        }
    )

    # ------------------------------------------------------------
    # Clean feature names.
    # ------------------------------------------------------------

    expression["feature_name"] = (
        expression["feature_name"]
        .astype("string")
        .str.strip()
    )

    expression = expression[
        expression["feature_name"].notna()
        & (expression["feature_name"] != "")
    ].copy()

    # ------------------------------------------------------------
    # Convert expression values to numeric.
    # ------------------------------------------------------------

    for column in selected_columns:

        expression[column] = pd.to_numeric(
            expression[column],
            errors="coerce",
        )

    expression[selected_columns] = (
        expression[selected_columns]
        .replace([float("inf"), float("-inf")], float("nan"))
        .fillna(0.0)
    )

    # ------------------------------------------------------------
    # Rename sample columns AFTER selecting them.
    # ------------------------------------------------------------

    rename_map = {}

    mapping_rows = []

    for original, normalized in sample_columns:

        rename_map[original] = normalized

        mapping_rows.append(
            {
                "source_workbook": workbook.name,
                "excel_column": str(original),
                "sample_id_normalized": normalized,
            }
        )

    expression = expression.rename(
        columns=rename_map
    )

    normalized_sample_columns = [
        normalized
        for _, normalized in sample_columns
    ]

    # ------------------------------------------------------------
    # Aggregate duplicate miRNA feature names.
    # ------------------------------------------------------------

    grouped = (
        expression
        .groupby(
            "feature_name",
            sort=True,
        )[normalized_sample_columns]
        .sum()
    )

    print(
        f"Features after aggregation: {grouped.shape[0]}"
    )

    # ------------------------------------------------------------
    # Transpose:
    #
    # rows    = patients
    # columns = miRNAs
    # ------------------------------------------------------------

    patient_matrix = grouped.T.copy()

    patient_matrix.index = [
        normalize_sample_id(x)
        for x in patient_matrix.index
    ]

    patient_matrix.index.name = (
        "sample_id_normalized"
    )

    mapping = pd.DataFrame(
        mapping_rows
    )

    # ------------------------------------------------------------
    # Safety checks.
    # ------------------------------------------------------------

    if len(patient_matrix.index) != len(set(patient_matrix.index)):

        raise RuntimeError(
            f"Duplicate normalized patient IDs detected in "
            f"{workbook.name} after renaming."
        )

    return patient_matrix, mapping

def build_expression_matrix(
    metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:

    workbooks = find_workbooks()

    valid_geo_ids = set(
        metadata["sample_id_normalized"]
    )

    matrices = []
    mappings = []

    for workbook in workbooks:

        matrix, mapping = read_expression_workbook(
            workbook=workbook,
            valid_ids=valid_geo_ids,
        )

        matrices.append(matrix)
        mappings.append(mapping)

    # ------------------------------------------------------------
    # Make sure sample IDs aren't duplicated across workbooks.
    # ------------------------------------------------------------

    all_ids = [
        sample_id
        for matrix in matrices
        for sample_id in matrix.index
    ]

    duplicates = sorted(
        {
            x
            for x in all_ids
            if all_ids.count(x) > 1
        }
    )

    if duplicates:

        raise RuntimeError(
            "Duplicate patients across expression workbooks:\n"
            + "\n".join(duplicates)
        )

    # ------------------------------------------------------------
    # Combine all patients.
    # ------------------------------------------------------------

    combined = pd.concat(
        matrices,
        axis=0,
        join="outer",
    )

    combined = combined.fillna(0.0)

    # ------------------------------------------------------------
    # Remove E-1056.
    # ------------------------------------------------------------

    if "e1056" in combined.index:

        print()
        print(
            "Removing E-1056 according to the original study."
        )

        combined = combined.drop(
            index="e1056"
        )

    # ------------------------------------------------------------
    # Restrict GEO metadata to expression patients.
    # ------------------------------------------------------------

    metadata = metadata[
        ~metadata["excluded_from_analysis"]
    ].copy()

    metadata_ids = set(
        metadata["sample_id_normalized"]
    )

    expression_ids = set(
        combined.index
    )

    missing_expression = sorted(
        metadata_ids - expression_ids
    )

    missing_metadata = sorted(
        expression_ids - metadata_ids
    )

    if missing_expression:

        print()
        print(
            "WARNING: GEO metadata samples missing from expression:"
        )

        for sample_id in missing_expression:
            print(" ", sample_id)

    if missing_metadata:

        print()
        print(
            "WARNING: Expression samples missing from GEO metadata:"
        )

        for sample_id in missing_metadata:
            print(" ", sample_id)

    # Only retain patients for which we have both expression and labels.
    common = [
        sample_id
        for sample_id in metadata[
            "sample_id_normalized"
        ]
        if sample_id in combined.index
    ]

    metadata = (
        metadata
        .set_index("sample_id_normalized")
        .loc[common]
        .reset_index()
    )

    combined = combined.loc[
        common
    ]

    return combined, metadata


# ================================================================
# DATASET AUDIT
# ================================================================

def audit_dataset(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
) -> dict:

    print()
    print("=" * 90)
    print("FINAL DATASET AUDIT")
    print("=" * 90)

    n_patients = int(matrix.shape[0])
    n_features = int(matrix.shape[1])

    diagnosis_counts = (
        metadata["diagnosis"]
        .value_counts()
        .to_dict()
    )

    cancer_count = int(
        (metadata["binary_label"] == 1).sum()
    )

    noncancer_count = int(
        (metadata["binary_label"] == 0).sum()
    )

    finite = bool(
        np.isfinite(
            matrix.to_numpy(dtype=float)
        ).all()
    )

    nonnegative = bool(
        (
            matrix.to_numpy(dtype=float) >= 0
        ).all()
    )

    unique_patients = (
        not matrix.index.duplicated().any()
    )

    aligned = (
        list(matrix.index)
        == list(metadata["sample_id_normalized"])
    )

    audit = {
        "dataset": "GSE94533",
        "n_patients": n_patients,
        "n_mirna_features": n_features,
        "diagnosis_counts": diagnosis_counts,
        "binary_counts": {
            "cancer": cancer_count,
            "non_cancer": noncancer_count,
        },
        "excluded_samples": [
            "E-1056",
        ],
        "matrix_checks": {
            "all_finite": finite,
            "all_nonnegative": nonnegative,
            "unique_patient_ids": unique_patients,
            "metadata_expression_aligned": aligned,
        },
        "endpoint_definition": {
            "1": "GEO diagnosis = Cancer",
            "0": (
                "GEO diagnosis = Benign, Borderline, or Control"
            ),
        },
        "important_provenance_note": (
            "Current GEO diagnosis annotations are preserved. "
            "They are not overwritten to force agreement with "
            "the diagnosis counts reported in the 2017 eLife paper."
        ),
    }

    print()
    print(f"Patients:          {n_patients}")
    print(f"miRNA features:    {n_features}")

    print()
    print("Diagnosis counts:")
    for key, value in diagnosis_counts.items():
        print(
            f"  {key:12s} {value}"
        )

    print()
    print("Binary counts:")
    print(
        f"  Cancer:        {cancer_count}"
    )
    print(
        f"  Non-cancer:    {noncancer_count}"
    )

    print()
    print("Matrix checks:")
    print(
        f"  Finite values:        {finite}"
    )
    print(
        f"  Non-negative values:   {nonnegative}"
    )
    print(
        f"  Unique patient IDs:    {unique_patients}"
    )
    print(
        f"  Metadata aligned:      {aligned}"
    )

    # ------------------------------------------------------------
    # Hard integrity requirements.
    # ------------------------------------------------------------

    if n_patients != 179:

        raise RuntimeError(
            f"Expected 179 eligible patients after E-1056 exclusion, "
            f"but found {n_patients}."
        )

    if not finite:

        raise RuntimeError(
            "Expression matrix contains NaN or infinite values."
        )

    if not nonnegative:

        raise RuntimeError(
            "Expression matrix contains negative TPM values."
        )

    if not unique_patients:

        raise RuntimeError(
            "Duplicate patient IDs remain."
        )

    if not aligned:

        raise RuntimeError(
            "Metadata and expression matrix are not aligned."
        )

    # We intentionally do NOT impose the historical 98/21/45/15
    # count constraint here, because the current GEO metadata are
    # the labels being used for this dataset.

    return audit


# ================================================================
# SAVE DATA
# ================================================================

def save_outputs(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    audit: dict,
) -> None:

    INTERIM_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    PROCESSED_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    REPORT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ------------------------------------------------------------
    # Metadata
    # ------------------------------------------------------------

    metadata.to_csv(
        METADATA_OUTPUT,
        index=False,
    )

    # ------------------------------------------------------------
    # Matrix
    # ------------------------------------------------------------

    matrix.to_csv(
        MATRIX_OUTPUT,
        compression="gzip",
    )

    # ------------------------------------------------------------
    # Audit
    # ------------------------------------------------------------

    AUDIT_OUTPUT.write_text(
        json.dumps(
            audit,
            indent=2,
        )
    )

    print()
    print("=" * 90)
    print("OUTPUT FILES")
    print("=" * 90)

    print()
    print(
        f"Metadata:\n  {METADATA_OUTPUT}"
    )

    print()
    print(
        f"Expression matrix:\n  {MATRIX_OUTPUT}"
    )

    print()
    print(
        f"Audit:\n  {AUDIT_OUTPUT}"
    )


# ================================================================
# MAIN
# ================================================================

def main() -> None:

    print("=" * 90)
    print("GSE94533 CLEAN PATIENT-LEVEL DATASET BUILDER")
    print("=" * 90)

    # ------------------------------------------------------------
    # Load metadata.
    # ------------------------------------------------------------

    metadata = load_metadata()

    print()
    print(
        f"GEO metadata records: {len(metadata)}"
    )

    print()
    print("Current GEO diagnosis distribution:")

    print(
        metadata["diagnosis"]
        .value_counts()
        .to_string()
    )

    # ------------------------------------------------------------
    # Build expression matrix.
    # ------------------------------------------------------------

    matrix, metadata = build_expression_matrix(
        metadata
    )

    # ------------------------------------------------------------
    # Audit.
    # ------------------------------------------------------------

    audit = audit_dataset(
        matrix=matrix,
        metadata=metadata,
    )

    # ------------------------------------------------------------
    # Save.
    # ------------------------------------------------------------

    save_outputs(
        matrix=matrix,
        metadata=metadata,
        audit=audit,
    )

    print()
    print("=" * 90)
    print("SUCCESS")
    print("=" * 90)

    print(
        """
The clean patient-level GSE94533 dataset is now ready.

No machine-learning model has been trained.

The dataset now has:

    patient
       |
       +---- GEO diagnosis
       +---- binary cancer endpoint
       +---- stage
       +---- histology
       |
       +---- serum miRNA TPM features

Next we can build:

    1. A faithful 14-miRNA historical benchmark
    2. A modern regularized miRNA baseline
    3. Stage-I/II analysis
    4. Benign-vs-cancer analysis
    5. Later, the multimodal architecture
    """
    )


if __name__ == "__main__":
    main()
