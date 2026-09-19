#!/usr/bin/env python3

from pathlib import Path
import csv
import gzip
import re

import numpy as np
import pandas as pd


ROOT = Path.home() / "ovarian-ai"

INPUT = (
    ROOT
    / "data"
    / "external"
    / "gse8841"
    / "GSE8841_series_matrix.txt.gz"
)

OUTPUT = (
    ROOT
    / "data"
    / "processed"
    / "gse8841_borderline_expression.csv.gz"
)

print("=" * 80)
print("BUILDING GSE8841 BORDERLINE EXPRESSION DATASET")
print("=" * 80)

# ------------------------------------------------------------------
# Read sample metadata without losing repeated characteristics rows
# ------------------------------------------------------------------

metadata = {}

with gzip.open(
    INPUT,
    "rt",
    encoding="utf-8",
    errors="replace",
) as f:

    for line in f:

        if line.startswith(
            "!series_matrix_table_begin"
        ):
            break

        if line.startswith("!Sample_"):

            parts = next(
                csv.reader(
                    [line.rstrip("\n")],
                    delimiter="\t",
                    quotechar='"',
                )
            )

            metadata.setdefault(
                parts[0],
                []
            ).append(parts[1:])


accessions = metadata[
    "!Sample_geo_accession"
][0]

titles = metadata[
    "!Sample_title"
][0]

char_rows = metadata[
    "!Sample_characteristics_ch1"
]

n_samples = len(accessions)

print()
print("Samples:", n_samples)

# ------------------------------------------------------------------
# Build per-sample metadata
# ------------------------------------------------------------------

sample_info = [
    {}
    for _ in range(n_samples)
]

for row in char_rows:

    for i, value in enumerate(row):

        if ":" not in value:
            continue

        key, val = value.split(
            ":",
            1,
        )

        sample_info[i][
            key.strip().lower()
        ] = val.strip()


records = []

for i in range(n_samples):

    info = sample_info[i]

    stage = info.get(
        "figo stage",
        "",
    )

    substage = info.get(
        "figo substage",
        "",
    )

    grade = info.get(
        "grade",
        "",
    )

    age = info.get(
        "age(years)",
        "",
    )

    os = info.get(
        "overall survival(days)",
        "",
    )

    pfs = info.get(
        "progression free survival(days)",
        "",
    )

    histotype = info.get(
        "histotype",
        "",
    )

    relapse = info.get(
        "relapsed",
        "",
    )

    status = info.get(
        "status",
        "",
    )

    # --------------------------------------------------------------
    # Define endpoint.
    #
    # Borderline is explicitly defined by Grade:Borderline.
    # Cancer is FIGO Stage:1 with non-borderline grade.
    # --------------------------------------------------------------

    if grade.strip().lower() == "borderline":

        label = "borderline"

    elif stage.strip() == "1":

        label = "stage_I_cancer"

    else:

        raise RuntimeError(
            f"Unexpected sample classification: "
            f"{accessions[i]} | {info}"
        )

    records.append({
        "sample_id": accessions[i],
        "title": titles[i],
        "label": label,
        "figo_stage": stage,
        "figo_substage": substage,
        "histotype": histotype,
        "grade": grade,
        "age": pd.to_numeric(
            age,
            errors="coerce",
        ),
        "overall_survival_days": pd.to_numeric(
            os,
            errors="coerce",
        ),
        "progression_free_survival_days": pd.to_numeric(
            pfs,
            errors="coerce",
        ),
        "relapsed": relapse,
        "status": status,
    })


clinical = pd.DataFrame(records)

# ------------------------------------------------------------------
# Parse expression matrix
# ------------------------------------------------------------------

print()
print("=" * 80)
print("READING EXPRESSION MATRIX")
print("=" * 80)

matrix_lines = []

with gzip.open(
    INPUT,
    "rt",
    encoding="utf-8",
    errors="replace",
) as f:

    inside = False

    for line in f:

        if line.startswith(
            "!series_matrix_table_begin"
        ):
            inside = True
            continue

        if line.startswith(
            "!series_matrix_table_end"
        ):
            break

        if inside and line.strip():
            matrix_lines.append(
                line
            )

if not matrix_lines:
    raise RuntimeError(
        "Expression matrix was not found."
    )

from io import StringIO

expr = pd.read_csv(
    StringIO(
        "".join(matrix_lines)
    ),
    sep="\t",
    quotechar='"',
)

print(
    "Expression rows:",
    len(expr)
)

print(
    "Expression columns:",
    len(expr.columns)
)

# ------------------------------------------------------------------
# Verify sample columns
# ------------------------------------------------------------------

sample_cols = [
    c
    for c in expr.columns
    if c in set(accessions)
]

if len(sample_cols) != n_samples:

    missing = sorted(
        set(accessions)
        - set(sample_cols)
    )

    raise RuntimeError(
        "Expression/sample mismatch. "
        f"Missing sample columns: {missing}"
    )

print(
    "Matched sample columns:",
    len(sample_cols)
)

# ------------------------------------------------------------------
# Convert expression values to numeric
# ------------------------------------------------------------------

feature_id_col = expr.columns[0]

expr[feature_id_col] = (
    expr[feature_id_col]
    .astype(str)
)

for c in sample_cols:

    expr[c] = pd.to_numeric(
        expr[c],
        errors="coerce",
    )

# ------------------------------------------------------------------
# Transpose: patients x expression features
# ------------------------------------------------------------------

expr_patient = expr[
    [feature_id_col] + sample_cols
].set_index(
    feature_id_col
).T.reset_index()

expr_patient = expr_patient.rename(
    columns={
        "index": "sample_id"
    }
)

expr_patient["sample_id"] = (
    expr_patient["sample_id"]
    .astype(str)
)

# ------------------------------------------------------------------
# Remove duplicate probe IDs if any
# ------------------------------------------------------------------

if expr_patient["sample_id"].duplicated().any():

    raise RuntimeError(
        "Duplicate sample IDs found."
    )

# ------------------------------------------------------------------
# Merge clinical + expression
# ------------------------------------------------------------------

df = clinical.merge(
    expr_patient,
    on="sample_id",
    how="inner",
    validate="one_to_one",
)

# ------------------------------------------------------------------
# Validate
# ------------------------------------------------------------------

expression_cols = [
    c
    for c in expr_patient.columns
    if c != "sample_id"
]

print()
print("=" * 80)
print("FINAL DATASET VALIDATION")
print("=" * 80)

print(
    "Patients:",
    len(df)
)

print(
    "Stage-I cancer:",
    int(
        (
            df["label"]
            == "stage_I_cancer"
        ).sum()
    )
)

print(
    "Borderline:",
    int(
        (
            df["label"]
            == "borderline"
        ).sum()
    )
)

print(
    "Expression features:",
    len(expression_cols)
)

print(
    "Age observed:",
    int(
        df["age"].notna().sum()
    )
)

print(
    "Expression missing values:",
    int(
        df[expression_cols]
        .isna()
        .sum()
        .sum()
    )
)

print()
print("Histotype counts:")
print(
    pd.crosstab(
        df["label"],
        df["histotype"],
    ).to_string()
)

# ------------------------------------------------------------------
# Save
# ------------------------------------------------------------------

OUTPUT.parent.mkdir(
    parents=True,
    exist_ok=True,
)

df.to_csv(
    OUTPUT,
    index=False,
    compression="gzip",
)

print()
print("=" * 80)
print("SAVED")
print("=" * 80)

print(OUTPUT)

print()
print("IMPORTANT:")
print(
    "FIGO stage and grade are retained as metadata,"
)
print(
    "but must NOT be used as predictive features."
)
print(
    "This dataset is tissue-expression data,"
)
print(
    "not serum-miRNA data."
)
