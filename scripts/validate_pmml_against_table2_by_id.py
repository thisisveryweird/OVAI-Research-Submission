#!/usr/bin/env python3

from __future__ import annotations

import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd
from sklearn.metrics import accuracy_score, confusion_matrix


ROOT = Path.home() / "ovarian-ai"

PREDICTIONS = (
    ROOT
    / "reports"
    / "tables"
    / "original_historical_clean_predictions.csv"
)

DOCX = (
    ROOT
    / "data"
    / "external"
    / "elife_28932"
    / "elife_poa_e28932_Supplementary_File_1.docx"
)

OUTPUT = (
    ROOT
    / "reports"
    / "tables"
    / "pmml_vs_table2_by_actual_sample_id.csv"
)


# ================================================================
# NORMALIZATION
# ================================================================

def norm_id(value) -> str:

    s = str(value).strip().lower()

    if re.fullmatch(r"\d+\.0+", s):
        s = s.split(".", 1)[0]

    return re.sub(r"[^a-z0-9]", "", s)


# ================================================================
# DOCX TABLE EXTRACTION
# ================================================================

def cell_text(cell) -> str:

    parts = []

    for node in cell.iter():

        if node.tag.endswith("}t"):
            parts.append(node.text or "")

    return "".join(parts).strip()


def extract_docx_tables(path: Path):

    with zipfile.ZipFile(path) as z:

        xml = z.read("word/document.xml")

    root = ET.fromstring(xml)

    tables = []

    for table in root.iter():

        if not table.tag.endswith("}tbl"):
            continue

        rows = []

        for row in table:

            if not row.tag.endswith("}tr"):
                continue

            cells = []

            for cell in row:

                if not cell.tag.endswith("}tc"):
                    continue

                cells.append(cell_text(cell))

            if cells:
                rows.append(cells)

        if rows:
            tables.append(rows)

    return tables


# ================================================================
# FIND TABLE 2
# ================================================================

def find_table2(tables):

    for i, table in enumerate(tables):

        if not table:
            continue

        header = [
            str(x).strip().lower()
            for x in table[0]
        ]

        has_id = "id" in header
        has_nn = "nn" in header

        if has_id and has_nn:

            return i + 1, table

    raise RuntimeError(
        "Could not identify published Table 2."
    )


# ================================================================
# LOAD PREDICTIONS
# ================================================================

pred = pd.read_csv(PREDICTIONS)

if len(pred) != 179:

    raise RuntimeError(
        f"Prediction file contains {len(pred)} patients, "
        "expected 179."
    )

pred["sample_id"] = pred["sample_id"].map(norm_id)

if pred["sample_id"].duplicated().any():

    raise RuntimeError(
        "Prediction file contains duplicate sample IDs."
    )


# ================================================================
# BUILD TABLE-2-ID -> ACTUAL SAMPLE-ID MAPPING
# ================================================================

def table2_id_to_actual_sample_id(table2_id: str) -> str:

    sid = norm_id(table2_id)

    # ------------------------------------------------------------
    # Cohort 2
    #
    # Published Table 2 uses patient numbers:
    #
    #     1 ... 55
    #     57 ... 120
    #
    # Patient 56 was excluded.
    #
    # The prediction file is already ordered as:
    #
    #     Table2 1 ... 55
    #     Table2 57 ... 120
    #     e-samples
    #
    # Therefore:
    #
    #     Table2 ID <=55
    #         -> prediction row ID-1
    #
    #     Table2 ID >=57
    #         -> prediction row ID-2
    # ------------------------------------------------------------

    if re.fullmatch(r"\d+", sid):

        number = int(sid)

        if 1 <= number <= 55:

            return norm_id(
                pred.iloc[number - 1]["sample_id"]
            )

        if 57 <= number <= 120:

            return norm_id(
                pred.iloc[number - 2]["sample_id"]
            )

        raise RuntimeError(
            f"Unexpected numeric Table 2 patient ID: {sid}"
        )

    # ------------------------------------------------------------
    # ERASMOS e-samples use their actual identifiers:
    #
    # e1001
    # e1002
    # ...
    # ------------------------------------------------------------

    if sid.startswith("e"):

        return sid

    raise RuntimeError(
        f"Unrecognized Table 2 patient ID: {table2_id!r}"
    )


# ================================================================
# READ TABLE 2
# ================================================================

print("=" * 90)
print("PMML VS PUBLISHED TABLE 2 — ACTUAL SAMPLE-ID VALIDATION")
print("=" * 90)

print()
print("Prediction patients:", len(pred))

tables = extract_docx_tables(DOCX)

table_number, table2 = find_table2(tables)

print()
print("Tables extracted:", len(tables))
print("Identified Table:", table_number)

header = [
    str(x).strip()
    for x in table2[0]
]

print()
print("Table 2 header:")
for i, value in enumerate(header):
    print(f"  {i}: {value}")


# ================================================================
# CONVERT TABLE 2 TO DATAFRAME
# ================================================================

rows = []

for row in table2[1:]:

    if len(row) < len(header):
        row = row + [""] * (
            len(header) - len(row)
        )

    if len(row) > len(header):
        row = row[:len(header)]

    if not str(row[0]).strip():
        continue

    rows.append(row)


published = pd.DataFrame(
    rows,
    columns=header,
)

published = published.iloc[:179].copy()

if len(published) != 179:

    raise RuntimeError(
        f"Published Table 2 contains "
        f"{len(published)} patient rows instead of 179."
    )

published["table2_id"] = (
    published["ID"]
    .map(norm_id)
)

published["actual_sample_id"] = (
    published["table2_id"]
    .map(table2_id_to_actual_sample_id)
)


# ================================================================
# VALIDATE MAPPING
# ================================================================

print()
print("=" * 90)
print("TABLE 2 -> ACTUAL SAMPLE-ID MAPPING")
print("=" * 90)

print()

print(
    published[
        [
            "table2_id",
            "actual_sample_id",
        ]
    ].head(10).to_string(index=False)
)

print()

print(
    published[
        [
            "table2_id",
            "actual_sample_id",
        ]
    ].tail(10).to_string(index=False)
)


if (
    published["actual_sample_id"]
    .duplicated()
    .any()
):

    duplicates = published.loc[
        published["actual_sample_id"].duplicated(
            keep=False
        ),
        [
            "table2_id",
            "actual_sample_id",
        ],
    ]

    print()
    print(duplicates.to_string(index=False))

    raise RuntimeError(
        "Table 2 mapping generated duplicate actual sample IDs."
    )


pred_ids = set(pred["sample_id"])
table_ids = set(
    published["actual_sample_id"]
)

missing_from_predictions = (
    table_ids - pred_ids
)

extra_predictions = (
    pred_ids - table_ids
)

print()
print("=" * 90)
print("ACTUAL SAMPLE-ID SET VALIDATION")
print("=" * 90)

print()
print("Prediction IDs:", len(pred_ids))
print("Table 2 IDs:  ", len(table_ids))

print()
print(
    "Missing from predictions:",
    len(missing_from_predictions),
)

if missing_from_predictions:
    print(
        sorted(missing_from_predictions)
    )

print()
print(
    "Extra prediction IDs:",
    len(extra_predictions),
)

if extra_predictions:
    print(
        sorted(extra_predictions)
    )

if missing_from_predictions or extra_predictions:

    raise RuntimeError(
        "Prediction and Table 2 sample-ID sets do not match."
    )


# ================================================================
# JOIN BY ACTUAL SAMPLE ID
# ================================================================

comparison = published.merge(
    pred[
        [
            "sample_id",
            "predicted_class",
            "cancer_activation",
            "controls_or_borderline_activation",
        ]
    ],
    left_on="actual_sample_id",
    right_on="sample_id",
    how="inner",
    validate="one_to_one",
)


if len(comparison) != 179:

    raise RuntimeError(
        f"Joined comparison contains {len(comparison)} "
        "rows instead of 179."
    )


# ================================================================
# NORMALIZE PUBLISHED NN LABEL
# ================================================================

comparison["published_nn"] = (
    comparison["NN"]
    .astype(str)
    .str.strip()
    .str.upper()
)

comparison["pmml_nn"] = comparison[
    "predicted_class"
].map(
    {
        "Ovarian Cancer": "CA",
        "Controls or Borderline": "BBC",
    }
)


comparison["exact_match"] = (
    comparison["published_nn"]
    == comparison["pmml_nn"]
)


# ================================================================
# RESULTS
# ================================================================

matches = int(
    comparison["exact_match"].sum()
)

mismatches = (
    len(comparison) - matches
)

agreement = (
    matches / len(comparison)
)

print()
print("=" * 90)
print("FINAL RESULTS")
print("=" * 90)

print()
print(f"Patients:            {len(comparison)}")
print(f"Exact matches:       {matches}")
print(f"Mismatches:          {mismatches}")
print(f"Prediction agreement:{agreement:.6f}")


# ================================================================
# CONFUSION MATRIX
# ================================================================

cm = confusion_matrix(
    comparison["published_nn"],
    comparison["pmml_nn"],
    labels=["BBC", "CA"],
)

print()
print("Confusion matrix:")
print()
print("                    PMML")
print("                 BBC     CA")
print(
    f"Published BBC    "
    f"{cm[0,0]:3d}   {cm[0,1]:3d}"
)
print(
    f"Published CA     "
    f"{cm[1,0]:3d}   {cm[1,1]:3d}"
)


# ================================================================
# MISMATCHES
# ================================================================

if mismatches:

    print()
    print("=" * 90)
    print("MISMATCHED PATIENTS")
    print("=" * 90)

    mismatch_cols = [
        "table2_id",
        "actual_sample_id",
        "published_nn",
        "pmml_nn",
        "cancer_activation",
        "controls_or_borderline_activation",
    ]

    print()

    print(
        comparison.loc[
            ~comparison["exact_match"],
            mismatch_cols,
        ].to_string(index=False)
    )

else:

    print()
    print("=" * 90)
    print("PERFECT PUBLISHED-PMML AGREEMENT")
    print("=" * 90)

    print()
    print(
        "Every one of the 179 published Table 2 "
        "predictions matches the PMML forward pass."
    )


# ================================================================
# SAVE
# ================================================================

comparison.to_csv(
    OUTPUT,
    index=False,
)

print()
print("Saved:")
print(f"  {OUTPUT}")


# ================================================================
# FINAL INTERPRETATION
# ================================================================

print()
print("=" * 90)

if mismatches == 0:

    print("HISTORICAL MODEL VALIDATION: PASSED")

else:

    print(
        "HISTORICAL MODEL VALIDATION: "
        f"{mismatches} MISMATCHES"
    )

print("=" * 90)

