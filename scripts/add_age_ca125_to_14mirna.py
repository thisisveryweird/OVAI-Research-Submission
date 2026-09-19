#!/usr/bin/env python3

from pathlib import Path
import re
import zipfile
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd


ROOT = Path.home() / "ovarian-ai"

BASE = (
    ROOT
    / "data"
    / "processed"
    / "gse94533_14mirna_benchmark.csv.gz"
)

PRED = (
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

OUT = (
    ROOT
    / "data"
    / "processed"
    / "gse94533_14mirna_clinical_benchmark.csv.gz"
)


def norm_id(value):
    s = str(value).strip().lower()

    if re.fullmatch(r"\d+\.0+", s):
        s = s.split(".", 1)[0]

    return re.sub(r"[^a-z0-9]", "", s)


def cell_text(cell):
    parts = []

    for node in cell.iter():
        if node.tag.endswith("}t"):
            parts.append(node.text or "")

    return "".join(parts).strip()


def extract_tables(path):
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml")

    root = ET.fromstring(xml)
    tables = []

    for tbl in root.iter():

        if not tbl.tag.endswith("}tbl"):
            continue

        rows = []

        for tr in tbl:

            if not tr.tag.endswith("}tr"):
                continue

            row = []

            for tc in tr:

                if tc.tag.endswith("}tc"):
                    row.append(
                        cell_text(tc)
                    )

            if row:
                rows.append(row)

        if rows:
            tables.append(rows)

    return tables


def parse_number(value):

    s = str(value).strip()

    if s in {"", "*", "‡", "nan"}:
        return np.nan

    try:
        return float(s)
    except ValueError:
        return np.nan


print("=" * 80)
print("ADDING AGE + CA-125 TO 14-miRNA BENCHMARK")
print("=" * 80)

base = pd.read_csv(BASE)
pred = pd.read_csv(PRED)

if len(base) != 179:
    raise RuntimeError(
        f"Expected 179 benchmark patients, found {len(base)}"
    )

if len(pred) != 179:
    raise RuntimeError(
        f"Expected 179 prediction patients, found {len(pred)}"
    )

pred["sample_key"] = pred["sample_id"].map(norm_id)

tables = extract_tables(DOCX)

if len(tables) < 2:
    raise RuntimeError("Could not find Table 2.")

table2 = tables[1]
header = table2[0]

id_col = header.index("ID")
age_col = header.index("Age")
ca_col = header.index("CA-125")

records = []

for row in table2[1:]:

    if len(row) <= max(
        id_col,
        age_col,
        ca_col,
    ):
        continue

    table_id = row[id_col].strip()

    if not table_id:
        continue

    # ---------------------------------------------------------------
    # Map Table 2 ID -> actual sample ID.
    # ---------------------------------------------------------------

    key = norm_id(table_id)

    if key.startswith("e"):

        actual_key = key

    else:

        position = int(float(table_id))

        if position == 56:
            # This patient is excluded from the 179-patient cohort.
            continue

        if position < 1 or position > 120:
            raise RuntimeError(
                f"Unexpected Table 2 numeric ID: {position}"
            )

        retained_index = (
            position - 1
            if position <= 55
            else position - 2
        )

        actual_key = pred.iloc[
            retained_index
        ]["sample_key"]

    records.append(
        {
            "sample_key": actual_key,
            "age": parse_number(row[age_col]),
            "ca125": parse_number(row[ca_col]),
        }
    )


clinical = pd.DataFrame(records)

clinical = clinical.drop_duplicates(
    subset="sample_key",
    keep="first",
)


if len(clinical) != 179:
    raise RuntimeError(
        f"Expected 179 clinical records, got {len(clinical)}"
    )


# ---------------------------------------------------------------
# Join onto the existing 14-miRNA benchmark.
# ---------------------------------------------------------------

base["sample_key"] = base["sample_id"].map(norm_id)

merged = base.merge(
    clinical,
    on="sample_key",
    how="left",
    validate="one_to_one",
)

if len(merged) != 179:
    raise RuntimeError(
        f"Expected 179 merged patients, got {len(merged)}"
    )


# ---------------------------------------------------------------
# Validation.
# ---------------------------------------------------------------

print()
print("Patients:", len(merged))

print(
    "Age observed:",
    int(merged["age"].notna().sum()),
)

print(
    "Age missing:",
    int(merged["age"].isna().sum()),
)

print(
    "CA-125 observed:",
    int(merged["ca125"].notna().sum()),
)

print(
    "CA-125 missing:",
    int(merged["ca125"].isna().sum()),
)

print()
print("Age summary:")
print(
    merged["age"].describe().to_string()
)

print()
print("CA-125 summary:")
print(
    merged["ca125"].describe().to_string()
)

# ---------------------------------------------------------------
# Keep the useful columns only.
# ---------------------------------------------------------------

keep = [
    "sample_id",
    "diagnosis",
    "binary_label",
    *[
        c
        for c in base.columns
        if c.startswith("hsa-")
    ],
    "age",
    "ca125",
]

merged[keep].to_csv(
    OUT,
    index=False,
    compression="gzip",
)

print()
print("=" * 80)
print("CREATED")
print("=" * 80)

print(OUT)
