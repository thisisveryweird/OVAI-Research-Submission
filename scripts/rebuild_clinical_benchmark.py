import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

ROOT = Path.home() / "ovarian-ai"

PMML = ROOT / "reports/tables/original_historical_clean_predictions.csv"

XLSX = ROOT / "data/external/elife_28932/elife_poa_e28932_Supplementary_File_6.xlsx"

DOCX = ROOT / "data/external/elife_28932/elife_poa_e28932_Supplementary_File_1.docx"

OUT = ROOT / "data/processed/gse94533_14mirna_clinical_benchmark.csv.gz"

FEATURES = [
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


def norm(x):
    s = str(x).strip().lower()
    if re.fullmatch(r"\d+\.0+", s):
        s = s.split(".", 1)[0]
    return re.sub(r"[^a-z0-9]", "", s)


def cell_text(cell):
    return "".join(
        node.text or ""
        for node in cell.iter()
        if node.tag.endswith("}t")
    ).strip()


def tables_from_docx(path):
    with zipfile.ZipFile(path) as z:
        root = ET.fromstring(z.read("word/document.xml"))

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
                    row.append(cell_text(tc))

            if row:
                rows.append(row)

        if rows:
            tables.append(rows)

    return tables


def number(x):
    s = str(x).strip()

    if s in {"", "*", "‡", "nan"}:
        return np.nan

    try:
        return float(s)
    except ValueError:
        return np.nan


print("=" * 80)
print("REBUILDING VERIFIED CLINICAL BENCHMARK")
print("=" * 80)

# ---------------------------------------------------------------
# 1. VERIFIED 179-PATIENT COHORT
# ---------------------------------------------------------------

pmml = pd.read_csv(PMML)

if len(pmml) != 179:
    raise RuntimeError(
        f"PMML prediction file has {len(pmml)} patients, expected 179."
    )

verified_ids = pmml["sample_id"].astype(str).tolist()
verified_keys = [norm(x) for x in verified_ids]

if "81195" in verified_keys:
    raise RuntimeError("81195 is unexpectedly in verified PMML cohort.")

if "e1056" not in verified_keys:
    raise RuntimeError("e1056 is unexpectedly missing from PMML cohort.")

# ---------------------------------------------------------------
# 2. DATASET 2 -> VERIFIED SAMPLE ORDER
# ---------------------------------------------------------------

raw = pd.read_excel(
    XLSX,
    sheet_name="Supplemental Dataset 2",
    header=None,
    engine="openpyxl",
)

if raw.shape[1] != 180:
    raise RuntimeError(
        f"Expected 180 Dataset 2 columns, got {raw.shape[1]}."
    )

values = {}

feature_names = (
    raw.iloc[:, 0]
    .astype(str)
    .str.strip()
)

for feature in FEATURES:

    matches = np.where(
        feature_names.values == feature
    )[0]

    if len(matches) != 1:
        raise RuntimeError(
            f"Could not uniquely locate {feature}."
        )

    row = pd.to_numeric(
        raw.iloc[matches[0], 1:180],
        errors="coerce",
    ).to_numpy(float)

    if len(row) != 179:
        raise RuntimeError(
            f"{feature}: expected 179 values, got {len(row)}."
        )

    if not np.isfinite(row).all():
        raise RuntimeError(
            f"{feature} contains non-finite values."
        )

    values[feature] = row

mirna = pd.DataFrame(values)
mirna["sample_id"] = verified_ids
mirna["sample_key"] = verified_keys

# ---------------------------------------------------------------
# 3. TABLE 2 -> AGE + CA-125
# ---------------------------------------------------------------

tables = tables_from_docx(DOCX)

if len(tables) < 2:
    raise RuntimeError("Could not find Table 2.")

table2 = tables[1]
header = table2[0]

id_col = header.index("ID")
age_col = header.index("Age")
ca_col = header.index("CA-125")

clinical_rows = []

for row in table2[1:]:

    if len(row) <= max(id_col, age_col, ca_col):
        continue

    table_id = row[id_col].strip()

    if not table_id:
        continue

    key = norm(table_id)

    # Explicit e-sample.
    if key.startswith("e"):

        actual_key = key

    else:

        position = int(float(table_id))

        # Original numeric position 56 = excluded sample 81195.
        if position == 56:
            continue

        if not (1 <= position <= 120):
            raise RuntimeError(
                f"Unexpected Table 2 numeric ID: {position}"
            )

        # Dataset 2 omits original position 56.
        retained_index = (
            position - 1
            if position <= 55
            else position - 2
        )

        actual_key = verified_keys[retained_index]

    clinical_rows.append(
        {
            "sample_key": actual_key,
            "age": number(row[age_col]),
            "ca125": number(row[ca_col]),
        }
    )

clinical = pd.DataFrame(clinical_rows)

clinical = clinical.drop_duplicates(
    "sample_key",
    keep="first",
)

if len(clinical) != 179:
    raise RuntimeError(
        f"Expected 179 clinical rows, got {len(clinical)}."
    )

# ---------------------------------------------------------------
# 4. MERGE
# ---------------------------------------------------------------

df = mirna.merge(
    clinical,
    on="sample_key",
    how="left",
    validate="one_to_one",
)

df = df.merge(
    pmml[
        ["sample_id", "diagnosis", "true_label"]
    ].assign(
        sample_key=lambda x: x["sample_id"].map(norm)
    ),
    on="sample_key",
    how="left",
    validate="one_to_one",
)

# Use the verified PMML sample identity.
df["sample_id"] = df["sample_id_y"]
df["diagnosis"] = df["diagnosis"]
df["binary_label"] = df["true_label"]

# ---------------------------------------------------------------
# 5. FINAL VALIDATION
# ---------------------------------------------------------------

final_keys = set(
    df["sample_id"].map(norm)
)

verified_set = set(verified_keys)

print()
print("=" * 80)
print("FINAL VALIDATION")
print("=" * 80)

print("Patients:", len(df))
print("81195 present:", "81195" in final_keys)
print("e1056 present:", "e1056" in final_keys)
print("IDs exactly match verified cohort:", final_keys == verified_set)

print()
print("Age observed:", int(df["age"].notna().sum()))
print("Age missing:", int(df["age"].isna().sum()))

print("CA-125 observed:", int(df["ca125"].notna().sum()))
print("CA-125 missing:", int(df["ca125"].isna().sum()))

print()
print("e1056:")
print(
    df.loc[
        df["sample_id"].astype(str).str.lower() == "e1056",
        ["sample_id", "diagnosis", "age", "ca125"]
    ].to_string(index=False)
)

if final_keys != verified_set:
    raise RuntimeError(
        "Final IDs do not exactly match verified PMML cohort."
    )

if "81195" in final_keys:
    raise RuntimeError("81195 is still present.")

if "e1056" not in final_keys:
    raise RuntimeError("e1056 is still missing.")

if df["age"].isna().any():
    missing = df.loc[
        df["age"].isna(),
        "sample_id"
    ].tolist()

    raise RuntimeError(
        f"Missing age values for: {missing}"
    )

# ---------------------------------------------------------------
# 6. SAVE
# ---------------------------------------------------------------

keep = [
    "sample_id",
    "diagnosis",
    "binary_label",
    *FEATURES,
    "age",
    "ca125",
]

df[keep].to_csv(
    OUT,
    index=False,
    compression="gzip",
)

print()
print("=" * 80)
print("SUCCESS")
print("=" * 80)
print("Saved:", OUT)
