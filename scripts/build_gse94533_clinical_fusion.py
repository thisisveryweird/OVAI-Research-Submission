#!/usr/bin/env python3

"""Attach eLife Age and CA-125 values to the verified GSE94533 patients."""

import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
MATRIX = PROJECT / "data" / "processed" / "gse94533_mirna_tpm.csv.gz"
METADATA = PROJECT / "data" / "processed" / "gse94533_metadata.csv"
SUPPLEMENT = PROJECT / "data" / "external" / "elife_28932" / "elife_poa_e28932_Supplementary_File_6.xlsx"
TABLE_DOCX = PROJECT / "data" / "external" / "elife_28932" / "elife_poa_e28932_Supplementary_File_1.docx"
OUTPUT = PROJECT / "data" / "processed" / "gse94533_14mirna_clinical.csv.gz"

FEATURES = [
    "hsa-miR-23b-3p", "hsa-miR-29a-3p", "hsa-miR-32-5p",
    "hsa-miR-92a-3p", "hsa-miR-150-5p", "hsa-miR-200a-3p",
    "hsa-miR-200c-3p", "hsa-miR-203a", "hsa-miR-320c",
    "hsa-miR-320d", "hsa-miR-335-5p", "hsa-miR-450b-5p",
    "hsa-miR-1246", "hsa-miR-1307-5p",
]


def normalize(value):
    text = str(value).strip().lower()
    if re.fullmatch(r"\d+\.0+", text):
        text = text.split(".", 1)[0]
    return re.sub(r"[^a-z0-9]", "", text)


def cell_text(cell):
    return "".join(node.text or "" for node in cell.iter() if node.tag.endswith("}t")).strip()


def docx_tables(path):
    with zipfile.ZipFile(path) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
    tables = []
    for table in root.iter():
        if not table.tag.endswith("}tbl"):
            continue
        rows = []
        for tr in table:
            if tr.tag.endswith("}tr"):
                rows.append([cell_text(tc) for tc in tr if tc.tag.endswith("}tc")])
        if rows:
            tables.append(rows)
    return tables


def number(value):
    text = str(value).strip()
    if text in {"", "*", "‡", "nan"}:
        return np.nan
    try:
        return float(text)
    except ValueError:
        return np.nan


def subject_by_number():
    sheet = pd.read_excel(SUPPLEMENT, sheet_name="Supplemental Dataset 6", header=None)
    labels = sheet.iloc[:, 0].astype(str).str.strip()
    subject_row = sheet.index[labels.eq("Subject ID")]
    number_row = sheet.index[labels.eq("Sample Number")]
    if len(subject_row) != 1 or len(number_row) != 1:
        raise RuntimeError("Could not identify Subject ID/Sample Number rows.")
    result = {}
    for sample_number, subject in zip(
        sheet.iloc[number_row[0], 1:], sheet.iloc[subject_row[0], 1:]
    ):
        if pd.notna(sample_number) and pd.notna(subject):
            result[int(float(sample_number))] = normalize(subject)
    return result


def clinical_by_sample():
    tables = docx_tables(TABLE_DOCX)
    table = next(
        table for table in tables
        if table and {str(x).strip() for x in table[0]} >= {"ID", "Age", "CA-125"}
    )
    header = table[0]
    id_col, age_col, ca_col = header.index("ID"), header.index("Age"), header.index("CA-125")
    number_map = subject_by_number()
    result = {}
    for row in table[1:]:
        if len(row) <= max(id_col, age_col, ca_col) or not row[id_col].strip():
            continue
        raw_id = row[id_col].strip()
        key = normalize(raw_id)
        if raw_id.replace(".", "", 1).isdigit():
            key = number_map.get(int(float(raw_id)), "")
        if key:
            result[key] = {"age": number(row[age_col]), "ca125": number(row[ca_col])}
    return result


def main():
    metadata = pd.read_csv(METADATA)
    matrix = pd.read_csv(MATRIX)
    missing = [feature for feature in FEATURES if feature not in matrix]
    if missing:
        raise RuntimeError(f"Missing historical miRNAs: {missing}")
    clinical = clinical_by_sample()
    base = metadata[["sample_id_normalized", "diagnosis", "binary_label"]].merge(
        matrix[["sample_id_normalized"] + FEATURES],
        on="sample_id_normalized", how="inner", validate="one_to_one",
    )
    clinical_frame = pd.DataFrame.from_dict(clinical, orient="index")
    clinical_frame.index.name = "sample_id_normalized"
    clinical_frame = clinical_frame.reset_index()
    output = base.merge(clinical_frame, on="sample_id_normalized", how="left", validate="one_to_one")
    if len(output) != 179:
        raise RuntimeError(f"Expected 179 patients, found {len(output)}")
    # CA-125 is legitimately missing for a subset of the published cohort;
    # the modeling pipeline imputes it inside each training fold.
    if output["age"].notna().sum() < 170 or output["ca125"].notna().sum() < 100:
        raise RuntimeError("Too few matched clinical values; refusing an uncertain join.")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUTPUT, index=False, compression="gzip")
    print(f"Patients: {len(output)}")
    print(f"Age observed: {int(output['age'].notna().sum())}")
    print(f"CA-125 observed: {int(output['ca125'].notna().sum())}")
    print(f"Saved: {OUTPUT}")


if __name__ == "__main__":
    main()