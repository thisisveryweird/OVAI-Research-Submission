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
    / "pmml_vs_published_table2_comparison.csv"
)


# ================================================================
# ID NORMALIZATION
# ================================================================

def norm_id(value) -> str:
    s = str(value).strip().lower()

    if re.fullmatch(r"\d+\.0+", s):
        s = s.split(".", 1)[0]

    return re.sub(r"[^a-z0-9]", "", s)


# ================================================================
# DOCX TABLE EXTRACTION
# ================================================================

def extract_docx_tables(path: Path):

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

            cells = []

            for tc in tr:

                if not tc.tag.endswith("}tc"):
                    continue

                texts = []

                for node in tc.iter():

                    if node.tag.endswith("}t"):
                        texts.append(node.text or "")

                cells.append(
                    "".join(texts).strip()
                )

            if cells:
                rows.append(cells)

        if rows:
            tables.append(rows)

    return tables


# ================================================================
# FIND TABLE 2
# ================================================================

def find_table2(tables):

    for table_number, table in enumerate(tables, 1):

        for row in table:

            normalized = {
                str(x).strip().lower()
                for x in row
            }

            if {
                "id",
                "study",
                "stage",
                "histology",
                "nn",
            }.issubset(normalized):

                print()
                print(
                    f"Identified Table 2 as extracted table "
                    f"{table_number}."
                )

                return table

    raise RuntimeError(
        "Could not identify Table 2 from its header contents."
    )


# ================================================================
# FIND HEADER ROW
# ================================================================

def find_header_row(table):

    required = {
        "id",
        "study",
        "stage",
        "grade",
        "histology",
        "age",
        "ca-125",
        "nn",
    }

    best = None

    for i, row in enumerate(table):

        normalized = {
            str(x).strip().lower()
            for x in row
        }

        overlap = len(required & normalized)

        if overlap >= 6:
            best = i
            break

    if best is None:
        print()
        print("Could not find Table 2 header.")
        print()
        print("First rows extracted from candidate table:")

        for i, row in enumerate(table[:10]):
            print(f"{i}: {row}")

        raise RuntimeError(
            "Could not locate Table 2 header row."
        )

    return best


# ================================================================
# PUBLISHED LABEL NORMALIZATION
# ================================================================

def published_label(value) -> int:

    s = str(value).strip().upper()

    if s == "CA":
        return 1

    if s == "BBC":
        return 0

    raise ValueError(
        f"Unexpected published NN prediction: {value!r}"
    )


# ================================================================
# MAIN
# ================================================================

def main():

    print("=" * 90)
    print("PMML vs PUBLISHED eLIFE TABLE 2 VALIDATION")
    print("=" * 90)

    # ------------------------------------------------------------
    # Load PMML predictions
    # ------------------------------------------------------------

    pred = pd.read_csv(PREDICTIONS)

    print()
    print("PMML prediction file:")
    print(f"  Patients: {len(pred)}")

    if len(pred) != 179:
        raise RuntimeError(
            f"Expected 179 PMML predictions, got {len(pred)}."
        )

    pred["sample_id"] = pred["sample_id"].map(norm_id)

    if pred["sample_id"].duplicated().any():
        raise RuntimeError(
            "PMML prediction file contains duplicate sample IDs."
        )

    # ------------------------------------------------------------
    # Historical order from our reconstructed PMML predictions.
    # ------------------------------------------------------------

    historical_ids = pred["sample_id"].tolist()

    # ------------------------------------------------------------
    # Read DOCX
    # ------------------------------------------------------------

    print()
    print("Reading Supplementary File 1...")

    tables = extract_docx_tables(DOCX)

    print(f"Tables extracted: {len(tables)}")

    table2 = find_table2(tables)

    header_idx = find_header_row(table2)

    headers = [
        str(x).strip()
        for x in table2[header_idx]
    ]

    print()
    print("Table 2 header:")

    for i, h in enumerate(headers):
        print(f"  {i}: {h}")

    # ------------------------------------------------------------
    # Parse rows following header.
    # ------------------------------------------------------------

    records = []

    for row in table2[header_idx + 1:]:

        if not row:
            continue

        identifier = str(row[0]).strip()

        # Valid publication IDs:
        #   1
        #   2
        #   ...
        #   120
        #   e1001
        #
        if not (
            re.fullmatch(r"\d+", identifier)
            or re.fullmatch(r"e\d+", identifier.lower())
        ):
            continue

        if len(row) < len(headers):
            row = row + [""] * (
                len(headers) - len(row)
            )

        records.append(row[:len(headers)])

    print()
    print(
        f"Published Table 2 patient rows: {len(records)}"
    )

    if len(records) != 179:
        raise RuntimeError(
            "Table 2 should contain exactly 179 patient rows."
        )

    published = pd.DataFrame(
        records,
        columns=headers,
    )

    # ------------------------------------------------------------
    # Locate NN column robustly.
    # ------------------------------------------------------------

    nn_col = None

    for c in published.columns:

        if str(c).strip().lower() == "nn":
            nn_col = c
            break

    if nn_col is None:
        raise RuntimeError(
            "Could not find NN column in Table 2."
        )

    # ------------------------------------------------------------
    # CRITICAL MAPPING
    #
    # Table 2 uses publication IDs.
    #
    # Our PMML predictions use real GSE sample IDs.
    #
    # Their correspondence is established by identical historical
    # patient order.
    #
    # Therefore:
    #
    # Table 2 row 1
    #    -> historical position 1
    #    -> GSE sample 50217
    #
    # Table 2 row 2
    #    -> historical position 2
    #    -> GSE sample 80675
    #
    # etc.
    # ------------------------------------------------------------

    published["historical_position"] = range(
        1,
        180,
    )

    published["sample_id"] = historical_ids

    published["published_nn"] = (
        published[nn_col]
        .map(published_label)
    )

    # ------------------------------------------------------------
    # PMML predicted labels
    # ------------------------------------------------------------

    pred["pmml_label"] = (
        pred["predicted_class"]
        .astype(str)
        .str.contains(
            "Ovarian Cancer",
            case=False,
            regex=False,
        )
        .astype(int)
    )

    # ------------------------------------------------------------
    # Compare
    # ------------------------------------------------------------

    comparison = published[
        [
            "historical_position",
            "ID",
            "sample_id",
            nn_col,
            "published_nn",
        ]
    ].merge(
        pred[
            [
                "sample_id",
                "cancer_activation",
                "controls_or_borderline_activation",
                "predicted_class",
                "pmml_label",
            ]
        ],
        on="sample_id",
        how="outer",
        validate="one_to_one",
        indicator=True,
    )

    if not (comparison["_merge"] == "both").all():

        print()
        print("FAILED TO MATCH ALL 179 PATIENTS")

        print(
            comparison[
                comparison["_merge"] != "both"
            ].to_string(index=False)
        )

        raise RuntimeError(
            "Patient mapping is incomplete."
        )

    comparison["exact_match"] = (
        comparison["published_nn"]
        == comparison["pmml_label"]
    )

    matches = int(
        comparison["exact_match"].sum()
    )

    mismatches = (
        len(comparison) - matches
    )

    agreement = accuracy_score(
        comparison["published_nn"],
        comparison["pmml_label"],
    )

    cm = confusion_matrix(
        comparison["published_nn"],
        comparison["pmml_label"],
        labels=[0, 1],
    )

    # ------------------------------------------------------------
    # Results
    # ------------------------------------------------------------

    print()
    print("=" * 90)
    print("RESULTS")
    print("=" * 90)

    print()
    print(f"Patients:                 {len(comparison)}")
    print(f"Exact matches:            {matches}")
    print(f"Mismatches:               {mismatches}")
    print(f"Prediction agreement:     {agreement:.6f}")

    print()
    print("Confusion matrix:")
    print()
    print("                    PMML")
    print("                 BBC     CA")
    print(
        f"Published BBC  {cm[0,0]:6d} {cm[0,1]:6d}"
    )
    print(
        f"Published CA   {cm[1,0]:6d} {cm[1,1]:6d}"
    )

    # ------------------------------------------------------------
    # Detailed mismatches
    # ------------------------------------------------------------

    if mismatches > 0:

        print()
        print("=" * 90)
        print("MISMATCHED PATIENTS")
        print("=" * 90)

        mismatched = comparison[
            ~comparison["exact_match"]
        ]

        print(
            mismatched[
                [
                    "historical_position",
                    "ID",
                    "sample_id",
                    nn_col,
                    "predicted_class",
                    "cancer_activation",
                    "controls_or_borderline_activation",
                ]
            ].to_string(index=False)
        )

    else:

        print()
        print("=" * 90)
        print("PERFECT AGREEMENT")
        print("=" * 90)

        print(
            "All 179 published neural-network predictions "
            "match the reconstructed PMML model."
        )

    # ------------------------------------------------------------
    # Save
    # ------------------------------------------------------------

    comparison.to_csv(
        OUTPUT,
        index=False,
    )

    print()
    print("Saved comparison:")
    print(f"  {OUTPUT}")


if __name__ == "__main__":
    main()
