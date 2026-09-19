#!/usr/bin/env python3

from __future__ import annotations

import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd


ROOT = Path.home() / "ovarian-ai"

SUPP = (
    ROOT
    / "data/external/elife_28932/"
    "elife_poa_e28932_Supplementary_File_6.xlsx"
)

DOCX = (
    ROOT
    / "data/external/elife_28932/"
    "elife_poa_e28932_Supplementary_File_1.docx"
)

PMML = (
    ROOT
    / "data/external/elife-28932-code1/"
    "Source-code1.xml"
)

GEO_METADATA = (
    ROOT
    / "data/interim/"
    "gse94533_geo_metadata_full.csv"
)


def norm_id(x):
    if pd.isna(x):
        return ""
    s = str(x).strip().lower()
    if re.fullmatch(r"\d+\.0+", s):
        s = s.split(".", 1)[0]
    return re.sub(r"[^a-z0-9]", "", s)


def get_sample_column(df):
    for col in df.columns:
        if str(col).strip().lower() in {"sample", "samples"}:
            return col

    for col in df.columns:
        name = str(col).strip().lower()
        if "sample" in name and "fastq" not in name:
            return col

    raise RuntimeError(
        f"Could not find sample column in columns: {list(df.columns)}"
    )


def read_sample_order(workbook):
    df = pd.read_excel(
        workbook,
        sheet_name="samples",
        header=1,
        engine="openpyxl",
    )

    col = get_sample_column(df)

    ids = [
        norm_id(x)
        for x in df[col].tolist()
        if norm_id(x)
    ]

    return ids


def recover_historical_ids():
    print("=" * 90)
    print("RECOVERING HISTORICAL SAMPLE ORDER")
    print("=" * 90)

    cohort2 = read_sample_order(
        ROOT
        / "data/raw/"
        "GSE94533_Processed_file_Cohort2.xlsx"
    )

    cohort1_all = read_sample_order(
        ROOT
        / "data/raw/"
        "GSE94533_Processed_file_Cohort1.xlsx"
    )

    print(f"Cohort 2 IDs: {len(cohort2)}")
    print(f"Cohort 1 IDs in samples sheet: {len(cohort1_all)}")

    if len(cohort2) != 120:
        raise RuntimeError(
            f"Expected 120 Cohort 2 samples, got {len(cohort2)}"
        )

    # ------------------------------------------------------------
    # Dataset 2 is the authoritative source for which ERASMOS
    # e-samples actually occur in the historical 179-sample matrix.
    #
    # Its header contains:
    #
    #   columns 1..120   -> positional numeric IDs
    #   columns 121..180 -> 60 actual e-sample IDs
    #
    # We read those 60 IDs directly.
    # ------------------------------------------------------------

    raw = pd.read_excel(
        SUPP,
        sheet_name="Supplemental Dataset 2",
        header=None,
        engine="openpyxl",
    )

    if raw.shape[1] != 180:
        raise RuntimeError(
            f"Expected Dataset 2 to have 180 columns, got {raw.shape[1]}"
        )

    header = raw.iloc[2, 1:].tolist()

    # Dataset 2 has 179 sample columns after the miRNA column.
    # The first 119 historical numeric patients represented after
    # exclusion are positions 1..55 and 57..120.
    #
    # The final 60 columns are the historical e-samples.
    e_sample_ids = [
        norm_id(x)
        for x in header[119:]
        if isinstance(x, str)
        and norm_id(x).startswith("e")
    ]

    if len(e_sample_ids) != 60:
        raise RuntimeError(
            "Expected exactly 60 e-samples in Dataset 2, got "
            f"{len(e_sample_ids)}:\n{e_sample_ids}"
        )

    # Verify every Dataset-2 e-sample exists in the Cohort 1
    # sample sheet, but do NOT require Cohort 1 itself to contain
    # only those 60 samples.
    cohort1_set = set(cohort1_all)

    missing_from_cohort1 = [
        x for x in e_sample_ids
        if x not in cohort1_set
    ]

    if missing_from_cohort1:
        raise RuntimeError(
            "Dataset 2 contains e-samples missing from the Cohort 1 "
            "sample sheet:\n"
            + "\n".join(missing_from_cohort1)
        )

    # The Dataset 2 header itself is authoritative for the e-sample
    # ordering. Preserve it exactly.
    cohort1_historical = e_sample_ids

    # ------------------------------------------------------------
    # Combine the two cohorts exactly as represented historically.
    # ------------------------------------------------------------

    combined = cohort2 + cohort1_historical

    if len(combined) != 180:
        raise RuntimeError(
            f"Expected 180 pre-exclusion samples, got {len(combined)}"
        )

    # Published Table 2 contains numeric patients:
    #
    #   1..55
    #   57..120
    #
    # Therefore numeric patient 56 is the omitted sample.
    excluded_position = 56
    excluded_id = combined[excluded_position - 1]

    print()
    print(
        f"Dataset 2 numeric position {excluded_position} -> "
        f"sample {excluded_id}"
    )

    if excluded_id != "81195":
        raise RuntimeError(
            "Expected published omitted patient 56 to map to 81195, "
            f"but recovered {excluded_id}"
        )

    historical = (
        combined[:55]
        + combined[56:]
    )

    if len(historical) != 179:
        raise RuntimeError(
            f"Expected 179 historical samples, got {len(historical)}"
        )

    print()
    print("Historical sample order:")
    print("  numeric positions 1-55")
    print("  numeric position 56 = 81195 EXCLUDED")
    print("  numeric positions 57-120")
    print("  all 60 Dataset 2 e-samples")

    print()
    print(f"Historical sample count: {len(historical)}")
    print(f"Historical e-sample count: {len(cohort1_historical)}")

    print()
    print("First 10:")
    for i, x in enumerate(historical[:10], 1):
        print(f"  {i:3d}: {x}")

    print()
    print("Around exclusion:")
    for i in range(53, 59):
        print(f"  {i:3d}: {historical[i - 1]}")

    print()
    print("Last 10:")
    start = len(historical) - 9
    for i, x in enumerate(historical[-10:], start):
        print(f"  {i:3d}: {x}")

    return historical



def load_dataset2(historical_ids):
    print()
    print("=" * 90)
    print("READING SUPPLEMENTAL DATASET 2")
    print("=" * 90)

    raw = pd.read_excel(
        SUPP,
        sheet_name="Supplemental Dataset 2",
        header=None,
        engine="openpyxl",
    )

    print(f"Raw shape: {raw.shape}")

    if raw.shape[1] != 180:
        raise RuntimeError(
            f"Expected 180 columns, got {raw.shape[1]}"
        )

    header = raw.iloc[2, 1:].tolist()

    print()
    print("Dataset 2 has 179 expression columns after the feature column?")

    # The workbook has 180 total columns:
    # column 0 = miRNA
    # columns 1..179 = 179 expression columns in the way pandas
    # presents this sheet.
    #
    # We therefore explicitly associate the 179 Dataset-2 columns
    # with the 179 recovered historical patients AFTER removing
    # numeric position 56 from the 180-sample biological order.
    #
    # Dataset 2 itself contains 179 columns here because the original
    # supplementary sheet already omits one numeric sample's values.

    sample_columns = list(range(1, 180))

    if len(sample_columns) != 179:
        raise RuntimeError("Dataset 2 expression column count is not 179.")

    # Verify e-sample tail.
    #
    # Dataset 2 contains 179 retained expression columns:
    #
    #   retained positions 1..119  = numeric cohort samples
    #   retained positions 120..179 = 60 e-samples
    #
    # The original biological numeric position 56 (sample 81195)
    # was excluded. Therefore the e-sample block begins at Python
    # index 119, not 120.
    #
    # `header` already contains the 179 Dataset-2 expression-column
    # headers, so there is no reason to reference undefined
    # `header_row` or `expression_columns`.

    if len(header) != 179:
        raise RuntimeError(
            "Dataset 2 header contains "
            f"{len(header)} expression columns instead of 179."
        )

    observed_tail = [
        norm_id(x)
        for x in header[119:]
    ]

    expected_tail = [
        norm_id(x)
        for x in historical_ids[119:]
    ]

    if observed_tail != expected_tail:
        raise RuntimeError(
            "Dataset 2 e-sample order does not match historical order.\\n"
            f"Observed: {observed_tail}\\n"
            f"Expected: {expected_tail}"
        )

    print()
    print("Dataset 2 e-sample order verified.")
    print(f"  First e-sample: {observed_tail[0]}")
    print(f"  Last e-sample:  {observed_tail[-1]}")

    # Reconstruct which biological sample each Dataset-2 column means.
    #
    # The 179 biological samples represented are:
    #   Cohort2 positions 1..55
    #   Cohort2 positions 57..120
    #   Cohort1 all 60
    #
    expected_ids = historical_ids

    features = [
        "hsa-miR-23b-3p",
        "hsa-miR-29a-3p",
        "hsa-miR-32-5p",
        "hsa-miR-92a-3p",
        "hsa-miR-150-5p",
        "hsa-miR-200a-3p",
        "hsa-miR-203a",
        "hsa-miR-320c",
        "hsa-miR-320d",
        "hsa-miR-335-5p",
        "hsa-miR-450b-5p",
        "hsa-miR-1246",
        "hsa-miR-1307-5p",
        "hsa-miR-200c-3p",
    ]

    # Preserve the exact PMML feature order separately later.
    available = {}

    for row_idx in range(3, len(raw)):
        feature = raw.iloc[row_idx, 0]

        if pd.isna(feature):
            continue

        feature = str(feature).strip()

        if feature in features:
            values = raw.iloc[row_idx, 1:180].tolist()

            if len(values) != 179:
                raise RuntimeError(
                    f"{feature}: expected 179 values, got {len(values)}"
                )

            available[feature] = values

    missing = [x for x in features if x not in available]

    if missing:
        raise RuntimeError(
            f"Missing Dataset 2 PMML features: {missing}"
        )

    matrix = pd.DataFrame(
        {
            feature: pd.to_numeric(
                available[feature],
                errors="coerce",
            )
            for feature in features
        },
        index=expected_ids,
    )

    print()
    print("Dataset 2 historical matrix:")
    print(f"  Patients: {matrix.shape[0]}")
    print(f"  Features: {matrix.shape[1]}")

    if matrix.isna().any().any():
        raise RuntimeError("Dataset 2 contains missing PMML feature values.")

    return matrix


def parse_pmml():
    print()
    print("=" * 90)
    print("READING ORIGINAL PMML")
    print("=" * 90)

    root = ET.parse(PMML).getroot()

    # PMML tags may have namespaces.
    def local(tag):
        return tag.split("}")[-1]

    norm = {}

    for elem in root.iter():
        if local(elem.tag) != "NormContinuous":
            continue

        field = elem.attrib.get("field")

        if not field or not field.startswith("hsa-miR-"):
            continue

        shift = float(elem.attrib["shift"])
        scale = float(elem.attrib["scale"])

        norm[field] = {
            "shift": shift,
            "scale": scale,
        }

    networks = [
        x for x in root.iter()
        if local(x.tag) == "NeuralNetwork"
    ]

    if len(networks) != 1:
        raise RuntimeError(
            f"Expected exactly one neural network, found {len(networks)}"
        )

    network = networks[0]

    layers = [
        x for x in network
        if local(x.tag) == "NeuralLayer"
    ]

    if len(layers) != 2:
        raise RuntimeError(
            f"Expected 2 neural layers, found {len(layers)}"
        )

    hidden = {}
    output = {}

    for layer_idx, layer in enumerate(layers):

        activation = layer.attrib.get(
            "activationFunction"
        )

        if activation != "tanh":
            raise RuntimeError(
                f"Expected tanh activation, got {activation}"
            )

        for neuron in layer:

            if local(neuron.tag) != "Neuron":
                continue

            nid = int(neuron.attrib["id"])
            bias = float(neuron.attrib["bias"])

            weights = {}

            for con in neuron:

                if local(con.tag) != "Con":
                    continue

                weights[int(con.attrib["from"])] = float(
                    con.attrib["weight"]
                )

            if layer_idx == 0:
                hidden[nid] = {
                    "bias": bias,
                    "weights": weights,
                }
            else:
                output[nid] = {
                    "bias": bias,
                    "weights": weights,
                }

    if len(hidden) != 7:
        raise RuntimeError(
            f"Expected 7 hidden neurons, found {len(hidden)}"
        )

    if len(output) != 2:
        raise RuntimeError(
            f"Expected 2 output neurons, found {len(output)}"
        )

    return norm, hidden, output


def forward_pass(matrix, norm, hidden, output, feature_order):
    x = matrix[feature_order].to_numpy(dtype=float)

    normalized = np.column_stack(
        [
            x[:, i] * norm[feature]["scale"]
            + norm[feature]["shift"]
            for i, feature in enumerate(feature_order)
        ]
    )

    hidden_ids = sorted(hidden)

    hidden_values = np.column_stack(
        [
            np.tanh(
                hidden[nid]["bias"]
                + normalized @ np.array(
                    [
                        hidden[nid]["weights"][i]
                        for i in range(len(feature_order))
                    ]
                )
            )
            for nid in hidden_ids
        ]
    )

    output_ids = sorted(output)

    output_values = np.column_stack(
        [
            np.tanh(
                output[nid]["bias"]
                + hidden_values @ np.array(
                    [
                        output[nid]["weights"][hid]
                        for hid in hidden_ids
                    ]
                )
            )
            for nid in output_ids
        ]
    )

    # neuron 21 = cancer
    # neuron 22 = controls/borderline
    predicted = np.where(
        output_values[:, 0] > output_values[:, 1],
        "CA",
        "BBC",
    )

    return output_values, predicted


def parse_table2():
    print()
    print("=" * 90)
    print("READING PUBLISHED TABLE 2")
    print("=" * 90)

    with zipfile.ZipFile(DOCX) as z:
        xml = z.read("word/document.xml")

    root = ET.fromstring(xml)

    ns = {
        "w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
    }

    tables = root.findall(".//w:tbl", ns)

    print(f"Tables found: {len(tables)}")

    if len(tables) < 2:
        raise RuntimeError("Could not find Table 2.")

    table = tables[1]

    rows = []

    for tr in table.findall("./w:tr", ns):

        cells = []

        for tc in tr.findall("./w:tc", ns):

            texts = tc.findall(".//w:t", ns)

            value = "".join(
                t.text or ""
                for t in texts
            ).strip()

            cells.append(value)

        rows.append(cells)

    header = rows[0]

    print()
    print("Table 2 header:")
    for i, h in enumerate(header):
        print(f"  {i}: {h}")

    required = {
        "ID",
        "Study",
        "NN",
    }

    if not required.issubset(set(header)):
        raise RuntimeError(
            f"Table 2 header missing required fields: {required}"
        )

    id_idx = header.index("ID")
    nn_idx = header.index("NN")

    records = []

    for row in rows[1:]:

        if len(row) <= max(id_idx, nn_idx):
            continue

        sid = norm_id(row[id_idx])
        nn = row[nn_idx].strip()

        if not sid:
            continue

        if nn not in {"CA", "BBC"}:
            continue

        records.append(
            {
                "sample_id": sid,
                "published_nn": nn,
            }
        )

    table2 = pd.DataFrame(records)

    print()
    print(f"Published Table 2 prediction rows: {len(table2)}")

    if len(table2) != 179:
        raise RuntimeError(
            f"Expected 179 Table 2 prediction rows, got {len(table2)}"
        )

    if table2["sample_id"].duplicated().any():
        raise RuntimeError("Duplicate IDs found in Table 2.")

    return table2


def main():
    historical_ids = recover_historical_ids()

    dataset2 = load_dataset2(historical_ids)

    norm, hidden, output = parse_pmml()

    feature_order = [
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

    missing_norm = [
        f for f in feature_order
        if f not in norm
    ]

    if missing_norm:
        raise RuntimeError(
            f"Missing PMML normalization fields: {missing_norm}"
        )

    output_values, predicted = forward_pass(
        dataset2,
        norm,
        hidden,
        output,
        feature_order,
    )

    predictions = pd.DataFrame(
        {
            "sample_id": dataset2.index,
            "pmml_predicted": predicted,
            "cancer_activation": output_values[:, 0],
            "bbc_activation": output_values[:, 1],
        }
    )

    table2 = parse_table2()

    merged = predictions.merge(
        table2,
        on="sample_id",
        how="outer",
        indicator=True,
    )

    if not (merged["_merge"] == "both").all():
        print()
        print("ID mismatch:")
        print(
            merged[
                merged["_merge"] != "both"
            ].to_string(index=False)
        )
        raise RuntimeError(
            "Historical sample IDs do not exactly match Table 2."
        )

    merged["match"] = (
        merged["pmml_predicted"]
        == merged["published_nn"]
    )

    agreement = merged["match"].mean()

    print()
    print("=" * 90)
    print("HISTORICAL PMML vs PUBLISHED TABLE 2")
    print("=" * 90)

    print()
    print(f"Patients compared: {len(merged)}")
    print(f"Exact prediction agreement: {agreement:.6f}")
    print(
        f"Matching predictions: "
        f"{merged['match'].sum()} / {len(merged)}"
    )

    mismatches = merged[
        ~merged["match"]
    ]

    print()
    print(f"Mismatches: {len(mismatches)}")

    if len(mismatches):
        print()
        print("FIRST MISMATCHES:")
        print(
            mismatches[
                [
                    "sample_id",
                    "published_nn",
                    "pmml_predicted",
                    "cancer_activation",
                    "bbc_activation",
                ]
            ]
            .head(30)
            .to_string(index=False)
        )

    out = (
        ROOT
        / "reports/tables/"
        "historical_pmml_vs_published_table2.csv"
    )

    merged.to_csv(out, index=False)

    print()
    print("Saved:")
    print(f"  {out}")

    if agreement == 1.0:
        print()
        print("=" * 90)
        print("SUCCESS: EXACT HISTORICAL REPRODUCTION")
        print("=" * 90)
        print(
            "The PMML forward pass reproduces every published "
            "Table 2 NN classification."
        )
    else:
        print()
        print("=" * 90)
        print("REPRODUCTION NOT YET EXACT")
        print("=" * 90)
        print(
            "The historical mapping is structurally correct, "
            "but the PMML/data representation still needs investigation."
        )


if __name__ == "__main__":
    main()
