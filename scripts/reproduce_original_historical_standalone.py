#!/usr/bin/env python3

"""
STANDALONE HISTORICAL REPRODUCTION OF THE ORIGINAL eLIFE
GSE94533 14-miRNA NEURAL NETWORK

This script intentionally does NOT use the previous
reproduce_original_pmml_historical_data.py script.

Reason:
    that script accumulated several incorrect assumptions about
    Supplemental Dataset 2 sample identifiers.

This standalone implementation directly establishes:

    Supplemental Dataset 2
            |
            v
    historical 179 sample columns
            |
            v
    first 120 positional columns -> Cohort 2 IDs
    final 59 explicit columns    -> explicit Cohort 1 IDs
            |
            v
    exact 14 PMML miRNAs
            |
            v
    PMML normalization
            |
            v
    14 -> 7 tanh -> 2 tanh
            |
            v
    cancer vs non-cancer prediction

NO MODEL TRAINING IS PERFORMED.

IMPORTANT:
------------
We deliberately preserve the historical supplementary dataset as it
is represented in Supplemental Dataset 2.

We do NOT remove E-1056 merely because an earlier cleaned dataset
removed it.

The historical representation must be reproduced independently from
our later cleaned modern dataset.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd


# =====================================================================
# PATHS
# =====================================================================

ROOT = Path.home() / "ovarian-ai"

SUPPLEMENTARY_XLSX = (
    ROOT
    / "data"
    / "external"
    / "elife_28932"
    / "elife_poa_e28932_Supplementary_File_6.xlsx"
)

PMML_PATH = (
    ROOT
    / "data"
    / "external"
    / "elife-28932-code1"
    / "Source-code1.xml"
)

GEO_METADATA = (
    ROOT
    / "data"
    / "interim"
    / "gse94533_geo_metadata_full.csv"
)

COHORT1_XLSX = (
    ROOT
    / "data"
    / "raw"
    / "GSE94533_Processed_file_Cohort1.xlsx"
)

COHORT2_XLSX = (
    ROOT
    / "data"
    / "raw"
    / "GSE94533_Processed_file_Cohort2.xlsx"
)

OUTPUT_DIR = ROOT / "reports" / "tables"

HISTORICAL_MATRIX_OUTPUT = (
    OUTPUT_DIR
    / "historical_dataset2_14mirna_matrix.csv"
)

PREDICTIONS_OUTPUT = (
    OUTPUT_DIR
    / "historical_original_pmml_predictions.csv"
)

AUDIT_OUTPUT = (
    OUTPUT_DIR
    / "historical_original_pmml_standalone_audit.json"
)


# =====================================================================
# EXACT PMML FEATURES
# =====================================================================

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


# =====================================================================
# SAMPLE ID NORMALIZATION
# =====================================================================

def normalize_id(value: object) -> str:

    if value is None:
        return ""

    if pd.isna(value):
        return ""

    s = str(value).strip().lower()

    # Excel numeric representation:
    # 50217.0 -> 50217
    if re.fullmatch(r"\d+\.0+", s):
        s = s.split(".", 1)[0]

    # E-1001 -> e1001
    s = re.sub(r"[^a-z0-9]", "", s)

    return s


# =====================================================================
# READ ORIGINAL PROCESSED WORKBOOK SAMPLE ORDER
# =====================================================================

def get_processed_sample_order():

    print()
    print("=" * 100)
    print("RECOVERING ORIGINAL PROCESSED-WORKBOOK SAMPLE ORDER")
    print("=" * 100)

    def read_order(path):

        header = pd.read_excel(
            path,
            sheet_name="mirna_tpm",
            header=2,
            nrows=0,
            engine="openpyxl",
        )

        result = []

        for column in header.columns:

            normalized = normalize_id(column)

            if not normalized:
                continue

            # Exclude metadata columns.
            if normalized in {
                "name",
                "sequence",
                "namecount",
                "identicalmirnas",
            }:
                continue

            # GSE94533 sample identifiers are either:
            #   e#### / E-####
            # or numeric identifiers.
            if normalized.startswith("e"):

                if normalized[1:].isdigit():
                    result.append(normalized)

            elif normalized.isdigit():

                result.append(normalized)

        return result

    cohort1 = read_order(COHORT1_XLSX)
    cohort2 = read_order(COHORT2_XLSX)

    print()
    print(f"Cohort 1 samples: {len(cohort1)}")
    print(f"Cohort 2 samples: {len(cohort2)}")

    print()
    print("Cohort 2 first 10:")
    for x in cohort2[:10]:
        print(f"  {x}")

    print()
    print("Cohort 1 first 10:")
    for x in cohort1[:10]:
        print(f"  {x}")

    if len(cohort2) != 120:
        raise RuntimeError(
            "Expected exactly 120 Cohort 2 samples, "
            f"found {len(cohort2)}."
        )

    if len(cohort1) != 60:
        raise RuntimeError(
            "Expected exactly 60 Cohort 1 samples, "
            f"found {len(cohort1)}."
        )

    return cohort2, cohort1


# =====================================================================
# LOAD SUPPLEMENTAL DATASET 2
# =====================================================================

def load_dataset2():

    print()
    print("=" * 100)
    print("READING SUPPLEMENTAL DATASET 2")
    print("=" * 100)

    raw = pd.read_excel(
        SUPPLEMENTARY_XLSX,
        sheet_name="Supplemental Dataset 2",
        header=None,
        engine="openpyxl",
    )

    print()
    print(f"Raw shape: {raw.shape}")

    # ---------------------------------------------------------------
    # Find the header row containing:
    #
    #   miRNA | 1 | 2 | 3 | ... | 120 | e1036 | e1048 | ...
    # ---------------------------------------------------------------

    header_row = None

    for r in range(min(10, len(raw))):

        values = [
            normalize_id(x)
            for x in raw.iloc[r].tolist()
        ]

        # Header row should contain many non-empty sample columns.
        numeric_count = 0
        total_identifier_like = 0

        for value in values:

            if not value:
                continue

            if value.isdigit():
                numeric_count += 1
                total_identifier_like += 1

            elif value.startswith("e") and value[1:].isdigit():
                total_identifier_like += 1

        if (
            numeric_count >= 100
            and total_identifier_like >= 150
        ):
            header_row = r
            break

    if header_row is None:

        raise RuntimeError(
            "Could not identify Supplemental Dataset 2 "
            "sample-header row."
        )

    print()
    print(
        f"Detected Dataset 2 header row: {header_row}"
    )

    header_values = raw.iloc[header_row].tolist()

    # First column is miRNA.
    feature_column = 0

    # All remaining sample columns should total 179.
    sample_columns = list(range(1, len(header_values)))

    sample_columns = [
        c
        for c in sample_columns
        if not pd.isna(header_values[c])
    ]

    print(
        f"Detected sample columns: {len(sample_columns)}"
    )

    if len(sample_columns) != 179:

        raise RuntimeError(
            "Expected exactly 179 Dataset 2 sample columns; "
            f"found {len(sample_columns)}."
        )

    raw_header_ids = [
        header_values[c]
        for c in sample_columns
    ]

    normalized_header_ids = [
        normalize_id(x)
        for x in raw_header_ids
    ]

    print()
    print("FIRST 20 RAW DATASET 2 HEADER VALUES:")

    for i, value in enumerate(
        raw_header_ids[:20],
        start=1,
    ):
        print(
            f"  column {i:3d}: {value!r}"
        )

    print()
    print("LAST 20 RAW DATASET 2 HEADER VALUES:")

    start = max(
        0,
        len(raw_header_ids) - 20,
    )

    for i in range(start, len(raw_header_ids)):

        print(
            f"  column {i + 1:3d}: "
            f"{raw_header_ids[i]!r}"
        )

    # ---------------------------------------------------------------
    # Recover the actual sample IDs.
    #
    # Dataset 2 uses:
    #
    #   first 120 columns:
    #       1, 2, ..., 120
    #
    #   final 59 columns:
    #       explicit Cohort 1 identifiers
    #
    # The first 120 are POSITION NUMBERS.
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    # ROBUST CHECK FOR DATASET 2 POSITIONAL HEADERS
    #
    # Excel/pandas may represent:
    #
    #     1
    #     2.0
    #     3.0
    #     ...
    #
    # as different Python/numpy types.
    #
    # Therefore we compare numerically rather than comparing strings.
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    # Robust Dataset 2 positional-header detection.
    #
    # The first 120 Dataset 2 columns may be represented by pandas
    # as Python ints, floats, or NumPy floating-point values.
    #
    # Do NOT use sample_ids here because that list is constructed
    # later in this function.
    # ---------------------------------------------------------------
    
    first_120_are_positions = False
    
    try:
        first_120_raw = [
            raw.iloc[header_row, column]
            for column in range(1, 121)
        ]
    
        first_120_are_positions = (
            len(first_120_raw) == 120
            and all(
                float(value) == float(position)
                for position, value
                in enumerate(
                    first_120_raw,
                    start=1,
                )
            )
        )
    
    except (TypeError, ValueError, IndexError):
        first_120_are_positions = False

    for i in range(120):

        value = raw_header_ids[i]

        try:
            numeric_value = float(value)

        except (TypeError, ValueError):
            first_120_are_positions = False
            break

        if not math.isfinite(numeric_value):
            first_120_are_positions = False
            break

        expected_value = float(i + 1)

        if not math.isclose(
            numeric_value,
            expected_value,
            rel_tol=0.0,
            abs_tol=1e-9,
        ):
            first_120_are_positions = False
            break

    if not first_120_are_positions:

        raise RuntimeError(
            "Dataset 2 does not have the expected first 120 "
            "positional headers 1..120."
        )

    print()
    print(
        "CONFIRMED:"
    )
    print(
        "  Dataset 2 columns 1-120 are positional sample numbers."
    )
    print(
        "  Dataset 2 columns 121-179 contain explicit IDs."
    )

    cohort2_order, cohort1_order = (
        get_processed_sample_order()
    )

    historical_sample_ids = (
        cohort2_order
        + normalized_header_ids[120:]
    )

    print()
    print(
        f"Recovered historical sample count: "
        f"{len(historical_sample_ids)}"
    )

    if len(historical_sample_ids) != 179:

        raise RuntimeError(
            "Historical sample list does not contain 179 samples."
        )

    if len(set(historical_sample_ids)) != 179:

        duplicates = []

        seen = set()

        for value in historical_sample_ids:

            if value in seen:
                duplicates.append(value)

            seen.add(value)

        raise RuntimeError(
            "Historical sample IDs are not unique.\n"
            + "\n".join(
                f"  {x}"
                for x in sorted(set(duplicates))
            )
        )

    # ---------------------------------------------------------------
    # IMPORTANT:
    #
    # We do NOT impose the previous E-1056 exclusion here.
    #
    # Dataset 2 is the historical representation being reproduced.
    # ---------------------------------------------------------------

    print()
    print(
        "Historical E-1056 present:"
        f" {'e1056' in historical_sample_ids}"
    )

    if "e1056" in historical_sample_ids:

        print()
        print(
            "NOTE: E-1056 is present in Supplemental Dataset 2."
        )

        print(
            "Therefore it is being retained for this historical"
        )

        print(
            "reproduction. We will not apply the later cleaned-dataset"
        )

        print(
            "exclusion rule here."
        )

    print()
    print(
        "FIRST 10 HISTORICAL SAMPLE MAPPINGS:"
    )

    for position in range(10):

        print(
            f"  Dataset2 position {position + 1:3d}"
            f" -> {historical_sample_ids[position]}"
        )

    print()
    print(
        "LAST 10 HISTORICAL SAMPLE MAPPINGS:"
    )

    for position in range(169, 179):

        print(
            f"  Dataset2 position {position + 1:3d}"
            f" -> {historical_sample_ids[position]}"
        )

    # ---------------------------------------------------------------
    # Locate exact feature rows.
    # ---------------------------------------------------------------

    feature_rows = {}

    for feature in FEATURES:

        matches = []

        for r in range(
            header_row + 1,
            len(raw),
        ):

            value = raw.iloc[
                r,
                feature_column,
            ]

            if (
                not pd.isna(value)
                and str(value).strip() == feature
            ):
                matches.append(r)

        if len(matches) != 1:

            raise RuntimeError(
                f"Feature {feature!r} was found "
                f"{len(matches)} times."
            )

        feature_rows[feature] = matches[0]

    print()
    print(
        f"Exact PMML features found: {len(feature_rows)}"
    )

    # ---------------------------------------------------------------
    # Construct matrix.
    #
    # Rows = patients
    # Columns = exact 14 miRNAs.
    # ---------------------------------------------------------------

    matrix = np.zeros(
        (179, 14),
        dtype=np.float64,
    )

    for j, feature in enumerate(FEATURES):

        row = feature_rows[feature]

        values = []

        for c in sample_columns:

            value = raw.iloc[row, c]

            if pd.isna(value):

                raise RuntimeError(
                    f"Missing Dataset 2 value for "
                    f"{feature}, column {c}."
                )

            values.append(float(value))

        matrix[:, j] = values

    # ---------------------------------------------------------------
    # Create patient-level dataframe.
    # ---------------------------------------------------------------

    df = pd.DataFrame(
        matrix,
        columns=FEATURES,
    )

    df.insert(
        0,
        "sample_id",
        historical_sample_ids,
    )

    print()
    print(
        "Historical representation:"
    )
    print(
        "  batch adjusted"
    )
    print(
        "  log10 transformed"
    )
    print(
        f"  patients: {len(df)}"
    )
    print(
        f"  PMML features: {len(FEATURES)}"
    )

    if not np.isfinite(matrix).all():

        raise RuntimeError(
            "Historical matrix contains non-finite values."
        )

    return df


# =====================================================================
# LOAD GEO DIAGNOSIS
# =====================================================================

def load_diagnosis():

    print()
    print("=" * 100)
    print("LOADING GEO DIAGNOSIS METADATA")
    print("=" * 100)

    metadata = pd.read_csv(
        GEO_METADATA
    )

    metadata["normalized_id"] = (
        metadata["title"].map(normalize_id)
    )

    metadata = metadata[
        metadata["normalized_id"] != ""
    ].copy()

    lookup = {}

    for _, row in metadata.iterrows():

        sid = row["normalized_id"]

        if sid in lookup:

            raise RuntimeError(
                "Duplicate GEO sample identifier: "
                f"{sid}"
            )

        lookup[sid] = row

    print()
    print(
        f"GEO diagnosis records: {len(lookup)}"
    )

    return lookup


# =====================================================================
# LOCAL NAME HELPER
# =====================================================================

def local_name(element):

    if "}" in element.tag:

        return element.tag.rsplit(
            "}",
            1,
        )[1]

    return element.tag


# =====================================================================
# PARSE PMML
# =====================================================================

def parse_pmml():

    print()
    print("=" * 100)
    print("PARSING ORIGINAL PMML")
    print("=" * 100)

    root = ET.parse(
        PMML_PATH
    ).getroot()

    neural_networks = []

    for element in root.iter():

        if local_name(element) == "NeuralNetwork":

            neural_networks.append(element)

    if len(neural_networks) != 1:

        raise RuntimeError(
            f"Expected exactly one NeuralNetwork; "
            f"found {len(neural_networks)}."
        )

    network = neural_networks[0]

    model_name = network.attrib.get(
        "modelName"
    )

    function_name = network.attrib.get(
        "functionName"
    )

    print()
    print(
        f"modelName: {model_name}"
    )

    print(
        f"functionName: {function_name}"
    )

    # ---------------------------------------------------------------
    # Inputs
    # ---------------------------------------------------------------

    pmml_inputs = []

    for element in network:

        if local_name(element) != "NeuralInputs":
            continue

        for neural_input in element:

            if local_name(neural_input) != "NeuralInput":
                continue

            norm_continuous = None

            for child in neural_input.iter():

                if local_name(child) == "NormContinuous":

                    norm_continuous = child
                    break

            if norm_continuous is None:

                continue

            field = norm_continuous.attrib[
                "field"
            ]

            shift = float(
                norm_continuous.attrib[
                    "shift"
                ]
            )

            scale = float(
                norm_continuous.attrib[
                    "scale"
                ]
            )

            linear_norms = []

            for child in norm_continuous:

                if local_name(child) == "LinearNorm":

                    linear_norms.append(
                        (
                            float(
                                child.attrib["orig"]
                            ),
                            float(
                                child.attrib["norm"]
                            ),
                        )
                    )

            pmml_inputs.append(
                {
                    "field": field,
                    "shift": shift,
                    "scale": scale,
                    "linear_norms": linear_norms,
                }
            )

    if len(pmml_inputs) != 14:

        raise RuntimeError(
            "Expected 14 PMML neural inputs; "
            f"found {len(pmml_inputs)}."
        )

    pmml_feature_names = [
        x["field"]
        for x in pmml_inputs
    ]

    print()
    print(
        "PMML input features:"
    )

    for i, feature in enumerate(
        pmml_feature_names,
        start=1,
    ):

        print(
            f"  {i:2d}. {feature}"
        )

    if pmml_feature_names != FEATURES:

        raise RuntimeError(
            "PMML feature order does not exactly match "
            "the recovered 14-feature order.\n"
            f"PMML: {pmml_feature_names}\n"
            f"Expected: {FEATURES}"
        )

    # ---------------------------------------------------------------
    # Neural layers
    # ---------------------------------------------------------------

    layers = []

    for element in network:

        if local_name(element) != "NeuralLayer":
            continue

        activation = element.attrib.get(
            "activation"
        )

        neurons = []

        for neuron in element:

            if local_name(neuron) != "Neuron":
                continue

            neuron_id = neuron.attrib.get(
                "id"
            )

            bias = float(
                neuron.attrib.get(
                    "bias",
                    "0",
                )
            )

            connections = {}

            for child in neuron:

                if local_name(child) != "Con":
                    continue

                source = child.attrib.get(
                    "from"
                )

                weight = float(
                    child.attrib.get(
                        "weight",
                        "0",
                    )
                )

                connections[source] = weight

            neurons.append(
                {
                    "id": neuron_id,
                    "bias": bias,
                    "connections": connections,
                }
            )

        layers.append(
            {
                "activation": activation,
                "neurons": neurons,
            }
        )

    if len(layers) != 2:

        raise RuntimeError(
            f"Expected 2 neural layers; found {len(layers)}."
        )

    print()
    print(
        "NETWORK ARCHITECTURE:"
    )

    for i, layer in enumerate(
        layers,
        start=1,
    ):

        print(
            f"  Layer {i}: "
            f"{len(layer['neurons'])} neurons "
            f"activation={layer['activation']}"
        )

    if (
        len(layers[0]["neurons"]) != 7
        or len(layers[1]["neurons"]) != 2
    ):

        raise RuntimeError(
            "Unexpected neural-network architecture."
        )

    # ---------------------------------------------------------------
    # Confirm output neuron IDs.
    # ---------------------------------------------------------------

    output_ids = [
        neuron["id"]
        for neuron in layers[1]["neurons"]
    ]

    print()
    print(
        "OUTPUT NEURONS:"
    )

    for neuron_id in output_ids:

        print(
            f"  {neuron_id}"
        )

    # From the audited PMML:
    #
    #   21 = Ovarian Cancer
    #   22 = Controls or Borderline

    if set(output_ids) != {"21", "22"}:

        raise RuntimeError(
            "Expected output neurons 21 and 22; "
            f"found {output_ids}."
        )

    pmml = {
        "model_name": model_name,
        "function_name": function_name,
        "inputs": pmml_inputs,
        "layers": layers,
        "cancer_output_neuron": "21",
        "noncancer_output_neuron": "22",
    }

    return pmml


# =====================================================================
# PMML FORWARD PASS
# =====================================================================

def forward_pass(
    X,
    pmml,
):

    # ---------------------------------------------------------------
    # Input normalization:
    #
    # The supplied PMML audit demonstrated that these particular
    # NormContinuous transforms behave as:
    #
    #     normalized = x * scale + shift
    #
    # and the LinearNorm endpoints verify that relationship.
    # ---------------------------------------------------------------

    normalized = np.zeros_like(
        X,
        dtype=np.float64,
    )

    for j, spec in enumerate(
        pmml["inputs"]
    ):

        normalized[:, j] = (
            X[:, j]
            * spec["scale"]
            + spec["shift"]
        )

        # Verify supplied LinearNorm endpoints when present.
        if len(spec["linear_norms"]) >= 2:

            orig_a, norm_a = (
                spec["linear_norms"][0]
            )

            orig_b, norm_b = (
                spec["linear_norms"][1]
            )

            calculated_a = (
                orig_a * spec["scale"]
                + spec["shift"]
            )

            calculated_b = (
                orig_b * spec["scale"]
                + spec["shift"]
            )

            if not math.isclose(
                calculated_a,
                norm_a,
                rel_tol=1e-7,
                abs_tol=1e-7,
            ):

                raise RuntimeError(
                    f"PMML normalization check failed "
                    f"for {spec['field']} at lower endpoint."
                )

            if not math.isclose(
                calculated_b,
                norm_b,
                rel_tol=1e-7,
                abs_tol=1e-7,
            ):

                raise RuntimeError(
                    f"PMML normalization check failed "
                    f"for {spec['field']} at upper endpoint."
                )

    # ---------------------------------------------------------------
    # Hidden layer
    # ---------------------------------------------------------------

    hidden_layer = pmml["layers"][0]

    if hidden_layer["activation"] != "tanh":

        raise RuntimeError(
            "Historical hidden activation is not tanh."
        )

    hidden = np.zeros(
        (
            X.shape[0],
            len(hidden_layer["neurons"]),
        ),
        dtype=np.float64,
    )

    input_names = [
        str(i)
        for i in range(
            X.shape[1]
        )
    ]

    for j, neuron in enumerate(
        hidden_layer["neurons"]
    ):

        z = np.full(
            X.shape[0],
            neuron["bias"],
            dtype=np.float64,
        )

        for source, weight in (
            neuron["connections"].items()
        ):

            source_index = int(source)

            z += (
                normalized[:, source_index]
                * weight
            )

        hidden[:, j] = np.tanh(z)

    # ---------------------------------------------------------------
    # Output layer
    # ---------------------------------------------------------------

    output_layer = pmml["layers"][1]

    if output_layer["activation"] != "tanh":

        raise RuntimeError(
            "Historical output activation is not tanh."
        )

    outputs = np.zeros(
        (
            X.shape[0],
            len(output_layer["neurons"]),
        ),
        dtype=np.float64,
    )

    hidden_index_by_id = {
        neuron["id"]: i
        for i, neuron in enumerate(
            hidden_layer["neurons"]
        )
    }

    for j, neuron in enumerate(
        output_layer["neurons"]
    ):

        z = np.full(
            X.shape[0],
            neuron["bias"],
            dtype=np.float64,
        )

        for source, weight in (
            neuron["connections"].items()
        ):

            if source not in hidden_index_by_id:

                raise RuntimeError(
                    "Output neuron references unknown "
                    f"hidden neuron {source}."
                )

            source_index = (
                hidden_index_by_id[source]
            )

            z += (
                hidden[:, source_index]
                * weight
            )

        outputs[:, j] = np.tanh(z)

    output_index = {
        neuron["id"]: i
        for i, neuron in enumerate(
            output_layer["neurons"]
        )
    }

    cancer_score = outputs[
        :,
        output_index[
            pmml["cancer_output_neuron"]
        ],
    ]

    noncancer_score = outputs[
        :,
        output_index[
            pmml["noncancer_output_neuron"]
        ],
    ]

    predicted_cancer = (
        cancer_score > noncancer_score
    ).astype(int)

    return (
        normalized,
        hidden,
        outputs,
        cancer_score,
        noncancer_score,
        predicted_cancer,
    )


# =====================================================================
# METRICS
# =====================================================================

def binary_metrics(
    y_true,
    y_pred,
    score,
):

    y_true = np.asarray(
        y_true,
        dtype=int,
    )

    y_pred = np.asarray(
        y_pred,
        dtype=int,
    )

    score = np.asarray(
        score,
        dtype=float,
    )

    tp = int(
        np.sum(
            (y_true == 1)
            & (y_pred == 1)
        )
    )

    tn = int(
        np.sum(
            (y_true == 0)
            & (y_pred == 0)
        )
    )

    fp = int(
        np.sum(
            (y_true == 0)
            & (y_pred == 1)
        )
    )

    fn = int(
        np.sum(
            (y_true == 1)
            & (y_pred == 0)
        )
    )

    sensitivity = (
        tp / (tp + fn)
        if (tp + fn) > 0
        else float("nan")
    )

    specificity = (
        tn / (tn + fp)
        if (tn + fp) > 0
        else float("nan")
    )

    accuracy = (
        (tp + tn)
        / max(
            1,
            len(y_true),
        )
    )

    auc = float("nan")

    try:

        from sklearn.metrics import roc_auc_score

        if len(np.unique(y_true)) == 2:

            auc = float(
                roc_auc_score(
                    y_true,
                    score,
                )
            )

    except Exception:

        pass

    return {
        "auc": auc,
        "accuracy": float(accuracy),
        "sensitivity": float(sensitivity),
        "specificity": float(specificity),
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
    }


# =====================================================================
# MAIN
# =====================================================================

def main():

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ---------------------------------------------------------------
    # Validate required inputs
    # ---------------------------------------------------------------

    required = [
        SUPPLEMENTARY_XLSX,
        PMML_PATH,
        GEO_METADATA,
        COHORT1_XLSX,
        COHORT2_XLSX,
    ]

    for path in required:

        if not path.exists():

            raise RuntimeError(
                f"Required file does not exist:\n{path}"
            )

    # ---------------------------------------------------------------
    # Historical Dataset 2
    # ---------------------------------------------------------------

    historical = load_dataset2()

    # Save the exact recovered historical 14-miRNA matrix.
    historical.to_csv(
        HISTORICAL_MATRIX_OUTPUT,
        index=False,
    )

    print()
    print(
        f"Saved historical matrix:\n"
        f"  {HISTORICAL_MATRIX_OUTPUT}"
    )

    # ---------------------------------------------------------------
    # GEO diagnosis
    # ---------------------------------------------------------------

    geo_lookup = load_diagnosis()

    diagnosis_raw = []
    diagnosis_binary = []
    stages = []
    histologies = []

    missing_geo = []

    for sample_id in historical[
        "sample_id"
    ]:

        record = geo_lookup.get(
            sample_id
        )

        if record is None:

            missing_geo.append(
                sample_id
            )

            diagnosis_raw.append(
                None
            )

            diagnosis_binary.append(
                None
            )

            stages.append(
                None
            )

            histologies.append(
                None
            )

            continue

        diagnosis = str(
            record["diagnosis"]
        ).strip().lower()

        diagnosis_raw.append(
            record["diagnosis"]
        )

        stages.append(
            record.get(
                "stage",
                None,
            )
        )

        histologies.append(
            record.get(
                "histology",
                None,
            )
        )

        if diagnosis == "cancer":

            diagnosis_binary.append(
                1
            )

        elif diagnosis in {
            "benign",
            "borderline",
            "control",
        }:

            diagnosis_binary.append(
                0
            )

        else:

            raise RuntimeError(
                "Unexpected GEO diagnosis "
                f"{record['diagnosis']!r} "
                f"for sample {sample_id}."
            )

    if missing_geo:

        raise RuntimeError(
            "Historical samples missing GEO metadata:\n"
            + "\n".join(
                f"  {x}"
                for x in missing_geo
            )
        )

    historical["diagnosis"] = (
        diagnosis_raw
    )

    historical["y"] = (
        diagnosis_binary
    )

    historical["stage"] = stages
    historical["histology"] = histologies

    # ---------------------------------------------------------------
    # Diagnosis distribution
    # ---------------------------------------------------------------

    print()
    print("=" * 100)
    print("HISTORICAL DIAGNOSIS DISTRIBUTION")
    print("=" * 100)

    counts = (
        historical["diagnosis"]
        .value_counts()
        .to_dict()
    )

    for key in [
        "Cancer",
        "Benign",
        "Borderline",
        "Control",
    ]:

        print(
            f"  {key:<12} "
            f"{counts.get(key, 0)}"
        )

    print()
    print(
        "Binary:"
    )

    print(
        "  Cancer:     "
        f"{int((historical['y'] == 1).sum())}"
    )

    print(
        "  Non-cancer: "
        f"{int((historical['y'] == 0).sum())}"
    )

    # ---------------------------------------------------------------
    # PMML
    # ---------------------------------------------------------------

    pmml = parse_pmml()

    # ---------------------------------------------------------------
    # Forward pass
    # ---------------------------------------------------------------

    X = historical[
        FEATURES
    ].to_numpy(
        dtype=np.float64
    )

    (
        normalized,
        hidden,
        outputs,
        cancer_score,
        noncancer_score,
        predicted_cancer,
    ) = forward_pass(
        X,
        pmml,
    )

    # ---------------------------------------------------------------
    # Metrics
    # ---------------------------------------------------------------

    y = historical[
        "y"
    ].to_numpy(
        dtype=int
    )

    metrics = binary_metrics(
        y,
        predicted_cancer,
        cancer_score,
    )

    historical[
        "cancer_score"
    ] = cancer_score

    historical[
        "noncancer_score"
    ] = noncancer_score

    historical[
        "predicted_cancer"
    ] = predicted_cancer

    historical[
        "correct"
    ] = (
        historical["y"]
        == historical["predicted_cancer"]
    )

    historical.to_csv(
        PREDICTIONS_OUTPUT,
        index=False,
    )

    # ---------------------------------------------------------------
    # Print results
    # ---------------------------------------------------------------

    print()
    print("=" * 100)
    print("ORIGINAL eLIFE PMML HISTORICAL REPRODUCTION")
    print("=" * 100)

    print()
    print(
        "Architecture:"
    )

    print(
        "  14 -> 7 tanh -> 2 tanh"
    )

    print()
    print(
        "Input representation:"
    )

    print(
        "  Supplemental Dataset 2"
    )

    print(
        "  batch-adjusted"
    )

    print(
        "  log10-transformed"
    )

    print()
    print(
        "Historical sample count:"
        f" {len(historical)}"
    )

    print()
    print(
        "Metrics on historical dataset:"
    )

    print(
        f"  AUC:         {metrics['auc']:.6f}"
    )

    print(
        f"  Accuracy:    {metrics['accuracy']:.6f}"
    )

    print(
        f"  Sensitivity: {metrics['sensitivity']:.6f}"
    )

    print(
        f"  Specificity: {metrics['specificity']:.6f}"
    )

    print()
    print(
        "Confusion matrix:"
    )

    print(
        f"  TP={metrics['tp']} "
        f"FN={metrics['fn']}"
    )

    print(
        f"  FP={metrics['fp']} "
        f"TN={metrics['tn']}"
    )

    print()
    print(
        "Output neuron scores:"
    )

    print(
        f"  Cancer neuron range: "
        f"{np.min(cancer_score):.6f} "
        f"to "
        f"{np.max(cancer_score):.6f}"
    )

    print(
        f"  Non-cancer neuron range: "
        f"{np.min(noncancer_score):.6f} "
        f"to "
        f"{np.max(noncancer_score):.6f}"
    )

    # ---------------------------------------------------------------
    # Audit JSON
    # ---------------------------------------------------------------

    audit = {
        "model": {
            "pmml": str(PMML_PATH),
            "model_name": pmml["model_name"],
            "function_name": pmml["function_name"],
            "architecture": "14-7-2",
            "hidden_activation": "tanh",
            "output_activation": "tanh",
            "cancer_output_neuron": "21",
            "noncancer_output_neuron": "22",
        },
        "historical_dataset": {
            "source": str(
                SUPPLEMENTARY_XLSX
            ),
            "sheet": "Supplemental Dataset 2",
            "representation": [
                "batch-adjusted",
                "log10-transformed",
            ],
            "patient_count": int(
                len(historical)
            ),
            "feature_count": 14,
            "e1056_present": (
                "e1056"
                in set(
                    historical["sample_id"]
                )
            ),
        },
        "diagnosis_counts": {
            str(key): int(value)
            for key, value in counts.items()
        },
        "binary_counts": {
            "cancer": int(
                (historical["y"] == 1).sum()
            ),
            "non_cancer": int(
                (historical["y"] == 0).sum()
            ),
        },
        "metrics": metrics,
        "files": {
            "historical_matrix": str(
                HISTORICAL_MATRIX_OUTPUT
            ),
            "predictions": str(
                PREDICTIONS_OUTPUT
            ),
            "audit": str(
                AUDIT_OUTPUT
            ),
        },
    }

    AUDIT_OUTPUT.write_text(
        json.dumps(
            audit,
            indent=2,
        )
    )

    print()
    print("=" * 100)
    print("OUTPUT FILES")
    print("=" * 100)

    print()
    print(
        f"Historical 14-miRNA matrix:\n"
        f"  {HISTORICAL_MATRIX_OUTPUT}"
    )

    print()
    print(
        f"Predictions:\n"
        f"  {PREDICTIONS_OUTPUT}"
    )

    print()
    print(
        f"Audit:\n"
        f"  {AUDIT_OUTPUT}"
    )

    print()
    print("=" * 100)
    print("SUCCESS")
    print("=" * 100)

    print()
    print(
        "NO MODEL TRAINING WAS PERFORMED."
    )

    print()
    print(
        "This reproduction now uses Supplemental Dataset 2 directly"
    )

    print(
        "and independently reconstructs its patient identifiers."
    )

    print()
    print(
        "Important:"
    )

    print(
        "The reported metrics are historical/apparent performance,"
    )

    print(
        "not independent validation."
    )


if __name__ == "__main__":
    main()
