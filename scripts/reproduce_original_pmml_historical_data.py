#!/usr/bin/env python3

"""
Reproduce the original eLife GSE94533 neural network using the
HISTORICAL DATA REPRESENTATION rather than our raw-TPM matrix.

Why this script exists
----------------------

The previous reproduction produced:

    AUC = 0.451682
    sensitivity = 0
    specificity = 1

That strongly suggested that the neural network itself was being
fed the wrong representation of the data.

The PMML normalization anchors are in ranges such as:

    hsa-miR-23b-3p:
        0.1816 -> 1.8782

and the eLife supplementary workbook identifies Supplemental Dataset 2
as:

    "Batch adjusted, log10-transformed miRNA expression data"

Those ranges are consistent with the PMML inputs.

Therefore the historical network should be evaluated on Supplemental
Dataset 2 rather than the raw TPM matrix.

This script:

    Supplemental Dataset 2
            |
            v
    179 historical sample positions
            |
            v
    map positions to actual GSE94533 sample IDs
            |
            v
    exact 14 PMML features
            |
            v
    PMML normalization
            |
            v
    14 -> 7 tanh -> 2 tanh
            |
            v
    historical predictions

NO MODEL TRAINING IS PERFORMED.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


# ===================================================================
# PATHS
# ===================================================================

ROOT = Path.home() / "ovarian-ai"

SUPPLEMENTARY_XLSX = (
    ROOT
    / "data"
    / "external"
    / "elife_28932"
    / "elife_poa_e28932_Supplementary_File_6.xlsx"
)

PMML_AUDIT = (
    ROOT
    / "reports"
    / "tables"
    / "original_pmml_deep_audit.json"
)

CURRENT_METADATA = (
    ROOT
    / "data"
    / "processed"
    / "gse94533_metadata.csv"
)

OUTPUT_DIR = (
    ROOT
    / "reports"
    / "tables"
)

HISTORICAL_MATRIX_OUTPUT = (
    ROOT
    / "data"
    / "processed"
    / "gse94533_historical_14mirna.csv.gz"
)

PREDICTIONS_OUTPUT = (
    OUTPUT_DIR
    / "original_pmml_historical_predictions.csv"
)

RESULTS_OUTPUT = (
    OUTPUT_DIR
    / "original_pmml_historical_results.json"
)


# ===================================================================
# EXACT PMML FEATURES
# ===================================================================

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


# ===================================================================
# CURRENT WORKBOOK SAMPLE ORDER
# ===================================================================

def normalize_id(value):

    if pd.isna(value):
        return ""

    s = str(value).strip().lower()

    # Numeric Excel identifiers such as 50217.0
    if s.endswith(".0"):
        try:
            s = str(int(float(s)))
        except Exception:
            pass

    # E-1001 -> e1001
    s = (
        s
        .replace("-", "")
        .replace("_", "")
        .replace(" ", "")
    )

    return s


def get_processed_sample_order():

    cohort1 = (
        ROOT
        / "data"
        / "raw"
        / "GSE94533_Processed_file_Cohort1.xlsx"
    )

    cohort2 = (
        ROOT
        / "data"
        / "raw"
        / "GSE94533_Processed_file_Cohort2.xlsx"
    )

    print()
    print("=" * 100)
    print("RECOVERING ORIGINAL SAMPLE ORDER")
    print("=" * 100)

    # ---------------------------------------------------------------
    # Cohort 2
    #
    # The supplementary historical dataset begins with the numeric
    # sample IDs from this workbook.
    # ---------------------------------------------------------------

    df2_header = pd.read_excel(
        cohort2,
        sheet_name="mirna_tpm",
        header=2,
        nrows=0,
        engine="openpyxl",
    )

    cohort2_ids = []

    for column in df2_header.columns:

        normalized = normalize_id(column)

        if normalized.isdigit():

            cohort2_ids.append(
                str(
                    int(
                        float(
                            column
                        )
                    )
                )
            )

    # ---------------------------------------------------------------
    # Cohort 1
    # ---------------------------------------------------------------

    df1_header = pd.read_excel(
        cohort1,
        sheet_name="mirna_tpm",
        header=2,
        nrows=0,
        engine="openpyxl",
    )

    cohort1_ids = []

    for column in df1_header.columns:

        normalized = normalize_id(column)

        if normalized.startswith("e"):

            cohort1_ids.append(
                normalized
            )

    # The historical primary cohort contains 120 numeric samples and
    # 60 E- samples, with E-1056 excluded.
    #
    # Therefore:
    #
    #     120 + 59 = 179
    #
    cohort1_ids = [
        x
        for x in cohort1_ids
        if x != "e1056"
    ]

    combined = (
        cohort2_ids
        + cohort1_ids
    )

    print()
    print(
        f"Cohort 2 sample IDs: {len(cohort2_ids)}"
    )

    print(
        f"Cohort 1 sample IDs after E-1056 removal: "
        f"{len(cohort1_ids)}"
    )

    print(
        f"Combined historical sample IDs: {len(combined)}"
    )

    if len(combined) != 179:

        raise RuntimeError(
            "Expected exactly 179 historical sample IDs."
        )

    if len(set(combined)) != 179:

        raise RuntimeError(
            "Historical sample IDs are not unique."
        )

    return combined


# ===================================================================
# LOAD SUPPLEMENTAL DATASET 2
# ===================================================================


def load_supplemental_dataset_2():

    print()
    print("=" * 100)
    print("READING SUPPLEMENTAL DATASET 2")
    print("=" * 100)

    if not SUPPLEMENTARY_XLSX.exists():

        raise RuntimeError(
            f"Supplementary workbook not found:\n"
            f"{SUPPLEMENTARY_XLSX}"
        )

    raw = pd.read_excel(
        SUPPLEMENTARY_XLSX,
        sheet_name="Supplemental Dataset 2",
        header=None,
        engine="openpyxl",
    )

    print()
    print(
        f"Raw shape: {raw.shape}"
    )

    # ---------------------------------------------------------------
    # IMPORTANT
    #
    # Supplemental Dataset 2 does NOT necessarily use:
    #
    #     1, 2, 3, ..., 179
    #
    # as its column headers.
    #
    # Depending on how Excel/pandas interprets the supplementary
    # workbook, the sample header row can contain the REAL GSE IDs:
    #
    #     e1036
    #     e1048
    #     ...
    #     50217
    #     80675
    #
    # Therefore we detect the sample-ID row rather than assuming
    # numeric sample positions.
    # ---------------------------------------------------------------

    def local_normalize_id(value):

        if pd.isna(value):
            return ""

        s = str(value).strip().lower()

        # Handle Excel numeric values such as 50217.0.
        if s.endswith(".0"):

            try:
                s = str(int(float(s)))
            except Exception:
                pass

        s = (
            s
            .replace("-", "")
            .replace("_", "")
            .replace(" ", "")
        )

        return s

    def looks_like_sample_id(value):

        s = local_normalize_id(value)

        if not s:
            return False

        # e1036-style Cohort 1 IDs.
        if s.startswith("e") and s[1:].isdigit():
            return True

        # Numeric Cohort 2 IDs.
        if s.isdigit():
            return True

        return False

    # ---------------------------------------------------------------
    # Detect the row containing the 179 actual sample IDs.
    #
    # We search the first several rows and count how many cells after
    # column zero look like real sample identifiers.
    # ---------------------------------------------------------------

    sample_header_candidates = []

    search_rows = min(
        20,
        len(raw),
    )

    for row_idx in range(
        search_rows
    ):

        values = raw.iloc[
            row_idx,
            1:
        ].tolist()

        ids_here = [
            local_normalize_id(x)
            for x in values
            if looks_like_sample_id(x)
        ]

        unique_ids_here = set(
            ids_here
        )

        sample_header_candidates.append(
            (
                row_idx,
                len(ids_here),
                len(unique_ids_here),
                ids_here,
            )
        )

    # Prefer a row with 179 distinct sample IDs.
    exact_candidates = [
        x
        for x in sample_header_candidates
        if x[2] == 179
    ]

    if exact_candidates:

        sample_header_row = exact_candidates[0][0]

    else:

        # Otherwise choose the row with the greatest number of
        # unique sample-like identifiers.
        ranked = sorted(
            sample_header_candidates,
            key=lambda x: (
                x[2],
                x[1],
            ),
            reverse=True,
        )

        if not ranked or ranked[0][2] < 100:

            print()
            print(
                "Could not identify the historical sample-ID row."
            )

            print()
            print(
                "Candidate rows:"
            )

            for row_idx, count, unique_count, ids_here in ranked:
                print(
                    f"  row={row_idx:2d} "
                    f"sample_like={count:3d} "
                    f"unique={unique_count:3d} "
                    f"first_ids={ids_here[:10]}"
                )

            raise RuntimeError(
                "Could not safely identify the Supplemental Dataset 2 "
                "sample-ID header row."
            )

        sample_header_row = ranked[0][0]

    # ---------------------------------------------------------------
    # Recover the actual sample IDs in the order they occur in the
    # supplementary workbook.
    # ---------------------------------------------------------------

    sample_ids = []

    sample_column_indices = []

    for col_idx in range(
        1,
        raw.shape[1],
    ):

        value = raw.iloc[
            sample_header_row,
            col_idx
        ]

        if looks_like_sample_id(value):

            normalized = local_normalize_id(
                value
            )

            sample_ids.append(
                normalized
            )

            sample_column_indices.append(
                col_idx
            )

    print()
    print(
        f"Detected sample-ID header row: {sample_header_row}"
    )

    print(
        f"Detected sample columns: {len(sample_ids)}"
    )

    print()
    print(
        "First sample IDs:"
    )

    for sample_id in sample_ids[:15]:

        print(
            f"  {sample_id}"
        )

    print()
    print(
        "Last sample IDs:"
    )

    for sample_id in sample_ids[-15:]:

        print(
            f"  {sample_id}"
        )

    # ---------------------------------------------------------------
    # IMPORTANT HISTORICAL DATASET DETAIL
    #
    # Supplemental Dataset 2 uses POSITIONAL sample numbers:
    #
    #     1, 2, 3, ..., 179
    #
    # These are NOT the actual GSE94533 sample identifiers.
    #
    # The original supplementary Dataset 1 makes clear that the
    # historical 179-sample matrix consists of:
    #
    #     first 120 = Cohort 2 numeric sample IDs
    #     next 59   = Cohort 1 E-#### IDs
    #
    # with E-1056 excluded from the historical 179-sample analysis.
    #
    # Therefore:
    #
    #     Dataset-2 position
    #             |
    #             v
    #     recovered original sample order
    #             |
    #             v
    #     actual GSE94533 sample ID
    #
    # We MUST make this translation before doing any patient-level
    # metadata matching.
    # ---------------------------------------------------------------

    positional_ids = [
        str(i)
        for i in range(
            1,
            180,
        )
    ]

    if sample_ids == positional_ids:

        print()
        print(
            "Detected positional sample headers 1..179."
        )

        print(
            "Translating historical positions to actual "
            "GSE94533 sample IDs..."
        )

        recovered_order = get_processed_sample_order()

        if len(recovered_order) != 179:

            raise RuntimeError(
                "Recovered original sample order contains "
                f"{len(recovered_order)} samples; expected 179."
            )

        sample_ids = [
            normalize_id(x)
            for x in recovered_order
        ]

    # ============================================================
    # HISTORICAL DATASET 2 SAMPLE-ID CORRECTION
    # ============================================================
    #
    # Supplemental Dataset 2 has a mixed header:
    #
    #     columns 1..120   -> 1, 2, 3, ..., 120
    #     columns 121..179 -> e1036, e1048, ..., e1021
    #
    # The first 120 values are historical POSITION NUMBERS.
    # They are not patient identifiers.
    #
    # The original sample order recovered from the processed
    # workbooks gives the real patient ID corresponding to each
    # historical position.
    #
    # Therefore:
    #
    #     Dataset 2 position
    #            |
    #            v
    #     original sample order
    #            |
    #            v
    #     real GSE94533 patient ID
    #
    # We only translate the first 120 positions.
    # The final 59 IDs are already explicit.
    # ============================================================

    if (
        len(sample_ids) == 179
        and all(
            str(sample_ids[i]).strip() == str(i + 1)
            for i in range(120)
        )
    ):

        print()
        print(
            "Detected mixed Dataset 2 sample-header format."
        )

        print(
            "Columns 1-120 are positional sample numbers."
        )

        print(
            "Columns 121-179 already contain explicit IDs."
        )

        print()
        print(
            "Translating historical positions into "
            "actual GSE94533 sample IDs..."
        )

        recovered_order = get_processed_sample_order()

        if len(recovered_order) != 179:
            raise RuntimeError(
                "Recovered original sample order has "
                f"{len(recovered_order)} samples; expected 179."
            )

        recovered_order = [
            normalize_id(x)
            for x in recovered_order
        ]

        # First 120 entries must correspond to Cohort 2.
        cohort2_ids = recovered_order[:120]

        # Final 59 entries in Dataset 2 are explicit Cohort 1 IDs.
        explicit_cohort1_ids = [
            normalize_id(x)
            for x in sample_ids[120:]
        ]

        # The recovered historical order also contains the expected
        # Cohort 1 sequence.
        expected_cohort1_ids = recovered_order[120:]

        if explicit_cohort1_ids != expected_cohort1_ids:

            print()
            print(
                "ERROR: Dataset 2 Cohort 1 ordering does not "
                "match the recovered historical order."
            )

            print()
            print(
                "Dataset 2 explicit Cohort 1 IDs:"
            )

            for position, value in enumerate(
                explicit_cohort1_ids,
                start=121,
            ):
                print(
                    f"  {position:3d} -> {value}"
                )

            print()
            print(
                "Recovered Cohort 1 IDs:"
            )

            for position, value in enumerate(
                expected_cohort1_ids,
                start=121,
            ):
                print(
                    f"  {position:3d} -> {value}"
                )

            raise RuntimeError(
                "Cohort 1 IDs in Supplemental Dataset 2 "
                "do not agree with the recovered sample order."
            )

        # Build the complete real-ID order.
        sample_ids = (
            cohort2_ids
            + explicit_cohort1_ids
        )

        # --------------------------------------------------------
        # Strong integrity checks
        # --------------------------------------------------------

        if len(sample_ids) != 179:
            raise RuntimeError(
                "Translated historical sample list does not "
                "contain exactly 179 samples."
            )

        if len(set(sample_ids)) != 179:
            raise RuntimeError(
                "Translated historical sample IDs are not unique."
            )

        print()
        print(
            "Historical position -> actual sample ID:"
        )

        for position in range(1, 11):

            print(
                f"  {position:3d} -> "
                f"{sample_ids[position - 1]}"
            )

        print(
            "  ..."
        )

        for position in range(171, 180):

            print(
                f"  {position:3d} -> "
                f"{sample_ids[position - 1]}"
            )

        print()
        print(
            f"Position 1   -> {sample_ids[0]}"
        )

        print(
            f"Position 120 -> {sample_ids[119]}"
        )

        print(
            f"Position 121 -> {sample_ids[120]}"
        )

        print(
            f"Position 179 -> {sample_ids[178]}"
        )

        print()
        print(
            "Historical position -> real sample ID mapping:"
        )

        for position, sample_id in list(
            enumerate(
                sample_ids,
                start=1,
            )
        )[:15]:

            print(
                f"  {position:3d} -> {sample_id}"
            )

        print(
            "  ..."
        )

        for position, sample_id in list(
            enumerate(
                sample_ids,
                start=1,
            )
        )[-15:]:

            print(
                f"  {position:3d} -> {sample_id}"
            )

        print()
        print(
            "Position 1:"
            f"  {sample_ids[0]}"
        )

        print(
            "Position 120:"
            f"  {sample_ids[119]}"
        )

        print(
            "Position 121:"
            f"  {sample_ids[120]}"
        )

        print(
            "Position 179:"
            f"  {sample_ids[178]}"
        )

    else:

        print()
        print(
            "Supplemental Dataset 2 appears to contain "
            "explicit sample identifiers rather than "
            "positional 1..179 headers."
        )

    # ---------------------------------------------------------------
    # Exactly 179 historical patients are expected.
    # ---------------------------------------------------------------

    if len(sample_ids) != 179:

        raise RuntimeError(
            "Supplemental Dataset 2 yielded "
            f"{len(sample_ids)} sample IDs; expected 179."
        )

    if len(set(sample_ids)) != 179:

        duplicates = (
            pd.Series(
                sample_ids
            )
            .value_counts()
        )

        duplicates = duplicates[
            duplicates > 1
        ]

        raise RuntimeError(
            "Supplemental Dataset 2 contains duplicate sample IDs:\n"
            + duplicates.to_string()
        )

    # ---------------------------------------------------------------
    # Find the rows corresponding to the exact 14 PMML features.
    #
    # We do NOT assume that all remaining rows are miRNAs.
    # We explicitly locate each requested feature by name.
    # ---------------------------------------------------------------

    first_column = (
        raw.iloc[
            :,
            0
        ]
        .astype(str)
        .str.strip()
    )

    feature_row_map = {}

    for feature in FEATURES:

        matches = np.where(
            first_column.eq(
                feature
            )
        )[0]

        if len(matches) == 0:

            raise RuntimeError(
                f"Could not find PMML feature "
                f"{feature!r} in Supplemental Dataset 2."
            )

        if len(matches) > 1:

            raise RuntimeError(
                f"Feature {feature!r} appears "
                f"{len(matches)} times in Supplemental Dataset 2."
            )

        feature_row_map[
            feature
        ] = int(
            matches[0]
        )

    print()
    print(
        "Exact PMML feature rows found:"
    )

    for feature in FEATURES:

        print(
            f"  {feature}: row "
            f"{feature_row_map[feature]}"
        )

    # ---------------------------------------------------------------
    # Build feature x sample values using the detected sample columns.
    # ---------------------------------------------------------------

    selected_rows = []

    for feature in FEATURES:

        row_idx = feature_row_map[
            feature
        ]

        values = []

        for col_idx in sample_column_indices:

            values.append(
                raw.iloc[
                    row_idx,
                    col_idx
                ]
            )

        selected_rows.append(
            values
        )

    expression = pd.DataFrame(
        selected_rows,
        index=FEATURES,
        columns=sample_ids,
    )

    expression = expression.apply(
        pd.to_numeric,
        errors="coerce",
    )

    # ---------------------------------------------------------------
    # Every selected feature must contain finite numeric values.
    # ---------------------------------------------------------------

    values = expression.to_numpy(
        dtype=float
    )

    if not np.isfinite(
        values
    ).all():

        bad = np.argwhere(
            ~np.isfinite(
                values
            )
        )

        first_bad_row = int(
            bad[0][0]
        )

        first_bad_col = int(
            bad[0][1]
        )

        raise RuntimeError(
            "Supplemental Dataset 2 contains non-finite "
            "values in the selected PMML matrix.\n"
            f"Feature: {FEATURES[first_bad_row]}\n"
            f"Sample: {sample_ids[first_bad_col]}"
        )

    # ---------------------------------------------------------------
    # Transpose:
    #
    #     sample x feature
    # ---------------------------------------------------------------

    selected = expression.T.copy()

    selected.columns = FEATURES

    # Preserve the historical positional column as a provenance field.
    selected[
        "historical_sample_number"
    ] = np.arange(
        1,
        len(selected) + 1,
    )

    # `historical_sample_number` is the original positional column
    # from Supplemental Dataset 2.
    #
    # `sample_id` below is the actual GSE94533 identifier recovered
    # from the original cohort ordering.
    #
    # Keeping both fields prevents us from confusing the two concepts.

    # Actual GSE94533 sample ID.
    selected[
        "sample_id"
    ] = selected.index.astype(
        str
    )

    # ---------------------------------------------------------------
    # Put identifiers first.
    # ---------------------------------------------------------------

    selected = selected[
        [
            "sample_id",
            "historical_sample_number",
        ]
        + FEATURES
    ]

    print()
    print(
        "Historical representation:"
    )

    print(
        "  Batch adjusted"
    )

    print(
        "  log10 transformed"
    )

    print(
        f"  Samples: {len(selected)}"
    )

    print(
        "  PMML features: 14"
    )

    return selected



def verify_against_processed_data(
    historical,
    sample_ids,
):

    print()
    print("=" * 100)
    print("VERIFYING HISTORICAL SAMPLE IDENTIFIERS")
    print("=" * 100)

    # ---------------------------------------------------------------
    # Expected GSE94533 IDs reconstructed from the processed
    # Cohort 1 + Cohort 2 workbooks.
    # ---------------------------------------------------------------

    expected_ids = set(
        get_processed_sample_order()
    )

    observed_ids = set(
        normalize_id(x)
        for x in historical[
            "sample_id"
        ]
    )

    print()
    print(
        f"Expected processed-workbook IDs: {len(expected_ids)}"
    )

    print(
        f"Observed historical IDs:          {len(observed_ids)}"
    )

    missing = sorted(
        expected_ids
        - observed_ids
    )

    extra = sorted(
        observed_ids
        - expected_ids
    )

    if missing:

        raise RuntimeError(
            "Historical Dataset 2 is missing sample IDs:\n"
            + "\n".join(
                f"  {x}"
                for x in missing
            )
        )

    if extra:

        raise RuntimeError(
            "Historical Dataset 2 contains unexpected sample IDs:\n"
            + "\n".join(
                f"  {x}"
                for x in extra
            )
        )

    # ---------------------------------------------------------------
    # Confirm the order supplied by this function is identical to
    # the actual order from the historical workbook.
    # ---------------------------------------------------------------

    if list(
        historical[
            "sample_id"
        ]
    ) != list(
        sample_ids
    ):

        raise RuntimeError(
            "Internal sample-ID order mismatch."
        )

    print()
    print(
        "Sample-ID set check: PASS"
    )

    print(
        "Sample-ID uniqueness: PASS"
    )

    print(
        "Sample-ID order preservation: PASS"
    )

    # ---------------------------------------------------------------
    # We deliberately DO NOT require the first sample to be 50217.
    #
    # Different supplementary datasets can store the same 179
    # patients in different orders.
    # ---------------------------------------------------------------

    print()
    print(
        "No artificial first-sample assumption applied."
    )

    return historical


def load_pmml_parameters():

    print()
    print("=" * 100)
    print("LOADING PMML PARAMETERS")
    print("=" * 100)

    audit = json.loads(
        PMML_AUDIT.read_text()
    )

    nn = audit[
        "neural_network"
    ]

    # ---------------------------------------------------------------
    # Input normalization.
    # ---------------------------------------------------------------

    normalization = {}

    for feature, neural_input in zip(
        FEATURES,
        nn[
            "neural_inputs"
        ],
    ):

        derived = neural_input[
            "expressions"
        ][0]

        norm = derived[
            "children"
        ][0]

        attrs = norm[
            "attributes"
        ]

        normalization[
            feature
        ] = {
            "shift": float(
                attrs[
                    "shift"
                ]
            ),
            "scale": float(
                attrs[
                    "scale"
                ]
            ),
        }

    # ---------------------------------------------------------------
    # Hidden layer
    # ---------------------------------------------------------------

    hidden = nn[
        "neural_layers"
    ][0]

    output = nn[
        "neural_layers"
    ][1]

    hidden_weights = np.zeros(
        (
            7,
            14,
        ),
        dtype=float,
    )

    hidden_bias = np.zeros(
        7,
        dtype=float,
    )

    for i, neuron in enumerate(
        hidden[
            "neurons"
        ]
    ):

        hidden_bias[i] = float(
            neuron[
                "bias"
            ]
        )

        for connection in neuron[
            "connections"
        ]:

            j = int(
                connection[
                    "from"
                ]
            )

            hidden_weights[
                i,
                j
            ] = float(
                connection[
                    "weight"
                ]
            )

    # ---------------------------------------------------------------
    # Output layer
    # ---------------------------------------------------------------

    output_weights = np.zeros(
        (
            2,
            7,
        ),
        dtype=float,
    )

    output_bias = np.zeros(
        2,
        dtype=float,
    )

    hidden_ids = [
        "14",
        "15",
        "16",
        "17",
        "18",
        "19",
        "20",
    ]

    for i, neuron in enumerate(
        output[
            "neurons"
        ]
    ):

        output_bias[i] = float(
            neuron[
                "bias"
            ]
        )

        for connection in neuron[
            "connections"
        ]:

            hidden_id = str(
                connection[
                    "from"
                ]
            )

            j = hidden_ids.index(
                hidden_id
            )

            output_weights[
                i,
                j
            ] = float(
                connection[
                    "weight"
                ]
            )

    return (
        normalization,
        hidden_weights,
        hidden_bias,
        output_weights,
        output_bias,
    )


# ===================================================================
# FORWARD PASS
# ===================================================================

def forward_pass(
    data,
    normalization,
    hidden_weights,
    hidden_bias,
    output_weights,
    output_bias,
):

    X = data[
        FEATURES
    ].to_numpy(
        dtype=float
    )

    # ---------------------------------------------------------------
    # PMML NormContinuous:
    #
    # normalized = x * scale + shift
    # ---------------------------------------------------------------

    X_norm = np.zeros_like(
        X,
        dtype=float,
    )

    for j, feature in enumerate(
        FEATURES
    ):

        scale = normalization[
            feature
        ][
            "scale"
        ]

        shift = normalization[
            feature
        ][
            "shift"
        ]

        X_norm[
            :,
            j
        ] = (
            X[
                :,
                j
            ]
            * scale
            + shift
        )

    # ---------------------------------------------------------------
    # Hidden layer.
    # ---------------------------------------------------------------

    H_linear = (
        X_norm
        @ hidden_weights.T
        + hidden_bias
    )

    H = np.tanh(
        H_linear
    )

    # ---------------------------------------------------------------
    # Output layer.
    # ---------------------------------------------------------------

    O_linear = (
        H
        @ output_weights.T
        + output_bias
    )

    O = np.tanh(
        O_linear
    )

    return (
        X_norm,
        H_linear,
        H,
        O_linear,
        O,
    )


# ===================================================================
# METRICS
# ===================================================================

def auc_rank(
    y_true,
    scores,
):

    y_true = np.asarray(
        y_true,
        dtype=int,
    )

    scores = np.asarray(
        scores,
        dtype=float,
    )

    positive = scores[
        y_true == 1
    ]

    negative = scores[
        y_true == 0
    ]

    wins = 0.0

    for p in positive:

        wins += np.sum(
            p > negative
        )

        wins += 0.5 * np.sum(
            p == negative
        )

    return float(
        wins
        / (
            len(positive)
            * len(negative)
        )
    )


def classification_metrics(
    y_true,
    y_pred,
):

    y_true = np.asarray(
        y_true,
        dtype=int,
    )

    y_pred = np.asarray(
        y_pred,
        dtype=int,
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

    return {
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
        "accuracy": (
            tp + tn
        )
        / len(y_true),
        "sensitivity": (
            tp
            / (tp + fn)
        )
        if tp + fn
        else None,
        "specificity": (
            tn
            / (tn + fp)
        )
        if tn + fp
        else None,
    }


# ===================================================================
# MAIN
# ===================================================================

def main():

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("=" * 100)
    print(
        "HISTORICAL eLIFE GSE94533 PMML REPRODUCTION"
    )
    print("=" * 100)

    # ---------------------------------------------------------------
    # Recover exact historical sample order.
    # ---------------------------------------------------------------

    sample_ids = (
        get_processed_sample_order()
    )

    # ---------------------------------------------------------------
    # Load historical data representation.
    # ---------------------------------------------------------------

    historical = (
        load_supplemental_dataset_2()
    )

    historical = (
        verify_against_processed_data(
            historical,
            sample_ids,
        )
    )

    # ---------------------------------------------------------------
    # Add actual sample IDs.
    # ---------------------------------------------------------------

    historical[
        "sample_id"
    ] = historical[
        "sample_id"
    ].map(
        normalize_id
    )

    # ---------------------------------------------------------------
    # Load current metadata.
    # ---------------------------------------------------------------

    metadata = pd.read_csv(
        CURRENT_METADATA
    )

    metadata[
        "sample_id"
    ] = (
        metadata[
            "sample_id"
        ]
        .astype(str)
        .map(
            normalize_id
        )
    )

    # ---------------------------------------------------------------
    # Normalize historical IDs.
    # ---------------------------------------------------------------

    historical[
        "sample_id"
    ] = historical[
        "sample_id"
    ].map(
        normalize_id
    )

    # ---------------------------------------------------------------
    # Join diagnosis.
    # ---------------------------------------------------------------

    keep_metadata_columns = [
        "sample_id",
        "diagnosis",
        "binary_label",
    ]

    available = [
        col
        for col in keep_metadata_columns
        if col in metadata.columns
    ]

    if set(
        [
            "sample_id",
            "binary_label",
        ]
    ) - set(available):

        raise RuntimeError(
            "Current metadata does not contain expected label columns."
        )

    merged = historical.merge(
        metadata[
            available
        ],
        on="sample_id",
        how="left",
        validate="one_to_one",
    )

    if merged[
        "binary_label"
    ].isna().any():

        missing_ids = merged.loc[
            merged[
                "binary_label"
            ].isna(),
            "sample_id",
        ].tolist()

        raise RuntimeError(
            "Could not map metadata for:\n"
            + "\n".join(
                missing_ids
            )
        )

    # ---------------------------------------------------------------
    # Remove metadata from feature matrix only where appropriate.
    # ---------------------------------------------------------------

    # Save historical 14-miRNA matrix first.
    historical[
        [
            "sample_id",
            "historical_sample_number",
        ]
        + FEATURES
    ].to_csv(
        HISTORICAL_MATRIX_OUTPUT,
        index=False,
    )

    # ---------------------------------------------------------------
    # Run PMML.
    # ---------------------------------------------------------------

    (
        normalization,
        hidden_weights,
        hidden_bias,
        output_weights,
        output_bias,
    ) = load_pmml_parameters()

    (
        X_norm,
        H_linear,
        H,
        O_linear,
        O,
    ) = forward_pass(
        merged,
        normalization,
        hidden_weights,
        hidden_bias,
        output_weights,
        output_bias,
    )

    cancer_score = O[
        :,
        0
    ]

    noncancer_score = O[
        :,
        1
    ]

    # PMML output neuron 21 = cancer.
    predicted = (
        cancer_score
        > noncancer_score
    ).astype(
        int
    )

    y_true = merged[
        "binary_label"
    ].astype(
        int
    ).to_numpy()

    auc = auc_rank(
        y_true,
        cancer_score,
    )

    metrics = classification_metrics(
        y_true,
        predicted,
    )

    # ---------------------------------------------------------------
    # Save predictions.
    # ---------------------------------------------------------------

    predictions = pd.DataFrame(
        {
            "sample_id": merged[
                "sample_id"
            ],
            "diagnosis": merged[
                "diagnosis"
            ],
            "historical_sample_number": merged[
                "historical_sample_number"
            ],
            "true_binary_label": y_true,
            "cancer_output_21": cancer_score,
            "noncancer_output_22": noncancer_score,
            "output_margin": (
                cancer_score
                - noncancer_score
            ),
            "predicted_binary_label": predicted,
            "predicted_class": np.where(
                predicted == 1,
                "Ovarian Cancer",
                "Controls or Borderline",
            ),
        }
    )

    predictions.to_csv(
        PREDICTIONS_OUTPUT,
        index=False,
    )

    # ---------------------------------------------------------------
    # Output results.
    # ---------------------------------------------------------------

    results = {

        "dataset": "GSE94533",

        "representation": (
            "eLife Supplemental Dataset 2: "
            "batch-adjusted log10-transformed miRNA expression"
        ),

        "n_samples": int(
            len(merged)
        ),

        "n_features": 14,

        "model": {
            "architecture": "14-7-2",
            "hidden_activation": "tanh",
            "output_activation": "tanh",
            "cancer_neuron": "21",
            "noncancer_neuron": "22",
        },

        "results": {
            "auc": auc,
            **metrics,
        },

        "diagnosis_counts": (
            merged[
                "diagnosis"
            ]
            .value_counts()
            .to_dict()
        ),

        "verification": {
            "historical_sample_count_is_179": (
                len(merged) == 179
            ),
            "features_are_exact_14": (
                len(FEATURES) == 14
            ),
            "all_inputs_finite": bool(
                np.isfinite(
                    merged[
                        FEATURES
                    ].to_numpy(
                        dtype=float
                    )
                ).all()
            ),
            "all_predictions_finite": bool(
                np.isfinite(
                    O
                ).all()
            ),
            "unique_sample_ids": (
                merged[
                    "sample_id"
                ].nunique()
                == 179
            ),
        },

        "interpretation": (
            "These are reproduction/resubstitution results using the "
            "historical data representation. They are not independent "
            "clinical validation."
        ),
    }

    RESULTS_OUTPUT.write_text(
        json.dumps(
            results,
            indent=2,
        )
    )

    # ---------------------------------------------------------------
    # Terminal output.
    # ---------------------------------------------------------------

    print()
    print("=" * 100)
    print(
        "HISTORICAL REPRODUCTION COMPLETE"
    )
    print("=" * 100)

    print()
    print(
        "Historical representation:"
    )

    print(
        "  Supplemental Dataset 2"
    )

    print(
        "  Batch-adjusted"
    )

    print(
        "  log10-transformed"
    )

    print()
    print(
        f"Patients: {len(merged)}"
    )

    print(
        "Features: 14"
    )

    print()
    print(
        "Sample mapping:"
    )

    print(
        f"  Unique IDs: "
        f"{merged['sample_id'].nunique()}"
    )

    print()
    print(
        "Diagnosis counts:"
    )

    for diagnosis, count in (
        merged[
            "diagnosis"
        ]
        .value_counts()
        .items()
    ):

        print(
            f"  {diagnosis}: {count}"
        )

    print()
    print(
        "============================================"
    )

    print(
        "HISTORICAL PMML PERFORMANCE"
    )

    print(
        "============================================"
    )

    print()
    print(
        f"AUC:         {auc:.6f}"
    )

    print(
        f"Accuracy:    {metrics['accuracy']:.6f}"
    )

    print(
        f"Sensitivity: {metrics['sensitivity']:.6f}"
        if metrics["sensitivity"] is not None
        else "Sensitivity: unavailable"
    )

    print(
        f"Specificity: {metrics['specificity']:.6f}"
        if metrics["specificity"] is not None
        else "Specificity: unavailable"
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
        "Files:"
    )

    print(
        f"  Historical matrix:"
        f"\n    {HISTORICAL_MATRIX_OUTPUT}"
    )

    print(
        f"  Predictions:"
        f"\n    {PREDICTIONS_OUTPUT}"
    )

    print(
        f"  Results:"
        f"\n    {RESULTS_OUTPUT}"
    )

    print()
    print("=" * 100)
    print(
        "NO MODEL TRAINING PERFORMED"
    )
    print("=" * 100)

    print(
        """
This run uses the historical expression representation required by
the original PMML.

The previous raw-TPM experiment should therefore NOT be used as the
historical benchmark.

The important comparison now becomes:

    historical published PMML
              vs
    modern machine-learning baseline
              vs
    multimodal model

Only the latter two will involve training.
"""
    )


if __name__ == "__main__":
    main()
