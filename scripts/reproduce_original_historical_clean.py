#!/usr/bin/env python3

"""
CLEAN HISTORICAL REPRODUCTION OF THE ORIGINAL eLIFE GSE94533 MODEL

This is a completely independent reproducer.

It does NOT use:
    reproduce_original_historical_standalone.py
    or any of its broken header-validation logic.

Historical data:
    eLife Supplementary Dataset 2
    "Batch adjusted, log10-transformed miRNA expression data"

Historical model:
    14 inputs
        ->
    7 tanh hidden neurons
        ->
    2 tanh output neurons

Output neurons:
    21 = Ovarian Cancer
    22 = Controls or Borderline

No training is performed.
"""


from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


# =====================================================================
# PATHS
# =====================================================================

ROOT = Path.home() / "ovarian-ai"

RAW_DIR = ROOT / "data" / "raw"
INTERIM_DIR = ROOT / "data" / "interim"
PROCESSED_DIR = ROOT / "data" / "processed"
EXTERNAL_DIR = ROOT / "data" / "external"
REPORT_DIR = ROOT / "reports" / "tables"

SUPPLEMENTARY_XLSX = (
    EXTERNAL_DIR
    / "elife_28932"
    / "elife_poa_e28932_Supplementary_File_6.xlsx"
)

METADATA_PATH = (
    INTERIM_DIR
    / "gse94533_geo_metadata_full.csv"
)

COHORT1_XLSX = (
    RAW_DIR
    / "GSE94533_Processed_file_Cohort1.xlsx"
)

COHORT2_XLSX = (
    RAW_DIR
    / "GSE94533_Processed_file_Cohort2.xlsx"
)

PREDICTIONS_OUTPUT = (
    REPORT_DIR
    / "original_historical_clean_predictions.csv"
)

RESULTS_OUTPUT = (
    REPORT_DIR
    / "original_historical_clean_results.json"
)

TRACE_OUTPUT = (
    REPORT_DIR
    / "original_historical_clean_forward_trace.csv"
)

MAPPING_OUTPUT = (
    REPORT_DIR
    / "original_historical_clean_sample_mapping.csv"
)


# =====================================================================
# EXACT 14 PMML FEATURES
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
# EXACT PMML INPUT NORMALIZATION
#
# PMML:
#
#     normalized = raw * scale + shift
#
# These values were recovered directly from Source-code1.xml.
# =====================================================================

SHIFT = np.array([
    -1.07027802107246e-001,
    -7.95408204565984e-001,
     3.24648647984534e-003,
    -1.89351735084109e+000,
    -7.21089703778619e-002,
     6.30181709930541e-003,
    -2.57784834438645e-002,
    -9.78764131194626e-002,
    -6.71320107138154e-001,
    -4.71454485707025e-001,
    -3.14034228530097e-001,
     5.66695208923269e-004,
    -4.35059758202962e-001,
    -8.37077538483651e-001,
], dtype=np.float64)


SCALE = np.array([
    5.89399659952541e-001,
    4.79357856973238e-001,
    4.41834266003779e-001,
    6.39113086349240e-001,
    4.00221482285109e-001,
    3.12751514591303e-001,
    3.51081463499536e-001,
    4.16744089273568e-001,
    5.18100475607077e-001,
    5.21636279355979e-001,
    6.61611621630445e-001,
    4.84347902320302e-001,
    4.09071005081325e-001,
    8.26531686726199e-001,
], dtype=np.float64)


PMML_MIN = np.array([
    1.81587824661900e-001,
    1.65932026146051e+000,
   -7.34774717499520e-003,
    2.96272661487389e+000,
    1.80172663311694e-001,
   -2.01495973809767e-002,
    7.34259313690555e-002,
    2.34859751196645e-001,
    1.29573343153477e+000,
    9.03799264669801e-001,
    4.74650411605838e-001,
   -1.17001685401852e-003,
    1.06353115424661e+000,
    1.01275916208273e+000,
], dtype=np.float64)


PMML_MAX = np.array([
    1.87822945502952e+000,
    3.74544440744656e+000,
    2.25594434432487e+000,
    4.52739493626945e+000,
    2.67878916508064e+000,
    3.17727696442730e+000,
    2.92176770946274e+000,
    2.63441387983015e+000,
    3.22586097837453e+000,
    2.82084384070009e+000,
    1.98611116487321e+000,
    2.06346161509779e+000,
    3.50809453707839e+000,
    2.22263413246758e+000,
], dtype=np.float64)


# =====================================================================
# EXACT HIDDEN LAYER
#
# 7 tanh neurons
# =====================================================================

HIDDEN_BIAS = np.array([
    -1.53490257231569,
     0.131657023685626,
     1.43853764533907,
    -1.86441750440977,
    -1.39007362347401,
    -2.27245624254360,
    -1.41238964397195,
], dtype=np.float64)


HIDDEN_WEIGHTS = np.array([
    [
         3.86636425127573,
         1.22136962869416,
        -2.72359339884494,
        -1.77230446723615,
        -1.46216194289540,
         3.24558345763305,
         3.06996215451877,
        -5.77021402108199,
         2.25881864469434,
         3.71517728820130,
        -0.309473949823228,
         0.726212042761451,
         0.522190869640098,
        -2.39263599373798,
    ],
    [
         2.15896216805703,
        -1.57953776828174,
        -2.01711047435665,
        -0.750351056198362,
        -1.92743140136514,
         0.299847359614649,
        -2.22202910077016,
        -0.0476928013017713,
        -0.759498711144356,
        -1.05109239442980,
         0.351623715392210,
         3.54202422454860,
         1.89440838853097,
        -2.41189170155559,
    ],
    [
        -0.922813177953713,
        -3.34649069356616,
        -1.97874371696066,
         4.33203158757785,
         4.31521010368907,
        -5.27616181448894,
        -7.83331802110645,
        -3.45801845603552,
        -5.46339923062738,
        -6.66285394819879,
        -2.55681011469760,
        -3.50785278810738,
        -5.24223013070959,
        -4.24405894617469,
    ],
    [
         1.77158646977566,
         4.83593638728132,
         2.93633574574327,
        -5.73872600600682,
        -5.70284050186134,
         7.35116939670637,
        10.69544458105640,
         4.89961793180250,
         7.12785345464860,
         8.87746558581506,
         3.85582936805749,
         5.06143738013748,
         7.13807633592428,
         6.01517428721663,
    ],
    [
         1.72255469608114,
        -1.86165729356541,
        -1.90033969603888,
         2.37957191920508,
         1.19976490900057,
         1.50396988828879,
         0.568345390367603,
        -2.13131273172663,
         2.41361647447767,
         2.28470508155423,
         0.733410779420942,
         2.02146659723445,
        -0.773587812850011,
        -1.53745590427323,
    ],
    [
        -2.90443799319791,
        -3.68208822544325,
        -3.46789504316015,
         0.837052647324236,
         0.765436481536820,
        -4.99581134896554,
        -6.86730880201060,
        -4.18010466512886,
        -4.79366603735720,
        -5.54970156820046,
        -3.81043212715729,
        -4.52104578830334,
        -5.40512941600826,
        -4.46376429310638,
    ],
    [
        -0.818814960793102,
        -1.82559898853851,
         2.16129711051170,
         2.32670843074396,
        -1.03267361747515,
        -2.32475363308934,
         0.662257264300905,
         5.73805413129751,
        -2.24874818795921,
        -4.25178800936449,
        -0.467668737217193,
         0.596299435874968,
         1.44591620076925,
         2.00807879352191,
    ],
], dtype=np.float64)


# =====================================================================
# EXACT OUTPUT LAYER
#
# 2 tanh neurons
# =====================================================================

OUTPUT_BIAS = np.array([
    -4.14151909582033,
     4.23842029656941,
], dtype=np.float64)


OUTPUT_WEIGHTS = np.array([
    [
         9.74027777331552,
        -1.90077878045881,
       -12.65596685335290,
         1.64281899413842,
        -4.11740620285455,
         5.74336157346222,
         7.28595015138478,
    ],
    [
        -9.77050594103395,
         1.90978476815711,
        12.56065499298740,
        -1.70819333611325,
         4.11135939489197,
        -5.75632057594169,
        -7.25766983514833,
    ],
], dtype=np.float64)


# =====================================================================
# HELPERS
# =====================================================================

def normalize_id(value: object) -> str:

    if value is None:
        return ""

    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass

    s = str(value).strip().lower()

    # Numeric Excel values:
    #
    # 50217.0 -> 50217
    # 2.0     -> 2
    #
    if re.fullmatch(r"\d+\.0+", s):
        s = s.split(".", 1)[0]

    # E-1001 -> e1001
    s = re.sub(r"[^a-z0-9]", "", s)

    return s


def ensure_files_exist() -> None:

    required = [
        SUPPLEMENTARY_XLSX,
        METADATA_PATH,
        COHORT1_XLSX,
        COHORT2_XLSX,
    ]

    missing = [
        str(p)
        for p in required
        if not p.exists()
    ]

    if missing:

        print()
        print("=" * 100)
        print("MISSING REQUIRED FILES")
        print("=" * 100)

        for item in missing:
            print(f"  {item}")

        raise RuntimeError(
            "Required historical reproduction files are missing."
        )


def workbook_sample_ids(
    workbook: Path,
    sheet_name: str = "mirna_tpm",
) -> list[str]:

    header = pd.read_excel(
        workbook,
        sheet_name=sheet_name,
        header=2,
        nrows=0,
        engine="openpyxl",
    )

    raw_columns = list(header.columns)

    # The first four columns are:
    #
    # Name
    # Sequence
    # Name count
    # Identical miRNAs
    #
    # Everything after that is a sample column.

    sample_columns = raw_columns[4:]

    result = []

    for value in sample_columns:

        normalized = normalize_id(value)

        if normalized:
            result.append(normalized)

    return result


# =====================================================================
# RECOVER ORIGINAL 179-SAMPLE ORDER
# =====================================================================

def recover_historical_sample_order() -> tuple[list[str], list[str]]:

    print()
    print("=" * 100)
    print("RECOVERING EXACT HISTORICAL SAMPLE ORDER")
    print("=" * 100)

    # ---------------------------------------------------------------
    # Cohort 2 workbook contains 120 samples.
    #
    # Dataset 2 contains 179 expression columns:
    #
    #   numeric headers corresponding to Cohort 2 positions
    #   +
    #   60 explicit e#### samples
    #
    # Critically, Dataset 2 already omits numeric position 56.
    # Therefore the historical cohort is:
    #
    #   Cohort 2 positions 1..55
    #   Cohort 2 positions 57..120
    #   all 60 Dataset-2 e-samples
    #
    # E-1056 MUST NOT be removed.
    # ---------------------------------------------------------------

    cohort2 = workbook_sample_ids(COHORT2_XLSX)
    cohort1 = workbook_sample_ids(COHORT1_XLSX)

    print()
    print(f"Cohort 2 workbook samples: {len(cohort2)}")
    print(f"Cohort 1 workbook samples: {len(cohort1)}")

    if len(cohort2) != 120:
        raise RuntimeError(
            f"Expected 120 Cohort 2 samples, found {len(cohort2)}."
        )

    if len(cohort1) != 60:
        raise RuntimeError(
            f"Expected 60 Cohort 1 samples, found {len(cohort1)}."
        )

    # ---------------------------------------------------------------
    # Read Dataset 2 header directly.
    # ---------------------------------------------------------------

    raw = pd.read_excel(
        SUPPLEMENTARY_XLSX,
        sheet_name="Supplemental Dataset 2",
        header=None,
        engine="openpyxl",
    )

    header_row = 2

    raw_header = list(
        raw.iloc[
            header_row,
            1:180,
        ]
    )

    if len(raw_header) != 179:
        raise RuntimeError(
            f"Expected 179 Dataset 2 expression columns, "
            f"found {len(raw_header)}."
        )

    historical = []

    numeric_positions = []
    observed_esamples = []

    for column_number, value in enumerate(raw_header, start=1):

        # Numeric Dataset 2 headers represent Cohort 2 positions.
        try:
            numeric = float(value)

            if numeric.is_integer():
                position = int(numeric)

                if not (1 <= position <= 120):
                    raise RuntimeError(
                        f"Dataset 2 numeric header {value!r} at "
                        f"column {column_number} is outside 1..120."
                    )

                sample_id = cohort2[position - 1]

                historical.append(sample_id)
                numeric_positions.append(position)
                continue

        except (TypeError, ValueError):
            pass

        # Non-numeric headers should be explicit e#### IDs.
        sample_id = normalize_id(value)

        if not sample_id.startswith("e"):
            raise RuntimeError(
                f"Unexpected Dataset 2 sample header at "
                f"column {column_number}: {value!r}"
            )

        historical.append(sample_id)
        observed_esamples.append(sample_id)

    # ---------------------------------------------------------------
    # Validate exact historical structure.
    # ---------------------------------------------------------------

    expected_numeric_positions = [
        *range(1, 56),
        *range(57, 121),
    ]

    if numeric_positions != expected_numeric_positions:
        raise RuntimeError(
            "Dataset 2 numeric sample positions do not match the "
            "historical structure.\n"
            f"Observed: {numeric_positions}\n"
            f"Expected: {expected_numeric_positions}"
        )

    if len(observed_esamples) != 60:
        raise RuntimeError(
            f"Expected 60 e-samples in Dataset 2, "
            f"found {len(observed_esamples)}."
        )

    expected_esamples = set(cohort1)
    observed_esample_set = set(observed_esamples)

    if observed_esample_set != expected_esamples:
        missing = sorted(expected_esamples - observed_esample_set)
        unexpected = sorted(observed_esample_set - expected_esamples)

        raise RuntimeError(
            "Dataset 2 e-sample set does not match Cohort 1.\n"
            f"Missing: {missing}\n"
            f"Unexpected: {unexpected}"
        )

    if len(historical) != 179:
        raise RuntimeError(
            f"Expected 179 historical samples, found {len(historical)}."
        )

    if len(set(historical)) != 179:
        raise RuntimeError(
            "Historical sample mapping contains duplicate sample IDs."
        )

    if cohort2[55] in historical:
        raise RuntimeError(
            "Excluded numeric Cohort 2 position 56 "
            f"({cohort2[55]}) is present in historical data."
        )

    if "e1056" not in historical:
        raise RuntimeError(
            "E-1056 is unexpectedly absent. "
            "It should be retained in the historical dataset."
        )

    print()
    print("Exact historical sample mapping recovered from Dataset 2 header:")

    print()
    print("First 10:")
    for i, sid in enumerate(historical[:10], start=1):
        print(f"  {i:3d}: {sid}")

    print()
    print("Around excluded numeric position 56:")
    for i in range(50, 62):
        if i <= len(historical):
            print(f"  {i:3d}: {historical[i - 1]}")

    print()
    print("Last 10:")
    start = len(historical) - 9
    for i, sid in enumerate(
        historical[-10:],
        start=start,
    ):
        print(f"  {i:3d}: {sid}")

    print()
    print("Historical sample count:", len(historical))
    print("Numeric samples retained:", len(numeric_positions))
    print("E-samples retained:", len(observed_esamples))
    print("Excluded numeric position: 56")
    print("Excluded mapped sample:", cohort2[55])
    print("E-1056 retained: YES")

    return historical, cohort2




# =====================================================================
# LOAD DATASET 2
# =====================================================================

def load_dataset2(
    historical_ids: list[str],
) -> pd.DataFrame:

    print()
    print("=" * 100)
    print("READING eLIFE SUPPLEMENTAL DATASET 2")
    print("=" * 100)

    raw = pd.read_excel(
        SUPPLEMENTARY_XLSX,
        sheet_name="Supplemental Dataset 2",
        header=None,
        engine="openpyxl",
    )

    print()
    print(f"Raw shape: {raw.shape}")

    header_row = 2

    raw_header = list(
        raw.iloc[
            header_row,
            1:180,
        ]
    )

    if len(raw_header) != 179:
        raise RuntimeError(
            "Dataset 2 did not contain exactly 179 sample columns."
        )

    print()
    print("FIRST 10 RAW HEADER VALUES:")

    for i, value in enumerate(
        raw_header[:10],
        start=1,
    ):

        print(
            f"  column {i:3d}: {value!r}"
        )

    print()
    print("LAST 10 RAW HEADER VALUES:")

    start = 170

    for i, value in enumerate(
        raw_header[-10:],
        start=start,
    ):

        print(
            f"  column {i:3d}: {value!r}"
        )

    # ---------------------------------------------------------------
    # HEADER AUDIT
    #
    # Dataset 2 contains 179 expression columns. Its header is a
    # mixed representation:
    #
    #   early columns use positional sample numbers
    #   later columns contain explicit E-#### sample IDs
    #
    # We do NOT use these header values to identify patients.
    #
    # The historical patient identity/order was independently recovered
    # from the two processed cohort workbooks above:
    #
    #   positions 1..120   -> Cohort 2 workbook order
    #   positions 121..179 -> Cohort 1 workbook order
    #
    # Therefore the Dataset 2 columns are treated as POSITIONAL.
    # This avoids depending on pandas' interpretation of mixed Excel
    # headers while preserving the original column ordering.
    # ---------------------------------------------------------------

    if len(raw_header) != 179:
        raise RuntimeError(
            "Dataset 2 does not contain exactly 179 sample columns."
        )

    # Audit the first 120 headers without making them a prerequisite
    # for the reconstruction.

    positional_header_ok = True
    first_bad_position = None
    first_bad_value = None

    for expected, observed in enumerate(
        raw_header[:120],
        start=1,
    ):

        try:
            numeric = float(observed)
        except (TypeError, ValueError):
            positional_header_ok = False
            first_bad_position = expected
            first_bad_value = observed
            break

        if not math.isclose(
            numeric,
            float(expected),
            rel_tol=0.0,
            abs_tol=1e-9,
        ):
            positional_header_ok = False
            first_bad_position = expected
            first_bad_value = observed
            break

    print()

    if positional_header_ok:

        print(
            "Dataset 2 header audit: first 120 headers "
            "match positions 1..120."
        )

    else:

        print(
            "Dataset 2 header audit: mixed/non-positional header "
            "detected."
        )

        print(
            f"  first differing position: {first_bad_position}"
        )

        print(
            f"  observed value: {first_bad_value!r}"
        )

        print(
            "  This does NOT affect reconstruction because "
            "sample identity is assigned by column position."
        )

    print()
    print(
        "Using independently recovered historical sample order "
        "for all 179 Dataset 2 columns."
    )

    # ---------------------------------------------------------------
    # Extract exact PMML features.
    #
    # First column is the miRNA name.
    # Columns 2..180 are the 179 samples.
    # ---------------------------------------------------------------

    feature_series = raw.iloc[
        :,
        0,
    ].astype(str).str.strip()

    rows = {}

    for feature in FEATURES:

        matches = np.where(
            feature_series.values == feature
        )[0]

        if len(matches) != 1:

            raise RuntimeError(
                f"Expected exactly one row for "
                f"{feature}, found {len(matches)}."
            )

        rows[feature] = int(matches[0])

    print()
    print("Exact PMML feature rows:")

    for feature in FEATURES:

        print(
            f"  {feature}: row {rows[feature]}"
        )

    # ---------------------------------------------------------------
    # Pull numeric matrix.
    # ---------------------------------------------------------------

    values = []

    for feature in FEATURES:

        row_index = rows[feature]

        row_values = raw.iloc[
            row_index,
            1:180,
        ]

        numeric = pd.to_numeric(
            row_values,
            errors="coerce",
        ).to_numpy(
            dtype=np.float64,
        )

        if len(numeric) != 179:
            raise RuntimeError(
                f"{feature} has {len(numeric)} values "
                "instead of 179."
            )

        if not np.all(
            np.isfinite(numeric)
        ):

            bad = np.where(
                ~np.isfinite(numeric)
            )[0]

            raise RuntimeError(
                f"{feature} contains non-finite values "
                f"at Dataset 2 columns {bad[:10].tolist()}."
            )

        values.append(numeric)

    # Features x samples -> samples x features
    X = np.vstack(values).T

    dataset = pd.DataFrame(
        X,
        columns=FEATURES,
    )

    dataset.insert(
        0,
        "sample_id",
        historical_ids,
    )

    print()
    print(
        f"Historical Dataset 2 matrix: "
        f"{dataset.shape[0]} patients x "
        f"{len(FEATURES)} PMML features"
    )

    return dataset


# =====================================================================
# LOAD GEO LABELS
# =====================================================================

def load_labels(
    historical_ids: list[str],
) -> pd.DataFrame:

    print()
    print("=" * 100)
    print("LOADING GEO DIAGNOSIS LABELS")
    print("=" * 100)

    metadata = pd.read_csv(
        METADATA_PATH
    )

    required = {
        "title",
        "diagnosis",
    }

    missing = required - set(metadata.columns)

    if missing:
        raise RuntimeError(
            "GEO metadata is missing required columns: "
            f"{sorted(missing)}\n"
            f"Available columns: {metadata.columns.tolist()}"
        )

    metadata = metadata.copy()

    # ---------------------------------------------------------------
    # The GEO metadata stores the original biological sample ID in
    # the `title` field.
    #
    # Examples:
    #   50217  -> "50217"
    #   E-1001 -> "E-1001"
    #
    # Normalize those identifiers so that they match the historical
    # IDs recovered from the supplementary expression workbook.
    # ---------------------------------------------------------------

    metadata["sample_id"] = (
        metadata["title"]
        .astype(str)
        .str.strip()
        .str.lower()
        .str.replace("-", "", regex=False)
        .str.replace(r"\.0+$", "", regex=True)
    )

    # ---------------------------------------------------------------
    # Verify uniqueness.
    # ---------------------------------------------------------------

    duplicate_ids = (
        metadata.loc[
            metadata["sample_id"].duplicated(keep=False),
            "sample_id",
        ]
        .drop_duplicates()
        .tolist()
    )

    if duplicate_ids:
        raise RuntimeError(
            "GEO metadata contains duplicate normalized sample IDs:\n"
            + "\n".join(sorted(duplicate_ids))
        )

    # ---------------------------------------------------------------
    # Index by biological sample ID.
    # ---------------------------------------------------------------

    indexed = metadata.set_index(
        "sample_id"
    )

    missing_ids = [
        sid
        for sid in historical_ids
        if sid not in indexed.index
    ]

    if missing_ids:
        raise RuntimeError(
            "Historical samples missing GEO metadata:\n"
            + "\n".join(missing_ids)
        )

    # ---------------------------------------------------------------
    # Recover labels in EXACT historical expression order.
    # ---------------------------------------------------------------

    selected = indexed.loc[
        historical_ids,
        [
            "diagnosis",
            "stage",
            "histology",
            "title",
        ],
    ].copy()

    selected.insert(
        0,
        "sample_id",
        historical_ids,
    )

    # ---------------------------------------------------------------
    # Binary endpoint:
    #
    #   Cancer = 1
    #   Benign/Borderline/Control = 0
    # ---------------------------------------------------------------

    diagnosis = (
        selected["diagnosis"]
        .astype(str)
        .str.strip()
        .str.title()
    )

    valid_diagnoses = {
        "Cancer",
        "Benign",
        "Borderline",
        "Control",
    }

    unexpected = sorted(
        set(diagnosis) - valid_diagnoses
    )

    if unexpected:
        raise RuntimeError(
            "Unexpected diagnosis values:\n"
            + "\n".join(unexpected)
        )

    selected["diagnosis"] = diagnosis

    selected["true_label"] = (
        selected["diagnosis"] == "Cancer"
    ).astype(int)

    # ---------------------------------------------------------------
    # Audit counts.
    # ---------------------------------------------------------------

    print()
    print("Historical diagnosis values:")

    counts = (
        selected["diagnosis"]
        .value_counts()
        .sort_index()
    )

    for diagnosis_name, count in counts.items():
        print(
            f"  {diagnosis_name:<12} {count}"
        )

    print()
    print(
        f"Cancer:      "
        f"{int(selected['true_label'].sum())}"
    )

    print(
        f"Non-cancer:  "
        f"{int((selected['true_label'] == 0).sum())}"
    )

    return selected






# =====================================================================
# VERIFY PMML NORMALIZATION
# =====================================================================

def verify_normalization() -> None:

    print()
    print("=" * 100)
    print("VERIFYING PMML NORMALIZATION PARAMETERS")
    print("=" * 100)

    calculated_scale = (
        1.0
        / (
            PMML_MAX
            - PMML_MIN
        )
    )

    calculated_shift = (
        -PMML_MIN
        * calculated_scale
    )

    scale_ok = np.allclose(
        calculated_scale,
        SCALE,
        rtol=0.0,
        atol=1e-12,
    )

    shift_ok = np.allclose(
        calculated_shift,
        SHIFT,
        rtol=0.0,
        atol=1e-12,
    )

    print()
    print(
        f"Scale agrees with PMML endpoints: {scale_ok}"
    )

    print(
        f"Shift agrees with PMML endpoints: {shift_ok}"
    )

    if not scale_ok or not shift_ok:

        raise RuntimeError(
            "PMML normalization constants failed verification."
        )

    # Verify both endpoints map to approximately 0 and 1.
    lower = (
        PMML_MIN
        * SCALE
        + SHIFT
    )

    upper = (
        PMML_MAX
        * SCALE
        + SHIFT
    )

    if not np.allclose(
        lower,
        0.0,
        atol=1e-12,
    ):

        raise RuntimeError(
            "PMML minimum normalization does not map to zero."
        )

    if not np.allclose(
        upper,
        1.0,
        atol=1e-12,
    ):

        raise RuntimeError(
            "PMML maximum normalization does not map to one."
        )

    print(
        "Endpoint normalization check: PASS"
    )


# =====================================================================
# FORWARD PASS
# =====================================================================

def forward_pass(
    X_raw: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:

    if X_raw.ndim != 2:
        raise ValueError(
            "X_raw must be a 2D matrix."
        )

    if X_raw.shape[1] != 14:
        raise ValueError(
            f"Expected 14 features, got {X_raw.shape[1]}."
        )

    # ---------------------------------------------------------------
    # PMML input transformation.
    # ---------------------------------------------------------------

    X = (
        X_raw
        * SCALE
        + SHIFT
    )

    # ---------------------------------------------------------------
    # Hidden layer.
    #
    # h = tanh(XW^T + bias)
    # ---------------------------------------------------------------

    hidden_pre = (
        X
        @ HIDDEN_WEIGHTS.T
        + HIDDEN_BIAS
    )

    hidden = np.tanh(
        hidden_pre
    )

    # ---------------------------------------------------------------
    # Output layer.
    # ---------------------------------------------------------------

    output_pre = (
        hidden
        @ OUTPUT_WEIGHTS.T
        + OUTPUT_BIAS
    )

    output = np.tanh(
        output_pre
    )

    return (
        X,
        hidden_pre,
        hidden,
        output,
    )


# =====================================================================
# SAVE SAMPLE MAPPING
# =====================================================================

def save_mapping(
    dataset: pd.DataFrame,
    metadata: pd.DataFrame,
) -> None:

    mapping = dataset[
        ["sample_id"]
    ].copy()

    mapping["diagnosis"] = [
        metadata.loc[
            sid,
            "diagnosis",
        ]
        for sid in mapping["sample_id"]
    ]

    mapping["binary_label"] = [
        int(
            metadata.loc[
                sid,
                "true_label",
            ]
        )
        for sid in mapping["sample_id"]
    ]

    mapping["historical_position"] = np.arange(
        1,
        len(mapping) + 1,
    )

    mapping[
        "historical_dataset_section"
    ] = [
        "cohort2"
        if position <= 120
        else "cohort1"
        for position
        in mapping["historical_position"]
    ]

    mapping.to_csv(
        MAPPING_OUTPUT,
        index=False,
    )


# =====================================================================
# MAIN
# =====================================================================

def main() -> None:

    print()
    print("=" * 100)
    print("CLEAN HISTORICAL eLIFE GSE94533 MODEL REPRODUCTION")
    print("=" * 100)

    print()
    print(
        "This script does NOT train a model."
    )

    print()
    print(
        f"Supplementary workbook:\n"
        f"  {SUPPLEMENTARY_XLSX}"
    )

    print()
    print(
        f"Metadata:\n"
        f"  {METADATA_PATH}"
    )

    ensure_files_exist()

    REPORT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ---------------------------------------------------------------
    # 1. Recover historical order.
    # ---------------------------------------------------------------

    historical_ids, cohort2_ids = (
        recover_historical_sample_order()
    )

    # ---------------------------------------------------------------
    # 2. Load historical Dataset 2.
    # ---------------------------------------------------------------

    dataset = load_dataset2(
        historical_ids
    )

    # ---------------------------------------------------------------
    # 3. Load labels.
    # ---------------------------------------------------------------

    metadata = load_labels(
        historical_ids
    )

    # ---------------------------------------------------------------
    # 4. Verify normalization.
    # ---------------------------------------------------------------

    verify_normalization()

    # ---------------------------------------------------------------
    # 5. Verify dimensions.
    # ---------------------------------------------------------------

    X_raw = dataset[
        FEATURES
    ].to_numpy(
        dtype=np.float64,
    )

    y = np.array(
        [
            int(
                metadata.loc[
                    sid,
                    "true_label",
                ]
            )
            for sid
            in dataset["sample_id"]
        ],
        dtype=np.int64,
    )

    if X_raw.shape != (
        179,
        14,
    ):

        raise RuntimeError(
            f"Expected X shape (179, 14), got {X_raw.shape}."
        )

    if y.shape != (179,):

        raise RuntimeError(
            f"Expected y shape (179,), got {y.shape}."
        )

    if not np.all(
        np.isfinite(X_raw)
    ):

        raise RuntimeError(
            "Historical expression matrix contains "
            "non-finite values."
        )

    # ---------------------------------------------------------------
    # 6. Run exact forward pass.
    # ---------------------------------------------------------------

    (
        X_normalized,
        hidden_pre,
        hidden,
        output,
    ) = forward_pass(
        X_raw
    )

    cancer_activation = (
        output[:, 0]
    )

    control_activation = (
        output[:, 1]
    )

    predicted_class = (
        cancer_activation
        > control_activation
    ).astype(int)

    # ---------------------------------------------------------------
    # 7. Metrics.
    # ---------------------------------------------------------------

    # Avoid requiring sklearn for the core reproduction.
    try:

        from sklearn.metrics import (
            accuracy_score,
            confusion_matrix,
            roc_auc_score,
            sensitivity_score,
            specificity_score,
        )

    except ImportError:

        accuracy_score = None
        confusion_matrix = None
        roc_auc_score = None
        sensitivity_score = None
        specificity_score = None

    if roc_auc_score is not None:

        auc = float(
            roc_auc_score(
                y,
                cancer_activation,
            )
        )

        accuracy = float(
            accuracy_score(
                y,
                predicted_class,
            )
        )

        cm = confusion_matrix(
            y,
            predicted_class,
        )

        tn, fp, fn, tp = (
            int(cm[0, 0]),
            int(cm[0, 1]),
            int(cm[1, 0]),
            int(cm[1, 1]),
        )

    else:

        auc = None

        tn = int(
            np.sum(
                (y == 0)
                & (predicted_class == 0)
            )
        )

        fp = int(
            np.sum(
                (y == 0)
                & (predicted_class == 1)
            )
        )

        fn = int(
            np.sum(
                (y == 1)
                & (predicted_class == 0)
            )
        )

        tp = int(
            np.sum(
                (y == 1)
                & (predicted_class == 1)
            )
        )

        accuracy = (
            (tp + tn)
            / len(y)
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

    ties = int(
        np.sum(
            np.isclose(
                cancer_activation,
                control_activation,
                rtol=0.0,
                atol=1e-12,
            )
        )
    )

    # ---------------------------------------------------------------
    # 8. Prediction table.
    # ---------------------------------------------------------------

    predictions = dataset[
        ["sample_id"]
    ].copy()

    predictions["diagnosis"] = [
        metadata.loc[
            sid,
            "diagnosis",
        ]
        for sid in predictions[
            "sample_id"
        ]
    ]

    predictions["true_label"] = y

    predictions[
        "cancer_activation"
    ] = cancer_activation

    predictions[
        "controls_or_borderline_activation"
    ] = control_activation

    predictions[
        "predicted_label"
    ] = predicted_class

    predictions[
        "predicted_class"
    ] = np.where(
        predicted_class == 1,
        "Ovarian Cancer",
        "Controls or Borderline",
    )

    predictions.to_csv(
        PREDICTIONS_OUTPUT,
        index=False,
    )

    # ---------------------------------------------------------------
    # 9. Forward-pass trace.
    # ---------------------------------------------------------------

    trace = dataset[
        ["sample_id"]
        + FEATURES
    ].copy()

    for i in range(14):

        trace[
            f"normalized_{i+1:02d}"
        ] = X_normalized[:, i]

    for i in range(7):

        trace[
            f"hidden_{i+14}"
        ] = hidden[:, i]

    trace[
        "output_21_cancer"
    ] = output[:, 0]

    trace[
        "output_22_controls_or_borderline"
    ] = output[:, 1]

    trace.to_csv(
        TRACE_OUTPUT,
        index=False,
    )

    # ---------------------------------------------------------------
    # 10. Mapping.
    # ---------------------------------------------------------------

    save_mapping(
        dataset,
        metadata,
    )

    # ---------------------------------------------------------------
    # 11. Results.
    # ---------------------------------------------------------------

    results = {
        "dataset": {
            "source": "eLife Supplemental Dataset 2",
            "representation": (
                "batch adjusted, log10-transformed "
                "miRNA expression"
            ),
            "patients": int(len(y)),
            "features": 14,
            "feature_names": FEATURES,
        },

        "sample_order": {
            "cohort2_count": 120,
            "cohort1_count_after_e1056_exclusion": 59,
            "total": 179,
            "excluded_sample": "e1056",
        },

        "labels": {
            "cancer": int(np.sum(y == 1)),
            "non_cancer": int(np.sum(y == 0)),
        },

        "model": {
            "architecture": "14-7-2",
            "hidden_activation": "tanh",
            "output_activation": "tanh",
            "output_neuron_21": "Ovarian Cancer",
            "output_neuron_22": (
                "Controls or Borderline"
            ),
            "probability_transformation_used": False,
        },

        "normalization": {
            "formula": (
                "normalized = raw * scale + shift"
            ),
            "verified_against_pmml_endpoints": True,
        },

        "metrics_on_historical_179_patients": {
            "auc": auc,
            "accuracy": accuracy,
            "sensitivity": sensitivity,
            "specificity": specificity,
            "tp": tp,
            "tn": tn,
            "fp": fp,
            "fn": fn,
            "ties": ties,
        },

        "files": {
            "predictions": str(
                PREDICTIONS_OUTPUT
            ),
            "forward_trace": str(
                TRACE_OUTPUT
            ),
            "sample_mapping": str(
                MAPPING_OUTPUT
            ),
        },

        "training_performed": False,
    }

    RESULTS_OUTPUT.write_text(
        json.dumps(
            results,
            indent=2,
        )
    )

    # =================================================================
    # FINAL REPORT
    # =================================================================

    print()
    print("=" * 100)
    print("HISTORICAL REPRODUCTION COMPLETE")
    print("=" * 100)

    print()
    print("Dataset:")
    print(f"  Patients:  {len(y)}")
    print("  Features:  14")
    print("  E-1056 excluded: YES")

    print()
    print("Labels:")
    print(f"  Cancer:      {int(np.sum(y == 1))}")
    print(f"  Non-cancer:  {int(np.sum(y == 0))}")

    print()
    print("Architecture:")
    print("  14 -> 7 tanh -> 2 tanh")

    print()
    print("Outputs:")
    print("  Neuron 21 = Ovarian Cancer")
    print("  Neuron 22 = Controls or Borderline")

    print()
    print("=" * 100)
    print("RESULTS")
    print("=" * 100)

    print()
    if auc is not None:
        print(
            f"ROC AUC:       {auc:.6f}"
        )
    else:
        print(
            "ROC AUC:       unavailable"
        )

    print(
        f"Accuracy:      {accuracy:.6f}"
    )

    print(
        f"Sensitivity:   {sensitivity:.6f}"
    )

    print(
        f"Specificity:   {specificity:.6f}"
    )

    print()
    print("Confusion matrix:")
    print(
        f"  TP={tp} FN={fn}"
    )
    print(
        f"  FP={fp} TN={tn}"
    )

    print()
    print(
        f"Prediction ties: {ties}"
    )

    print()
    print("=" * 100)
    print("FILES CREATED")
    print("=" * 100)

    print()
    print(
        "Predictions:"
    )
    print(
        f"  {PREDICTIONS_OUTPUT}"
    )

    print()
    print(
        "Forward trace:"
    )
    print(
        f"  {TRACE_OUTPUT}"
    )

    print()
    print(
        "Sample mapping:"
    )
    print(
        f"  {MAPPING_OUTPUT}"
    )

    print()
    print(
        "Results JSON:"
    )
    print(
        f"  {RESULTS_OUTPUT}"
    )

    print()
    print("=" * 100)
    print("NO MODEL TRAINING WAS PERFORMED")
    print("=" * 100)


if __name__ == "__main__":
    main()
