#!/usr/bin/env python3

import gzip
import csv
from io import StringIO

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler


ROOT = "data"

DISCOVERY = (
    f"{ROOT}/processed/"
    "gse8841_borderline_expression.csv.gz"
)

TEST = (
    f"{ROOT}/external/gse12172/"
    "GSE12172_series_matrix.txt.gz"
)

GPL5689 = (
    f"{ROOT}/external/gse_platforms/"
    "GPL5689.annot.gz"
)

GPL570 = (
    f"{ROOT}/external/gse_platforms/"
    "GPL570.txt"
)

# These are the 12 genes that successfully transferred
# from GPL5689 to GPL570.
GENES = [
    "VPS72",
    "ATP6V0B",
    "AHCY",
    "ITPKB",
    "SERPINF1",
    "NBR1",
    "PHKB",
    "POLR2F",
    "RSU1",
    "SRPX2",
    "CNN2",
    "HDGF",
]


# ================================================================
# PLATFORM READERS
# ================================================================

def load_gpl5689(path):

    rows = []

    with gzip.open(
        path,
        "rt",
        encoding="utf-8",
        errors="replace",
    ) as f:

        inside = False

        for line in f:

            line = line.rstrip("\n")

            if line.startswith(
                "!platform_table_begin"
            ):
                inside = True
                continue

            if line.startswith(
                "!platform_table_end"
            ):
                break

            if inside and line:
                rows.append(line)

    header = rows[0].split("\t")

    tab = pd.DataFrame(
        [r.split("\t") for r in rows[1:]],
        columns=header,
    )

    return dict(
        zip(
            tab["ID"].astype(str).str.strip(),
            tab["Gene symbol"]
            .fillna("")
            .astype(str)
            .str.strip()
            .str.upper(),
        )
    )


def load_gpl570(path):

    rows = []

    with open(
        path,
        "r",
        encoding="utf-8",
        errors="replace",
    ) as f:

        inside = False

        for line in f:

            line = line.rstrip("\n")

            if line.startswith(
                "!platform_table_begin"
            ):
                inside = True
                continue

            if line.startswith(
                "!platform_table_end"
            ):
                break

            if inside and line:
                rows.append(line)

    header = rows[0].split("\t")

    tab = pd.DataFrame(
        [r.split("\t") for r in rows[1:]],
        columns=header,
    )

    return dict(
        zip(
            tab["ID"].astype(str).str.strip(),
            tab["Gene Symbol"]
            .fillna("")
            .astype(str)
            .str.strip()
            .str.upper(),
        )
    )


# ================================================================
# GEO READERS
# ================================================================

def load_geo_metadata(path):

    meta = {}

    with gzip.open(
        path,
        "rt",
        encoding="utf-8",
        errors="replace",
    ) as f:

        for line in f:

            if line.startswith(
                "!series_matrix_table_begin"
            ):
                break

            if line.startswith(
                "!Sample_geo_accession"
            ):

                parts = next(
                    csv.reader(
                        [line.rstrip("\n")],
                        delimiter="\t",
                        quotechar='"',
                    )
                )

                meta["accessions"] = parts[1:]

            elif line.startswith(
                "!Sample_characteristics_ch1"
            ):

                parts = next(
                    csv.reader(
                        [line.rstrip("\n")],
                        delimiter="\t",
                        quotechar='"',
                    )
                )

                meta.setdefault(
                    "chars",
                    [],
                ).append(parts[1:])

    return meta


def load_geo_expression(path):

    lines = []

    with gzip.open(
        path,
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
                lines.append(line)

    return pd.read_csv(
        StringIO("".join(lines)),
        sep="\t",
        quotechar='"',
        low_memory=False,
    )


# ================================================================
# START
# ================================================================

print("=" * 80)
print("FROZEN SIGNATURE: BORDERLINE vs STAGE-I CANCER")
print("=" * 80)


# ================================================================
# 1. GSE8841 DISCOVERY COHORT
# ================================================================

disc = pd.read_csv(
    DISCOVERY,
    low_memory=False,
)

# Same strict cohort used to discover the signature.
disc = disc[
    disc["histotype"]
    .astype(str)
    .str.lower()
    .isin([
        "serous",
        "mucinous",
    ])
].copy()

metadata_cols = [
    "sample_id",
    "title",
    "label",
    "figo_stage",
    "figo_substage",
    "histotype",
    "grade",
    "age",
    "overall_survival_days",
    "progression_free_survival_days",
    "relapsed",
    "status",
]

disc_probes = [
    c
    for c in disc.columns
    if c not in metadata_cols
]

probe_to_gene_8841 = load_gpl5689(
    GPL5689
)

disc_gene = {}

for gene in GENES:

    matching_probes = [
        p
        for p in disc_probes
        if probe_to_gene_8841.get(
            str(p),
            "",
        ) == gene
    ]

    if not matching_probes:

        raise RuntimeError(
            f"No GSE8841 probe found for {gene}"
        )

    disc_gene[gene] = (
        disc[matching_probes]
        .apply(
            pd.to_numeric,
            errors="coerce",
        )
        .median(axis=1)
    )

disc_X = pd.DataFrame(
    disc_gene
)

if disc_X.isna().any().any():

    raise RuntimeError(
        "GSE8841 frozen signature contains missing values."
    )

y_disc = (
    disc["label"]
    == "stage_I_cancer"
).astype(int).to_numpy()

print()
print("DISCOVERY COHORT")
print("Patients:", len(disc_X))
print(
    "Stage-I cancer:",
    int(y_disc.sum())
)
print(
    "Borderline:",
    int((y_disc == 0).sum())
)


# ================================================================
# 2. GSE12172 METADATA
# ================================================================

meta = load_geo_metadata(
    TEST
)

accessions = meta["accessions"]
char_rows = meta["chars"]

print()
print("=" * 80)
print("GSE12172 METADATA")
print("=" * 80)

print(
    "Samples:",
    len(accessions)
)

print(
    "Characteristic rows:",
    len(char_rows)
)


# ================================================================
# 3. DIAGNOSIS ROW
# ================================================================

diagnosis_row = None

for row_number, row in enumerate(
    char_rows,
    start=1,
):

    normalized = {
        str(x)
        .strip()
        .lower()
        for x in row
    }

    if (
        "lmp" in normalized
        and "malignant" in normalized
    ):

        diagnosis_row = row

        print(
            "Diagnosis row:",
            row_number,
        )

        break

if diagnosis_row is None:

    raise RuntimeError(
        "Could not locate LMP/Malignant diagnosis row."
    )


labels = np.array([
    (
        "borderline"
        if str(x).strip().lower() == "lmp"
        else
        "malignant"
        if str(x).strip().lower() == "malignant"
        else
        "other"
    )
    for x in diagnosis_row
])


print()
print("DIAGNOSIS COUNTS")

for value, count in zip(
    *np.unique(
        labels,
        return_counts=True,
    )
):

    print(
        f"  {value}: {count}"
    )


# ================================================================
# 4. FIGO STAGE ROW
# ================================================================
#
# Your previous audit established that the FIGO stage is the
# 4th characteristic row (index 3).
#
# Values include:
#   Stage IA
#   Stage IB
#   Stage IC
#   Stage IIA
#   ...
#
# We use ONLY the stage information to define the test subset.
# Stage itself is NOT used as a model feature.
# ================================================================

if len(char_rows) < 4:

    raise RuntimeError(
        "GSE12172 has fewer than 4 characteristic rows."
    )

stage_row = char_rows[3]

print()
print("Using FIGO stage row: 4")


def normalize_stage(value):

    value = (
        str(value)
        .strip()
        .upper()
    )

    if value.startswith("STAGE IA"):
        return "IA"

    if value.startswith("STAGE IB"):
        return "IB"

    if value.startswith("STAGE IC"):
        return "IC"

    if value.startswith("STAGE IIA"):
        return "IIA"

    if value.startswith("STAGE IIB"):
        return "IIB"

    if value.startswith("STAGE IIC"):
        return "IIC"

    if value.startswith("STAGE IIIA"):
        return "IIIA"

    if value.startswith("STAGE IIIB"):
        return "IIIB"

    if value.startswith("STAGE IIIC"):
        return "IIIC"

    if value.startswith("STAGE IV"):
        return "IV"

    return ""


stages = np.array([
    normalize_stage(x)
    for x in stage_row
])


# ================================================================
# 5. DEFINE THE EXTERNAL STAGE-I TEST
# ================================================================

is_borderline = (
    labels == "borderline"
)

is_stage1_cancer = (
    (labels == "malignant")
    & np.isin(
        stages,
        [
            "IA",
            "IB",
            "IC",
        ],
    )
)

keep = (
    is_borderline
    | is_stage1_cancer
)

test_ids = (
    np.array(accessions)[keep]
)

test_labels = (
    labels[keep]
)

test_stages = (
    stages[keep]
)


print()
print("=" * 80)
print("GSE12172 STAGE-I TEST COHORT")
print("=" * 80)

print(
    "Borderline/LMP:",
    int(
        (test_labels == "borderline").sum()
    )
)

print(
    "Stage-I malignant:",
    int(
        (test_labels == "malignant").sum()
    )
)

print()
print("Stage-I malignant breakdown:")

malignant_stages = test_stages[
    test_labels == "malignant"
]

for stage, count in zip(
    *np.unique(
        malignant_stages,
        return_counts=True,
    )
):

    print(
        f"  {stage}: {count}"
    )

if int(
    (test_labels == "borderline").sum()
) != 30:

    raise RuntimeError(
        "Expected 30 borderline/LMP samples."
    )

if int(
    (test_labels == "malignant").sum()
) != 16:

    raise RuntimeError(
        "Expected 16 Stage-I malignant samples."
    )


# ================================================================
# 6. GSE12172 EXPRESSION
# ================================================================

expr = load_geo_expression(
    TEST
)

probe_col = expr.columns[0]

expr[probe_col] = (
    expr[probe_col]
    .astype(str)
    .str.strip()
)

probe_to_gene_12172 = load_gpl570(
    GPL570
)


# ================================================================
# 7. MAP THE SAME FROZEN GENES
# ================================================================

test_gene = {}

for gene in GENES:

    matching_probes = [
        p
        for p in expr[probe_col]
        if probe_to_gene_12172.get(
            p,
            "",
        ) == gene
    ]

    print(
        f"{gene:15} probes={len(matching_probes)}"
    )

    if not matching_probes:

        raise RuntimeError(
            f"Frozen gene {gene} "
            "cannot be measured on GPL570."
        )

    sub = expr[
        expr[probe_col]
        .isin(matching_probes)
    ]

    test_gene[gene] = (
        sub[test_ids]
        .apply(
            pd.to_numeric,
            errors="coerce",
        )
        .median(axis=0)
        .to_numpy(
            dtype=float
        )
    )


test_X = pd.DataFrame(
    test_gene,
    index=test_ids,
)

if test_X.isna().any().any():

    raise RuntimeError(
        "Missing external gene-expression values."
    )


# ================================================================
# 8. FINAL TRAINING ON GSE8841 ONLY
# ================================================================

scaler = StandardScaler()

X_disc = scaler.fit_transform(
    disc_X[GENES]
)

X_test = scaler.transform(
    test_X[GENES]
)

clf = LogisticRegression(
    C=0.1,
    class_weight="balanced",
    max_iter=5000,
    solver="liblinear",
)

clf.fit(
    X_disc,
    y_disc,
)


# ================================================================
# 9. EXTERNAL TEST
# ================================================================

prob = clf.predict_proba(
    X_test
)[:, 1]

# 1 = Stage-I malignant
y_test = (
    test_labels
    == "malignant"
).astype(int)


auc = roc_auc_score(
    y_test,
    prob,
)


# ================================================================
# RESULT
# ================================================================

print()
print("=" * 80)
print("FINAL STAGE-I EXTERNAL RESULT")
print("=" * 80)

print()
print(
    "Patients:",
    len(y_test)
)

print(
    "Borderline/LMP:",
    int((y_test == 0).sum())
)

print(
    "Stage-I malignant:",
    int((y_test == 1).sum())
)

print()
print(
    f"ROC AUC: {auc:.6f}"
)

print()
print("=" * 80)
print("IMPORTANT")
print("=" * 80)

print(
    "The 12 genes were frozen before GSE12172 testing."
)

print(
    "No feature selection was performed on GSE12172."
)

print(
    "No threshold was tuned on GSE12172."
)

print(
    "This is tissue-expression external validation."
)

print(
    "It is NOT validation of the serum-miRNA model."
)

print()
print("=" * 80)
print("DONE")
print("=" * 80)
