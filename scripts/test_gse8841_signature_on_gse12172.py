#!/usr/bin/env python3

import gzip
import csv
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.preprocessing import StandardScaler


ROOT = Path.home() / "ovarian-ai"

DISCOVERY = (
    ROOT / "data/processed/"
    "gse8841_borderline_expression.csv.gz"
)

TEST = (
    ROOT / "data/external/gse12172/"
    "GSE12172_series_matrix.txt.gz"
)

GPL5689 = (
    ROOT / "data/external/gse_platforms/"
    "GPL5689.annot.gz"
)

GPL570 = (
    ROOT / "data/external/gse_platforms/"
    "GPL570.txt"
)

SIGNATURE_FILE = (
    ROOT / "data/processed/"
    "gse8841_frozen_borderline_genes.txt"
)

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
    "IFI30",
]


def read_platform_gpl5689(path):
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

    table = pd.DataFrame(
        [r.split("\t") for r in rows[1:]],
        columns=header,
    )

    table["ID"] = (
        table["ID"]
        .astype(str)
        .str.strip()
    )

    table["Gene symbol"] = (
        table["Gene symbol"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.upper()
    )

    return table


def read_platform_gpl570(path):
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

    table = pd.DataFrame(
        [r.split("\t") for r in rows[1:]],
        columns=header,
    )

    table["ID"] = (
        table["ID"]
        .astype(str)
        .str.strip()
    )

    table["Gene Symbol"] = (
        table["Gene Symbol"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.upper()
    )

    return table


def read_gse_matrix_metadata(path):
    metadata = {}

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

            if line.startswith("!Sample_"):

                parts = next(
                    csv.reader(
                        [line.rstrip("\n")],
                        delimiter="\t",
                        quotechar='"',
                    )
                )

                metadata.setdefault(
                    parts[0],
                    [],
                ).append(parts[1:])

    return metadata


def read_gse_matrix(path):
    matrix_lines = []

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
                matrix_lines.append(line)

    return pd.read_csv(
        StringIO(
            "".join(matrix_lines)
        ),
        sep="\t",
        quotechar='"',
    )


print("=" * 80)
print("FROZEN GSE8841 SIGNATURE -> GSE12172 EXTERNAL VALIDATION")
print("=" * 80)

# ------------------------------------------------------------------
# Load frozen signature
# ------------------------------------------------------------------

if SIGNATURE_FILE.exists():

    frozen = [
        x.strip().upper()
        for x in SIGNATURE_FILE.read_text().splitlines()
        if x.strip()
    ]

else:
    frozen = GENES

GENES = [
    g for g in GENES
    if g in frozen
]

if not GENES:
    raise RuntimeError(
        "Frozen signature file contains none of the expected genes."
    )

print()
print("Frozen genes:", len(GENES))

for g in GENES:
    print(" ", g)

# ------------------------------------------------------------------
# GSE8841 discovery data
# ------------------------------------------------------------------

disc = pd.read_csv(DISCOVERY)

disc = disc[
    disc["histotype"]
    .astype(str)
    .str.lower()
    .isin(["serous", "mucinous"])
].copy()

disc_meta = [
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
    c for c in disc.columns
    if c not in disc_meta
]

# ------------------------------------------------------------------
# GSE8841 probe -> gene
# ------------------------------------------------------------------

ann8841 = read_platform_gpl5689(GPL5689)

probe_to_gene_8841 = dict(
    zip(
        ann8841["ID"],
        ann8841["Gene symbol"],
    )
)

# ------------------------------------------------------------------
# Build gene-level discovery matrix
# Median across probes mapping to each gene.
# ------------------------------------------------------------------

disc_gene = {}

print()
print("=" * 80)
print("GSE8841 GENE MAPPING")
print("=" * 80)

for gene in GENES:

    probes = [
        p for p in disc_probes
        if probe_to_gene_8841.get(
            str(p),
            "",
        ) == gene
    ]

    print(
        f"{gene:15} probes={len(probes)}"
    )

    if not probes:
        raise RuntimeError(
            f"No GSE8841 probe mapped to {gene}"
        )

    values = (
        disc[probes]
        .apply(
            pd.to_numeric,
            errors="coerce",
        )
    )

    disc_gene[gene] = values.median(axis=1)

disc_X = pd.DataFrame(disc_gene)

if disc_X.isna().any().any():
    raise RuntimeError(
        "Missing values in GSE8841 gene signature."
    )

y_disc = (
    disc["label"]
    == "stage_I_cancer"
).astype(int)

print()
print(
    "Discovery patients:",
    len(disc_X)
)

print(
    "Stage-I cancer:",
    int(y_disc.sum())
)

print(
    "Borderline:",
    int((y_disc == 0).sum())
)

# ------------------------------------------------------------------
# Load GSE12172 metadata and expression matrix
# ------------------------------------------------------------------

meta12172 = read_gse_matrix_metadata(TEST)

expr12172 = read_gse_matrix(TEST)

probe_col = expr12172.columns[0]

expr12172[probe_col] = (
    expr12172[probe_col]
    .astype(str)
    .str.strip()
)

sample_accessions = (
    meta12172["!Sample_geo_accession"][0]
)

# ------------------------------------------------------------------
# Construct GSE12172 labels.
# LMP = borderline
# Malignant = cancer
# ------------------------------------------------------------------

char_rows = meta12172[
    "!Sample_characteristics_ch1"
]

labels = []

for i in range(
    len(sample_accessions)
):

    text = " | ".join(
        row[i]
        for row in char_rows
        if i < len(row)
    ).lower()

    if "lmp" in text:

        labels.append(
            "borderline"
        )

    elif "malignant" in text:

        labels.append(
            "malignant"
        )

    else:

        labels.append(
            "other"
        )

labels = np.array(labels)

keep = np.isin(
    labels,
    [
        "borderline",
        "malignant",
    ],
)

test_accessions = np.array(
    sample_accessions
)[keep]

test_labels = labels[keep]

print()
print("=" * 80)
print("GSE12172 TEST COHORT")
print("=" * 80)

print(
    "Total:",
    len(labels)
)

print(
    "Borderline/LMP:",
    int(
        (test_labels == "borderline").sum()
    )
)

print(
    "Malignant:",
    int(
        (test_labels == "malignant").sum()
    )
)

if (
    (test_labels == "borderline").sum()
    != 30
):
    raise RuntimeError(
        "Expected 30 LMP/borderline samples."
    )

if (
    (test_labels == "malignant").sum()
    != 60
):
    raise RuntimeError(
        "Expected 60 malignant samples."
    )

# ------------------------------------------------------------------
# GPL570 probe -> gene
# ------------------------------------------------------------------

ann12172 = read_platform_gpl570(GPL570)

probe_to_gene_12172 = dict(
    zip(
        ann12172["ID"],
        ann12172["Gene Symbol"],
    )
)

# ------------------------------------------------------------------
# Build the EXACT SAME 13 gene features.
# NO gene selection here.
# ------------------------------------------------------------------

test_gene = {}

print()
print("=" * 80)
print("GSE12172 GENE MAPPING")
print("=" * 80)

for gene in GENES:

    probes = [
        p for p in expr12172[probe_col]
        if probe_to_gene_12172.get(
            p,
            "",
        ) == gene
    ]

    print(
        f"{gene:15} probes={len(probes)}"
    )

    if not probes:
        raise RuntimeError(
            f"No GSE12172 probe mapped to {gene}"
        )

    sub = expr12172[
        expr12172[probe_col].isin(probes)
    ]

    values = (
        sub[test_accessions]
        .apply(
            pd.to_numeric,
            errors="coerce",
        )
        .median(axis=0)
    )

    test_gene[gene] = values

test_X = pd.DataFrame(
    test_gene,
    index=test_accessions,
)

if test_X.isna().any().any():
    raise RuntimeError(
        "Missing values in GSE12172 signature."
    )

# ------------------------------------------------------------------
# Standardize using DISCOVERY cohort ONLY.
# ------------------------------------------------------------------

scaler = StandardScaler()

X_disc = scaler.fit_transform(
    disc_X[GENES]
)

X_test = scaler.transform(
    test_X[GENES]
)

# ------------------------------------------------------------------
# Train final classifier on ALL GSE8841 discovery patients.
#
# The gene list was frozen BEFORE this external test.
# ------------------------------------------------------------------

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

prob = clf.predict_proba(
    X_test
)[:, 1]

pred = (
    prob >= 0.5
).astype(int)

# 1 = malignant
y_test = (
    test_labels == "malignant"
).astype(int)

# ------------------------------------------------------------------
# Metrics
# ------------------------------------------------------------------

auc = roc_auc_score(
    y_test,
    prob,
)

accuracy = accuracy_score(
    y_test,
    pred,
)

balanced_accuracy = (
    balanced_accuracy_score(
        y_test,
        pred,
    )
)

tn, fp, fn, tp = confusion_matrix(
    y_test,
    pred,
).ravel()

sensitivity = (
    tp / (tp + fn)
    if tp + fn
    else np.nan
)

specificity = (
    tn / (tn + fp)
    if tn + fp
    else np.nan
)

print()
print("=" * 80)
print("GSE12172 EXTERNAL VALIDATION RESULT")
print("=" * 80)

print()
print(
    f"ROC AUC:           {auc:.6f}"
)

print(
    f"Accuracy:           {accuracy:.6f}"
)

print(
    f"Balanced accuracy: {balanced_accuracy:.6f}"
)

print(
    f"Sensitivity:       {sensitivity:.6f}"
)

print(
    f"Specificity:       {specificity:.6f}"
)

print()
print("Confusion matrix:")
print(
    f"  TN={tn} FP={fp}"
)

print(
    f"  FN={fn} TP={tp}"
)

print()
print("=" * 80)
print("INTERPRETATION")
print("=" * 80)

if auc < 0.60:
    print(
        "The frozen signature shows weak external discrimination."
    )
elif auc < 0.70:
    print(
        "The frozen signature shows modest external discrimination."
    )
elif auc < 0.80:
    print(
        "The frozen signature shows promising external discrimination."
    )
else:
    print(
        "The frozen signature shows strong exploratory external discrimination."
    )

print()
print("=" * 80)
print("IMPORTANT")
print("=" * 80)

print(
    "GSE8841 = discovery cohort."
)

print(
    "GSE12172 = completely separate external test cohort."
)

print(
    "No gene selection was performed using GSE12172."
)

print(
    "This tests a tissue-expression signature."
)

print(
    "It does NOT validate the serum-miRNA diagnostic model."
)

print()
print("=" * 80)
print("DONE")
print("=" * 80)
