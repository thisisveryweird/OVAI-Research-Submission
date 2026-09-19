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

DISCOVERY = ROOT / "data/processed/gse8841_borderline_expression.csv.gz"
TEST = ROOT / "data/external/gse12172/GSE12172_series_matrix.txt.gz"
GPL5689 = ROOT / "data/external/gse_platforms/GPL5689.annot.gz"
GPL570 = ROOT / "data/external/gse_platforms/GPL570.txt"

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


def load_gpl5689(path):
    rows = []

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        inside = False

        for line in f:
            line = line.rstrip("\n")

            if line.startswith("!platform_table_begin"):
                inside = True
                continue

            if line.startswith("!platform_table_end"):
                break

            if inside and line:
                rows.append(line)

    header = rows[0].split("\t")

    df = pd.DataFrame(
        [r.split("\t") for r in rows[1:]],
        columns=header,
    )

    df["ID"] = df["ID"].astype(str).str.strip()
    df["Gene symbol"] = (
        df["Gene symbol"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.upper()
    )

    return dict(
        zip(df["ID"], df["Gene symbol"])
    )


def load_gpl570(path):
    rows = []

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        inside = False

        for line in f:
            line = line.rstrip("\n")

            if line.startswith("!platform_table_begin"):
                inside = True
                continue

            if line.startswith("!platform_table_end"):
                break

            if inside and line:
                rows.append(line)

    header = rows[0].split("\t")

    df = pd.DataFrame(
        [r.split("\t") for r in rows[1:]],
        columns=header,
    )

    df["ID"] = df["ID"].astype(str).str.strip()
    df["Gene Symbol"] = (
        df["Gene Symbol"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.upper()
    )

    return dict(
        zip(df["ID"], df["Gene Symbol"])
    )


def load_geo_metadata(path):
    metadata = {}

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:

            if line.startswith("!series_matrix_table_begin"):
                break

            if line.startswith("!Sample_geo_accession"):
                parts = next(
                    csv.reader(
                        [line.rstrip("\n")],
                        delimiter="\t",
                        quotechar='"',
                    )
                )
                metadata["accessions"] = parts[1:]

            elif line.startswith("!Sample_characteristics_ch1"):
                parts = next(
                    csv.reader(
                        [line.rstrip("\n")],
                        delimiter="\t",
                        quotechar='"',
                    )
                )
                metadata.setdefault(
                    "characteristics",
                    []
                ).append(parts[1:])

    return metadata


def load_expression(path):
    lines = []

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        inside = False

        for line in f:

            if line.startswith("!series_matrix_table_begin"):
                inside = True
                continue

            if line.startswith("!series_matrix_table_end"):
                break

            if inside and line.strip():
                lines.append(line)

    return pd.read_csv(
        StringIO("".join(lines)),
        sep="\t",
        quotechar='"',
    )


print("=" * 80)
print("CLEAN FROZEN-SIGNATURE EXTERNAL VALIDATION")
print("=" * 80)

# ================================================================
# GSE8841 DISCOVERY
# ================================================================

disc = pd.read_csv(DISCOVERY, low_memory=False)

disc = disc[
    disc["histotype"]
    .astype(str)
    .str.lower()
    .isin(["serous", "mucinous"])
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
    c for c in disc.columns
    if c not in metadata_cols
]

probe_to_gene_8841 = load_gpl5689(GPL5689)

disc_gene = {}

for gene in GENES:

    probes = [
        p for p in disc_probes
        if probe_to_gene_8841.get(str(p), "") == gene
    ]

    if not probes:
        raise RuntimeError(
            f"GSE8841 has no probe mapping for {gene}"
        )

    vals = (
        disc[probes]
        .apply(pd.to_numeric, errors="coerce")
        .median(axis=1)
    )

    disc_gene[gene] = vals

disc_X = pd.DataFrame(disc_gene)

if disc_X.isna().any().any():
    raise RuntimeError(
        "Missing values in GSE8841 signature"
    )

y_disc = (
    disc["label"] == "stage_I_cancer"
).astype(int)

print()
print("DISCOVERY:")
print("  patients:", len(disc_X))
print("  stage-I cancer:", int(y_disc.sum()))
print("  borderline:", int((y_disc == 0).sum()))

# ================================================================
# GSE12172 METADATA
# ================================================================

meta = load_geo_metadata(TEST)

accessions = meta["accessions"]
char_rows = meta["characteristics"]

print()
print("GSE12172 characteristic rows:", len(char_rows))

# IMPORTANT:
# Your audit established that ROW 1 contains exactly:
#   Malignant (60)
#   LMP (30)

diagnosis_row = char_rows[0]

print(
    "Diagnosis row unique values:",
    sorted(set(str(x).strip() for x in diagnosis_row))
)

labels = []

for value in diagnosis_row:

    value = str(value).strip().lower()

    if value == "lmp":
        labels.append("borderline")

    elif value == "malignant":
        labels.append("malignant")

    else:
        labels.append("other")

labels = np.array(labels)

print()
print("GSE12172 LABEL COUNTS:")

for label, count in zip(
    *np.unique(labels, return_counts=True)
):
    print(f"  {label}: {count}")

if int((labels == "borderline").sum()) != 30:
    raise RuntimeError(
        "Expected exactly 30 LMP samples."
    )

if int((labels == "malignant").sum()) != 60:
    raise RuntimeError(
        "Expected exactly 60 malignant samples."
    )

# ================================================================
# GSE12172 EXPRESSION
# ================================================================

expr = load_expression(TEST)

probe_col = expr.columns[0]
expr[probe_col] = (
    expr[probe_col]
    .astype(str)
    .str.strip()
)

# Only the 90 desired samples.
test_accessions = [
    a for a, label in zip(accessions, labels)
    if label in {"borderline", "malignant"}
]

test_labels = np.array([
    label for label in labels
    if label in {"borderline", "malignant"}
])

probe_to_gene_12172 = load_gpl570(GPL570)

test_gene = {}

print()
print("=" * 80)
print("GSE12172 GENE MAPPING")
print("=" * 80)

for gene in GENES:

    probes = [
        p for p in expr[probe_col]
        if probe_to_gene_12172.get(p, "") == gene
    ]

    print(
        f"{gene:15} probes={len(probes)}"
    )

    if not probes:
        raise RuntimeError(
            f"GSE12172 has no probe for {gene}"
        )

    sub = expr[
        expr[probe_col].isin(probes)
    ]

    vals = (
        sub[test_accessions]
        .apply(pd.to_numeric, errors="coerce")
        .median(axis=0)
    )

    test_gene[gene] = vals.to_numpy(dtype=float)

test_X = pd.DataFrame(
    test_gene,
    index=test_accessions,
)

if test_X.isna().any().any():
    raise RuntimeError(
        "Missing values in GSE12172 signature"
    )

# ================================================================
# SCALE USING DISCOVERY ONLY
# ================================================================

scaler = StandardScaler()

X_disc = scaler.fit_transform(
    disc_X[GENES]
)

X_test = scaler.transform(
    test_X[GENES]
)

# ================================================================
# FINAL TRAINING
# ================================================================

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

# 1 = malignant
prediction = (
    prob >= 0.5
).astype(int)

y_test = (
    test_labels == "malignant"
).astype(int)

# ================================================================
# EXTERNAL METRICS
# ================================================================

auc = roc_auc_score(
    y_test,
    prob,
)

accuracy = accuracy_score(
    y_test,
    prediction,
)

balanced_accuracy = balanced_accuracy_score(
    y_test,
    prediction,
)

tn, fp, fn, tp = confusion_matrix(
    y_test,
    prediction,
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
print("GSE12172 EXTERNAL VALIDATION")
print("=" * 80)

print()
print("Patients:", len(y_test))
print("Borderline/LMP:", int((y_test == 0).sum()))
print("Malignant:", int((y_test == 1).sum()))

print()
print(f"ROC AUC:            {auc:.6f}")
print(f"Accuracy:            {accuracy:.6f}")
print(f"Balanced accuracy:   {balanced_accuracy:.6f}")
print(f"Sensitivity:         {sensitivity:.6f}")
print(f"Specificity:         {specificity:.6f}")

print()
print("Confusion matrix:")
print(f"  TN={tn} FP={fp}")
print(f"  FN={fn} TP={tp}")

print()
print("=" * 80)
print("IMPORTANT")
print("=" * 80)

print(
    "13 genes were frozen from GSE8841 BEFORE testing GSE12172."
)

print(
    "No gene selection was performed on GSE12172."
)

print(
    "This is tissue-expression external validation,"
)

print(
    "not validation of the serum-miRNA diagnostic model."
)

print()
print("=" * 80)
print("DONE")
print("=" * 80)
