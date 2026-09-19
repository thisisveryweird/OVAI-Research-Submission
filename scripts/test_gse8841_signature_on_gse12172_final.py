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
SIGNATURE = ROOT / "data/processed/gse8841_frozen_borderline_genes.txt"


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
    table = pd.DataFrame(
        [r.split("\t") for r in rows[1:]],
        columns=header,
    )

    return dict(
        zip(
            table["ID"].astype(str).str.strip(),
            table["Gene symbol"].fillna("")
            .astype(str)
            .str.strip()
            .str.upper(),
        )
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
    table = pd.DataFrame(
        [r.split("\t") for r in rows[1:]],
        columns=header,
    )

    return dict(
        zip(
            table["ID"].astype(str).str.strip(),
            table["Gene Symbol"].fillna("")
            .astype(str)
            .str.strip()
            .str.upper(),
        )
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
                    "characteristics", []
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
        low_memory=False,
    )


print("=" * 80)
print("FINAL FROZEN-SIGNATURE EXTERNAL VALIDATION")
print("=" * 80)

# ================================================================
# Frozen genes
# ================================================================

genes = [
    x.strip().upper()
    for x in SIGNATURE.read_text().splitlines()
    if x.strip()
]

print()
print("Frozen discovery genes:", len(genes))

for gene in genes:
    print(" ", gene)

# ================================================================
# GSE8841 discovery data
# ================================================================

disc = pd.read_csv(
    DISCOVERY,
    low_memory=False,
)

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

ann8841 = load_gpl5689(GPL5689)

disc_gene_values = {}

for gene in genes:

    probes = [
        p
        for p in disc_probes
        if ann8841.get(str(p), "") == gene
    ]

    if not probes:
        raise RuntimeError(
            f"Discovery gene {gene} has no GSE8841 probe."
        )

    disc_gene_values[gene] = (
        disc[probes]
        .apply(pd.to_numeric, errors="coerce")
        .median(axis=1)
    )

disc_X = pd.DataFrame(
    disc_gene_values
)

if disc_X.isna().any().any():
    raise RuntimeError(
        "Missing values in GSE8841 signature."
    )

y_disc = (
    disc["label"] == "stage_I_cancer"
).astype(int)

print()
print("Discovery patients:", len(disc_X))
print("Stage-I cancer:", int(y_disc.sum()))
print("Borderline:", int((y_disc == 0).sum()))

# ================================================================
# GSE12172
# ================================================================

meta = load_geo_metadata(TEST)
expr = load_expression(TEST)

accessions = meta["accessions"]
char_rows = meta["characteristics"]

print()
print("GSE12172 samples:", len(accessions))
print("Characteristic rows:", len(char_rows))

# Diagnosis is the row containing exact LMP/Malignant values.
diagnosis_row = None

for row in char_rows:
    normalized = {
        str(x).strip().lower()
        for x in row
    }

    if "lmp" in normalized and "malignant" in normalized:
        diagnosis_row = row
        break

if diagnosis_row is None:
    raise RuntimeError(
        "Could not locate GSE12172 diagnosis row."
    )

labels = np.array([
    "borderline" if str(x).strip().lower() == "lmp"
    else "malignant" if str(x).strip().lower() == "malignant"
    else "other"
    for x in diagnosis_row
])

print()
print("GSE12172 labels:")
for value, count in zip(
    *np.unique(labels, return_counts=True)
):
    print(f"  {value}: {count}")

keep = np.isin(
    labels,
    ["borderline", "malignant"],
)

test_accessions = np.array(accessions)[keep]
test_labels = labels[keep]

if len(test_accessions) != 90:
    raise RuntimeError(
        f"Expected 90 usable test samples, got {len(test_accessions)}."
    )

# ================================================================
# GPL570 mapping
# ================================================================

ann570 = load_gpl570(GPL570)

probe_col = expr.columns[0]

expr[probe_col] = (
    expr[probe_col]
    .astype(str)
    .str.strip()
)

test_gene_values = {}
available_genes = []
missing_genes = []

for gene in genes:

    probes = [
        p
        for p in expr[probe_col]
        if ann570.get(p, "") == gene
    ]

    if not probes:
        missing_genes.append(gene)
        print(
            f"{gene:15} probes=0  -> EXCLUDED"
        )
        continue

    available_genes.append(gene)

    print(
        f"{gene:15} probes={len(probes)}"
    )

    sub = expr[
        expr[probe_col].isin(probes)
    ]

    test_gene_values[gene] = (
        sub[test_accessions]
        .apply(pd.to_numeric, errors="coerce")
        .median(axis=0)
        .to_numpy(dtype=float)
    )

print()
print("=" * 80)
print("TRANSFERABLE SIGNATURE")
print("=" * 80)

print("Available genes:", len(available_genes))
print("Excluded genes:", missing_genes)

if len(available_genes) < 10:
    raise RuntimeError(
        "Fewer than 10 frozen genes are available on GPL570."
    )

# ================================================================
# Build matrices
# ================================================================

test_X = pd.DataFrame(
    {
        gene: test_gene_values[gene]
        for gene in available_genes
    },
    index=test_accessions,
)

disc_X_use = disc_X[
    available_genes
]

if disc_X_use.isna().any().any():
    raise RuntimeError(
        "Missing discovery values."
    )

if test_X.isna().any().any():
    raise RuntimeError(
        "Missing external values."
    )

# ================================================================
# Standardization
#
# Discovery cohort defines scale.
# ================================================================

scaler = StandardScaler()

X_disc = scaler.fit_transform(
    disc_X_use
)

X_test = scaler.transform(
    test_X
)

# ================================================================
# Frozen classifier
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

# Probability of malignant/stage-I cancer
prob = clf.predict_proba(
    X_test
)[:, 1]

pred = (
    prob >= 0.5
).astype(int)

y_test = (
    test_labels == "malignant"
).astype(int)

# ================================================================
# External metrics
# ================================================================

auc = roc_auc_score(
    y_test,
    prob,
)

accuracy = accuracy_score(
    y_test,
    pred,
)

balanced_accuracy = balanced_accuracy_score(
    y_test,
    pred,
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
print("Patients:", len(y_test))
print("Borderline/LMP:", int((y_test == 0).sum()))
print("Malignant:", int((y_test == 1).sum()))

print()
print(f"Transferable genes:   {len(available_genes)}")
print(f"Excluded genes:       {missing_genes}")

print()
print(f"ROC AUC:              {auc:.6f}")
print(f"Accuracy:             {accuracy:.6f}")
print(f"Balanced accuracy:    {balanced_accuracy:.6f}")
print(f"Sensitivity:          {sensitivity:.6f}")
print(f"Specificity:          {specificity:.6f}")

print()
print("Confusion matrix:")
print(f"  TN={tn} FP={fp}")
print(f"  FN={fn} TP={tp}")

print()
print("=" * 80)
print("DONE")
print("=" * 80)

print()
print("IMPORTANT:")
print(
    "Genes were frozen using GSE8841 only."
)
print(
    "GSE12172 was not used for feature selection."
)
print(
    "IFI30 was excluded because GPL570 did not map a probe to it."
)
print(
    "This is tissue-expression external validation."
)
