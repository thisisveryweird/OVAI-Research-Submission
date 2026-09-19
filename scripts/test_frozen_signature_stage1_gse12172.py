import gzip
import csv
from io import StringIO
import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    roc_auc_score,
    balanced_accuracy_score,
    confusion_matrix,
)
from sklearn.preprocessing import StandardScaler


ROOT = "data"

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

DISCOVERY = f"{ROOT}/processed/gse8841_borderline_expression.csv.gz"
TEST = f"{ROOT}/external/gse12172/GSE12172_series_matrix.txt.gz"
GPL5689 = f"{ROOT}/external/gse_platforms/GPL5689.annot.gz"
GPL570 = f"{ROOT}/external/gse_platforms/GPL570.txt"


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
    tab = pd.DataFrame(
        [r.split("\t") for r in rows[1:]],
        columns=header,
    )

    return dict(
        zip(
            tab["ID"].astype(str).str.strip(),
            tab["Gene symbol"].fillna("").astype(str)
            .str.strip().str.upper(),
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
    tab = pd.DataFrame(
        [r.split("\t") for r in rows[1:]],
        columns=header,
    )

    return dict(
        zip(
            tab["ID"].astype(str).str.strip(),
            tab["Gene Symbol"].fillna("").astype(str)
            .str.strip().str.upper(),
        )
    )


def load_meta(path):
    meta = {}

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
                meta["accessions"] = parts[1:]

            elif line.startswith("!Sample_characteristics_ch1"):
                parts = next(
                    csv.reader(
                        [line.rstrip("\n")],
                        delimiter="\t",
                        quotechar='"',
                    )
                )
                meta.setdefault("chars", []).append(parts[1:])

    return meta


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
print("FROZEN SIGNATURE: BORDERLINE vs STAGE-I CANCER")
print("=" * 80)

# ================================================================
# DISCOVERY
# ================================================================

disc = pd.read_csv(
    DISCOVERY,
    low_memory=False,
)

disc = disc[
    disc["histotype"].astype(str).str.lower().isin(
        ["serous", "mucinous"]
    )
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

probe_cols = [
    c for c in disc.columns
    if c not in metadata_cols
]

ann8841 = load_gpl5689(GPL5689)

disc_gene = {}

for gene in GENES:

    probes = [
        p for p in probe_cols
        if ann8841.get(str(p), "") == gene
    ]

    if not probes:
        raise RuntimeError(
            f"Missing discovery probe for {gene}"
        )

    disc_gene[gene] = (
        disc[probes]
        .apply(pd.to_numeric, errors="coerce")
        .median(axis=1)
    )

disc_X = pd.DataFrame(disc_gene)

if disc_X.isna().any().any():
    raise RuntimeError(
        "Missing discovery values."
    )

y_disc = (
    disc["label"] == "stage_I_cancer"
).astype(int).to_numpy()

print()
print("Discovery patients:", len(disc))
print("Stage-I cancer:", int(y_disc.sum()))
print("Borderline:", int((y_disc == 0).sum()))

# ================================================================
# GSE12172 METADATA
# ================================================================

meta = load_meta(TEST)

accessions = meta["accessions"]
chars = meta["chars"]

# Diagnosis row
diagnosis_row = None

for row in chars:
    vals = {
        str(x).strip().lower()
        for x in row
    }

    if "lmp" in vals and "malignant" in vals:
        diagnosis_row = row
        break

if diagnosis_row is None:
    raise RuntimeError(
        "Could not find diagnosis row."
    )

# FIGO stage row
#
# GSE12172 audit established that the 4th characteristic row
# contains the FIGO stage values:
#   Stage IA, Stage IB, Stage IC, Stage IIA, ...
#
# Use that known row directly rather than guessing from text.

stage_row = chars[3]

def normalize_stage(value):

    value = str(value).strip().upper()

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

# ---------------------------------------------------------------
# Keep:
#   ALL 30 LMP/borderline
#   ONLY Stage IA / IB / IC malignant
# ---------------------------------------------------------------

keep = (
    labels == "borderline"
) | (
    (labels == "malignant")
    & np.isin(
        stages,
        ["I"],
    )
)

test_ids = np.array(accessions)[keep]
test_labels = labels[keep]
test_stages = stages[keep]

print()
print("=" * 80)
print("GSE12172 STAGE-I TEST COHORT")
print("=" * 80)

print(
    "Borderline/LMP:",
    int((test_labels == "borderline").sum())
)

print(
    "Stage-I malignant:",
    int(
        (
            test_labels == "malignant"
        ).sum()
    )
)

print()
print("Stage distribution among malignant:")
for s, n in zip(
    *np.unique(
        test_stages[
            test_labels == "malignant"
        ],
        return_counts=True,
    )
):
    print(
        f"  {s}: {n}"
    )

expected_borderline = 30
expected_stage1 = 16

if int(
    (test_labels == "borderline").sum()
) != expected_borderline:
    raise RuntimeError(
        "Expected 30 borderline/LMP samples."
    )

if int(
    (test_labels == "malignant").sum()
) != expected_stage1:
    raise RuntimeError(
        "Expected 16 Stage-I malignant samples."
    )

# ================================================================
# EXPRESSION
# ================================================================

expr = load_expression(TEST)

probe_col = expr.columns[0]

expr[probe_col] = (
    expr[probe_col]
    .astype(str)
    .str.strip()
)

ann570 = load_gpl570(GPL570)

test_gene = {}

for gene in GENES:

    probes = [
        p
        for p in expr[probe_col]
        if ann570.get(p, "") == gene
    ]

    if not probes:
        continue

    sub = expr[
        expr[probe_col].isin(probes)
    ]

    test_gene[gene] = (
        sub[test_ids]
        .apply(pd.to_numeric, errors="coerce")
        .median(axis=0)
        .to_numpy(dtype=float)
    )

available = [
    g
    for g in GENES
    if g in test_gene
]

print()
print(
    "Transferable genes:",
    len(available)
)

print(
    "Excluded:",
    [
        g for g in GENES
        if g not in available
    ]
)

if len(available) < 10:
    raise RuntimeError(
        "Too few transferable genes."
    )

test_X = pd.DataFrame(
    {
        g: test_gene[g]
        for g in available
    },
    index=test_ids,
)

disc_use = disc_X[available]

if test_X.isna().any().any():
    raise RuntimeError(
        "Missing external expression values."
    )

# ================================================================
# TRAIN ONLY ON GSE8841
# ================================================================

scaler = StandardScaler()

X_disc = scaler.fit_transform(
    disc_use
)

X_test = scaler.transform(
    test_X[available]
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

prob = clf.predict_proba(
    X_test
)[:, 1]

# 1 = malignant
y_test = (
    test_labels == "malignant"
).astype(int)

# ================================================================
# AUC ONLY — NO THRESHOLD TUNING
# ================================================================

auc = roc_auc_score(
    y_test,
    prob,
)

print()
print("=" * 80)
print("FINAL STAGE-I EXTERNAL RESULT")
print("=" * 80)

print()
print(
    "Patients:", len(y_test)
)

print(
    "Borderline:", int((y_test == 0).sum())
)

print(
    "Stage-I malignant:", int((y_test == 1).sum())
)

print()
print(
    f"ROC AUC: {auc:.6f}"
)

print()
print(
    "This is a frozen external discrimination test."
)

print(
    "No genes were selected using GSE12172."
)

print(
    "No threshold was selected using GSE12172."
)

print("=" * 80)
print("DONE")
print("=" * 80)
