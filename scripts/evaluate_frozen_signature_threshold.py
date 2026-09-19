import gzip
import csv
from io import StringIO

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    balanced_accuracy_score,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
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
            tab["Gene symbol"]
            .fillna("")
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


def geo_meta(path):
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


def geo_expression(path):
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
print("FROZEN SIGNATURE THRESHOLD ANALYSIS")
print("=" * 80)

# ================================================================
# GSE8841 DISCOVERY
# ================================================================

disc = pd.read_csv(
    f"{ROOT}/processed/gse8841_borderline_expression.csv.gz",
    low_memory=False,
)

disc = disc[
    disc["histotype"]
    .astype(str)
    .str.lower()
    .isin(["serous", "mucinous"])
].copy()

meta_cols = [
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

probes = [
    c for c in disc.columns
    if c not in meta_cols
]

probe_map = load_gpl5689(
    f"{ROOT}/external/gse_platforms/GPL5689.annot.gz"
)

disc_genes = {}

for gene in GENES:

    gene_probes = [
        p for p in probes
        if probe_map.get(str(p), "") == gene
    ]

    if not gene_probes:
        print(
            f"WARNING: no discovery probe for {gene}"
        )
        continue

    disc_genes[gene] = (
        disc[gene_probes]
        .apply(pd.to_numeric, errors="coerce")
        .median(axis=1)
    )

disc_X = pd.DataFrame(disc_genes)

available_disc = [
    g for g in GENES
    if g in disc_X.columns
]

if len(available_disc) < 10:
    raise RuntimeError(
        "Too few discovery genes available."
    )

if disc_X[available_disc].isna().any().any():
    raise RuntimeError(
        "Missing discovery values."
    )

y_disc = (
    disc["label"] == "stage_I_cancer"
).astype(int).to_numpy()

print()
print("DISCOVERY")
print("Patients:", len(disc_X))
print("Stage-I cancer:", int(y_disc.sum()))
print("Borderline:", int((y_disc == 0).sum()))

# ================================================================
# DISCOVERY OOF PROBABILITIES
# ================================================================

cv = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=12345,
)

oof_prob = np.zeros(
    len(disc_X),
    dtype=float,
)

for train_idx, val_idx in cv.split(
    disc_X,
    y_disc,
):

    scaler = StandardScaler()

    X_train = scaler.fit_transform(
        disc_X.iloc[train_idx][available_disc]
    )

    X_val = scaler.transform(
        disc_X.iloc[val_idx][available_disc]
    )

    clf = LogisticRegression(
        C=0.1,
        class_weight="balanced",
        max_iter=5000,
        solver="liblinear",
    )

    clf.fit(
        X_train,
        y_disc[train_idx],
    )

    oof_prob[val_idx] = clf.predict_proba(
        X_val
    )[:, 1]

discovery_auc = roc_auc_score(
    y_disc,
    oof_prob,
)

# ================================================================
# CHOOSE THRESHOLD USING DISCOVERY ONLY
# ================================================================

thresholds = np.linspace(
    0.01,
    0.99,
    197,
)

best_threshold = None
best_balanced = -np.inf

for threshold in thresholds:

    pred = (
        oof_prob >= threshold
    ).astype(int)

    score = balanced_accuracy_score(
        y_disc,
        pred,
    )

    if score > best_balanced:

        best_balanced = score
        best_threshold = threshold

print()
print("=" * 80)
print("DISCOVERY THRESHOLD")
print("=" * 80)

print(
    f"Discovery OOF AUC:       {discovery_auc:.6f}"
)

print(
    f"Frozen threshold:        {best_threshold:.3f}"
)

print(
    f"Discovery balanced acc:  {best_balanced:.6f}"
)

# ================================================================
# FINAL DISCOVERY MODEL
# ================================================================

scaler = StandardScaler()

X_disc = scaler.fit_transform(
    disc_X[available_disc]
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
# GSE12172
# ================================================================

TEST = (
    f"{ROOT}/external/gse12172/"
    "GSE12172_series_matrix.txt.gz"
)

meta = geo_meta(TEST)

accessions = meta["accessions"]
chars = meta["chars"]

# Find exact diagnosis row.
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
        "Could not find GSE12172 diagnosis row."
    )

labels = np.array([
    "borderline"
    if str(x).strip().lower() == "lmp"
    else
    "malignant"
    if str(x).strip().lower() == "malignant"
    else
    "other"
    for x in diagnosis_row
])

keep = np.isin(
    labels,
    ["borderline", "malignant"],
)

test_ids = np.array(accessions)[keep]
test_labels = labels[keep]

print()
print("=" * 80)
print("GSE12172 EXTERNAL COHORT")
print("=" * 80)

print(
    "Borderline:",
    int((test_labels == "borderline").sum())
)

print(
    "Malignant:",
    int((test_labels == "malignant").sum())
)

if len(test_ids) != 90:
    raise RuntimeError(
        "Expected 90 GSE12172 test samples."
    )

# ================================================================
# GSE12172 EXPRESSION + GPL570
# ================================================================

expr = geo_expression(TEST)

probe_map_test = load_gpl570(
    f"{ROOT}/external/gse_platforms/GPL570.txt"
)

probe_col = expr.columns[0]

expr[probe_col] = (
    expr[probe_col]
    .astype(str)
    .str.strip()
)

test_genes = {}

for gene in available_disc:

    gene_probes = [
        p
        for p in expr[probe_col]
        if probe_map_test.get(p, "") == gene
    ]

    if not gene_probes:
        continue

    sub = expr[
        expr[probe_col].isin(gene_probes)
    ]

    test_genes[gene] = (
        sub[test_ids]
        .apply(pd.to_numeric, errors="coerce")
        .median(axis=0)
        .to_numpy(dtype=float)
    )

available_test = [
    g for g in available_disc
    if g in test_genes
]

excluded = [
    g for g in available_disc
    if g not in available_test
]

print()
print(
    "Transferable genes:",
    len(available_test)
)

print(
    "Excluded:",
    excluded
)

if len(available_test) < 10:
    raise RuntimeError(
        "Too few transferable genes."
    )

test_X = pd.DataFrame(
    {
        g: test_genes[g]
        for g in available_test
    },
    index=test_ids,
)

# Refit using the COMMON genes only.
common_scaler = StandardScaler()

X_disc_common = common_scaler.fit_transform(
    disc_X[available_test]
)

X_test = common_scaler.transform(
    test_X[available_test]
)

common_clf = LogisticRegression(
    C=0.1,
    class_weight="balanced",
    max_iter=5000,
    solver="liblinear",
)

common_clf.fit(
    X_disc_common,
    y_disc,
)

prob = common_clf.predict_proba(
    X_test
)[:, 1]

y_test = (
    test_labels == "malignant"
).astype(int)

# ================================================================
# EXTERNAL VALIDATION
# ================================================================

auc = roc_auc_score(
    y_test,
    prob,
)

pred = (
    prob >= best_threshold
).astype(int)

balanced = balanced_accuracy_score(
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
print("GSE12172 EXTERNAL VALIDATION")
print("=" * 80)

print()
print(
    f"ROC AUC:              {auc:.6f}"
)

print(
    f"Frozen threshold:     {best_threshold:.3f}"
)

print(
    f"Balanced accuracy:    {balanced:.6f}"
)

print(
    f"Sensitivity:          {sensitivity:.6f}"
)

print(
    f"Specificity:          {specificity:.6f}"
)

print()
print("Confusion matrix:")
print(f"  TN={tn} FP={fp}")
print(f"  FN={fn} TP={tp}")

print()
print("=" * 80)
print("DONE")
print("=" * 80)

print()
print("The threshold was chosen from GSE8841 OOF predictions only.")
print("GSE12172 was not used to select the threshold or genes.")
