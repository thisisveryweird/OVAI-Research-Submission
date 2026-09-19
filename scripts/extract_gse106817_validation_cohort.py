from pathlib import Path
import csv
import gzip
import re

import numpy as np
import pandas as pd


ROOT = Path.home() / "ovarian-ai"

MATRIX = (
    ROOT
    / "data"
    / "external"
    / "gse106817"
    / "GSE106817_series_matrix.txt.gz"
)

MAPPING = (
    ROOT
    / "reports"
    / "tables"
    / "historical_mirna_to_mimat.csv"
)

OUT = (
    ROOT
    / "data"
    / "processed"
    / "gse106817_external_validation.csv.gz"
)


def parse_quoted_fields(line):
    parts = next(
        csv.reader(
            [line.rstrip("\n")],
            delimiter="\t",
            quotechar='"',
        )
    )
    return parts[1:]


print("=" * 80)
print("BUILDING GSE106817 EXTERNAL VALIDATION COHORT")
print("=" * 80)

# ------------------------------------------------------------------
# Read GEO metadata
# ------------------------------------------------------------------

metadata = {}

with gzip.open(
    MATRIX,
    "rt",
    encoding="utf-8",
    errors="replace",
) as f:

    for line in f:

        if line.startswith("!series_matrix_table_begin"):
            break

        if line.startswith("!Sample_"):

            key = line.split("\t", 1)[0]

            metadata[key] = parse_quoted_fields(line)


titles = metadata["!Sample_title"]
accessions = metadata["!Sample_geo_accession"]
characteristics = metadata["!Sample_characteristics_ch1"]

if not (
    len(titles)
    == len(accessions)
    == len(characteristics)
):
    raise RuntimeError(
        "Sample metadata lengths do not agree."
    )

n = len(titles)

print()
print("Total GEO samples:", n)

# ------------------------------------------------------------------
# Classify samples from title
# ------------------------------------------------------------------

def classify(title):

    t = title.strip().lower()

    if t.startswith("ovarian cancer"):
        return "cancer"

    if t.startswith("borderline ovarian tumor"):
        return "borderline"

    if t.startswith("benign ovarian disease"):
        return "benign"

    if t.startswith("non-cancer"):
        return "control"

    return None


labels = [
    classify(x)
    for x in titles
]

groups = pd.Series(labels).value_counts()

print()
print("=" * 80)
print("VALIDATION GROUPS")
print("=" * 80)

for group in [
    "cancer",
    "borderline",
    "benign",
    "control",
]:

    print(
        f"{group:12} "
        f"{int(groups.get(group, 0))}"
    )

# ------------------------------------------------------------------
# Build age parser
# ------------------------------------------------------------------

def parse_age(text):

    t = str(text)

    patterns = [
        r"age[^:]*:\s*([0-9]+(?:\.[0-9]+)?)",
        r"age[^:]*=\s*([0-9]+(?:\.[0-9]+)?)",
        r"years?[^:]*:\s*([0-9]+(?:\.[0-9]+)?)",
    ]

    for pattern in patterns:

        m = re.search(
            pattern,
            t,
            flags=re.I,
        )

        if m:
            return float(m.group(1))

    return np.nan


ages = [
    parse_age(x)
    for x in characteristics
]

print()
print("=" * 80)
print("AGE")
print("=" * 80)

print(
    "Observed:",
    sum(pd.notna(x) for x in ages),
)

print(
    "Missing:",
    sum(pd.isna(x) for x in ages),
)

# ------------------------------------------------------------------
# Select relevant samples
# ------------------------------------------------------------------

keep = [
    label in {
        "cancer",
        "borderline",
        "benign",
        "control",
    }
    for label in labels
]

selected_accessions = [
    accessions[i]
    for i in range(n)
    if keep[i]
]

selected_titles = [
    titles[i]
    for i in range(n)
    if keep[i]
]

selected_labels = [
    labels[i]
    for i in range(n)
    if keep[i]
]

selected_ages = [
    ages[i]
    for i in range(n)
    if keep[i]
]

print()
print(
    "Relevant validation samples:",
    len(selected_accessions),
)

# ------------------------------------------------------------------
# Load mapped miRNAs
# ------------------------------------------------------------------

mapping = pd.read_csv(MAPPING)

mapping = mapping[
    mapping["matched"] == True
].copy()

mimat_by_name = dict(
    zip(
        mapping["historical_name"],
        mapping["mimat"].astype(str),
    )
)

historical_names = list(
    mapping["historical_name"]
)

# ------------------------------------------------------------------
# Find data table
# ------------------------------------------------------------------

with gzip.open(
    MATRIX,
    "rt",
    encoding="utf-8",
    errors="replace",
) as f:

    for line in f:

        if line.startswith(
            "!series_matrix_table_begin"
        ):
            break

    rows = []

    for line in f:

        if line.startswith(
            "!series_matrix_table_end"
        ):
            break

        rows.append(line)

from io import StringIO

matrix = pd.read_csv(
    StringIO("".join(rows)),
    sep="\t",
    index_col=0,
)

matrix.index = (
    matrix.index
    .astype(str)
    .str.strip()
    .str.replace('"', "", regex=False)
)

print()
print(
    "Expression feature rows:",
    len(matrix),
)

# ------------------------------------------------------------------
# Match the 172 MIMAT IDs
# ------------------------------------------------------------------

available = set(
    matrix.index
)

matched_names = []
matched_mimat = []

for name in historical_names:

    mimat = mimat_by_name[name]

    if mimat in available:

        matched_names.append(name)
        matched_mimat.append(mimat)

print()
print("=" * 80)
print("FEATURE MATCH")
print("=" * 80)

print(
    "Mapped historical miRNAs:",
    len(historical_names),
)

print(
    "Available in GSE106817:",
    len(matched_names),
)

# ------------------------------------------------------------------
# Extract selected samples
# ------------------------------------------------------------------

sample_to_col = {
    str(col).strip().strip('"'): col
    for col in matrix.columns
}

missing_samples = [
    x
    for x in selected_accessions
    if x not in sample_to_col
]

if missing_samples:

    raise RuntimeError(
        "Some GEO accessions were not found "
        "in the expression matrix: "
        + str(missing_samples[:10])
    )

cols = [
    sample_to_col[x]
    for x in selected_accessions
]

X = matrix.loc[
    matched_mimat,
    cols,
].apply(
    pd.to_numeric,
    errors="coerce",
).T

X.columns = matched_names

X.insert(
    0,
    "sample_id",
    selected_accessions,
)

X.insert(
    1,
    "title",
    selected_titles,
)

X.insert(
    2,
    "group",
    selected_labels,
)

X.insert(
    3,
    "age",
    selected_ages,
)

# ------------------------------------------------------------------
# Binary labels
#
# Primary diagnostic endpoint:
# cancer = 1
# benign/borderline/control = 0
# ------------------------------------------------------------------

X["binary_label"] = (
    X["group"] == "cancer"
).astype(int)

# ------------------------------------------------------------------
# Missingness
# ------------------------------------------------------------------

feature_cols = matched_names

bad = X[feature_cols].isna().sum()

print()
print("=" * 80)
print("FINAL DATA VALIDATION")
print("=" * 80)

print(
    "Patients:",
    len(X),
)

print(
    "Features:",
    len(feature_cols),
)

print(
    "Cancer:",
    int(X["binary_label"].sum()),
)

print(
    "Non-cancer:",
    int((X["binary_label"] == 0).sum()),
)

print(
    "Age observed:",
    int(X["age"].notna().sum()),
)

print(
    "Age missing:",
    int(X["age"].isna().sum()),
)

print(
    "Samples with any missing miRNA:",
    int(
        X[feature_cols]
        .isna()
        .any(axis=1)
        .sum()
    ),
)

print()
print("Group counts:")

print(
    X["group"]
    .value_counts()
    .to_string()
)

# ------------------------------------------------------------------
# Save
# ------------------------------------------------------------------

OUT.parent.mkdir(
    parents=True,
    exist_ok=True,
)

X.to_csv(
    OUT,
    index=False,
    compression="gzip",
)

print()
print("=" * 80)
print("SAVED")
print("=" * 80)

print(OUT)
