#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


ROOT = Path.home() / "ovarian-ai"

INPUT = (
    ROOT
    / "data"
    / "processed"
    / "tcga_ov_304_paired_molecular.csv.gz"
)

OUTPUT = (
    ROOT
    / "data"
    / "processed"
    / "tcga_ov_304_clean_multimodal.csv.gz"
)

MISSINGNESS_THRESHOLD = 0.20


print("=" * 80)
print("BUILDING CLEAN TCGA-OV MULTIMODAL DATASET")
print("=" * 80)

df = pd.read_csv(INPUT)

mirna_cols = [
    c for c in df.columns
    if c.startswith("hsa-")
]

rna_cols = [
    "PAX8",
    "EPCAM",
    "KRT7",
    "KRT8",
    "KRT18",
    "MUC16",
    "WFDC2",
    "FOLR1",
    "WT1",
    "MKI67",
]

print()
print("Input patients:", len(df))
print("Input miRNAs:", len(mirna_cols))
print("RNA features:", len(rna_cols))

# ------------------------------------------------------------------
# Remove extremely sparse miRNAs.
# ------------------------------------------------------------------

missing_fraction = df[mirna_cols].isna().mean()

keep_mirnas = [
    c
    for c in mirna_cols
    if missing_fraction[c] <= MISSINGNESS_THRESHOLD
]

removed_mirnas = [
    c
    for c in mirna_cols
    if missing_fraction[c] > MISSINGNESS_THRESHOLD
]

print()
print("=" * 80)
print("FEATURE FILTER")
print("=" * 80)

print(
    "Missingness threshold:",
    f"{MISSINGNESS_THRESHOLD:.0%}",
)

print(
    "miRNAs kept:",
    len(keep_mirnas),
)

print(
    "miRNAs removed:",
    len(removed_mirnas),
)

# ------------------------------------------------------------------
# Keep only required columns.
# ------------------------------------------------------------------

required = (
    ["sample_id", "age"]
    + keep_mirnas
    + rna_cols
)

missing_required = [
    c
    for c in required
    if c not in df.columns
]

if missing_required:
    raise RuntimeError(
        "Missing required columns: "
        + str(missing_required)
    )

clean = df[required].copy()

# ------------------------------------------------------------------
# Report missingness BEFORE imputation.
# ------------------------------------------------------------------

before_missing = int(
    clean[keep_mirnas]
    .isna()
    .sum()
    .sum()
)

print()
print(
    "Missing retained-miRNA values before imputation:",
    before_missing,
)

# ------------------------------------------------------------------
# Median imputation for sparse remaining missing values.
#
# This is ONLY dataset preparation. When we eventually train a
# predictive model, imputation must be performed within training
# folds to prevent leakage.
# ------------------------------------------------------------------

for col in keep_mirnas:

    if clean[col].isna().any():

        median = clean[col].median()

        if pd.isna(median):
            raise RuntimeError(
                f"Cannot impute {col}: median is NaN."
            )

        clean[col] = clean[col].fillna(
            median
        )

after_missing = int(
    clean[keep_mirnas]
    .isna()
    .sum()
    .sum()
)

# ------------------------------------------------------------------
# Validate RNA and age.
# ------------------------------------------------------------------

rna_missing = int(
    clean[rna_cols]
    .isna()
    .sum()
    .sum()
)

age_missing = int(
    clean["age"].isna().sum()
)

print()
print("=" * 80)
print("FINAL VALIDATION")
print("=" * 80)

print("Patients:", len(clean))
print("miRNAs:", len(keep_mirnas))
print("RNA features:", len(rna_cols))
print("Age missing:", age_missing)
print(
    "Retained-miRNA missing before:",
    before_missing,
)
print(
    "Retained-miRNA missing after:",
    after_missing,
)
print("RNA missing:", rna_missing)

if after_missing != 0:
    raise RuntimeError(
        "Missing miRNA values remain after imputation."
    )

if age_missing != 0:
    raise RuntimeError(
        "Missing age values remain."
    )

if rna_missing != 0:
    raise RuntimeError(
        "Missing RNA values remain."
    )

# ------------------------------------------------------------------
# Save.
# ------------------------------------------------------------------

OUTPUT.parent.mkdir(
    parents=True,
    exist_ok=True,
)

clean.to_csv(
    OUTPUT,
    index=False,
    compression="gzip",
)

print()
print("=" * 80)
print("SAVED")
print("=" * 80)

print(OUTPUT)

print()
print("FINAL DATASET:")
print(
    f"  {len(clean)} patients"
)
print(
    f"  {len(keep_mirnas)} miRNAs"
)
print(
    f"  {len(rna_cols)} RNA features"
)
print("  age")
print()
print("NO MODEL WAS TRAINED.")
