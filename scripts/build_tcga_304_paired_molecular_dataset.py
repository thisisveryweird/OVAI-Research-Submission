#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
import xenaPython as xena


ROOT = Path.home() / "ovarian-ai"

HUB = "https://tcga.xenahubs.net"

RNA_DATASET = (
    "TCGA.OV.sampleMap/HiSeqV2"
)

MIRNA_DATASET = (
    "TCGA.OV.sampleMap/miRNA_HiSeq_gene"
)

CLINICAL_DATASET = (
    "TCGA.OV.sampleMap/OV_clinicalMatrix"
)

MIRNA_AUDIT = (
    ROOT
    / "reports"
    / "tables"
    / "tcga_ov_172mirna_age_audit.csv"
)

OUT = (
    ROOT
    / "data"
    / "processed"
    / "tcga_ov_304_paired_molecular.csv.gz"
)


RNA_GENES = [
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


print("=" * 80)
print("BUILDING TCGA-OV 304-PATIENT PAIRED MOLECULAR DATASET")
print("=" * 80)

# ------------------------------------------------------------------
# Get sample overlap
# ------------------------------------------------------------------

rna_samples = list(map(
    str,
    xena.dataset_samples(
        HUB,
        RNA_DATASET,
        None,
    ),
))

mirna_samples = list(map(
    str,
    xena.dataset_samples(
        HUB,
        MIRNA_DATASET,
        None,
    ),
))

clinical_samples = set(map(
    str,
    xena.dataset_samples(
        HUB,
        CLINICAL_DATASET,
        None,
    ),
))

paired = [
    s for s in rna_samples
    if s in set(mirna_samples)
    and s in clinical_samples
]

print()
print("RNA samples:", len(rna_samples))
print("miRNA samples:", len(mirna_samples))
print("Clinical samples:", len(clinical_samples))
print("Paired samples:", len(paired))

if len(paired) != 304:
    print()
    print(
        "WARNING: expected approximately 304 paired samples."
    )

# ------------------------------------------------------------------
# Load already-audited miRNA data
# ------------------------------------------------------------------

mirna = pd.read_csv(
    MIRNA_AUDIT
)

mirna["sample_id"] = (
    mirna["sample_id"]
    .astype(str)
)

mirna = mirna[
    mirna["sample_id"].isin(paired)
].copy()

print()
print(
    "Paired miRNA rows:",
    len(mirna),
)

# ------------------------------------------------------------------
# Extract RNA genes
# ------------------------------------------------------------------

rna_fields = [
    str(x)
    for x in xena.dataset_field(
        HUB,
        RNA_DATASET,
    )
]

rna_field_set = set(rna_fields)

missing_genes = [
    gene
    for gene in RNA_GENES
    if gene not in rna_field_set
]

if missing_genes:
    raise RuntimeError(
        "Missing RNA genes: "
        + str(missing_genes)
    )

print()
print(
    "RNA genes requested:",
    len(RNA_GENES),
)

print(
    "RNA genes available:",
    len(RNA_GENES)
    - len(missing_genes),
)

# ------------------------------------------------------------------
# Query RNA values
# ------------------------------------------------------------------

print()
print("Querying RNA expression...")

_, values = xena.dataset_probe_values(
    HUB,
    RNA_DATASET,
    paired,
    RNA_GENES,
)

rna_df = pd.DataFrame(
    np.asarray(values).T,
    columns=RNA_GENES,
)

rna_df.insert(
    0,
    "sample_id",
    paired,
)

# ------------------------------------------------------------------
# Clinical data
# ------------------------------------------------------------------

clinical_fields = set(
    str(x)
    for x in xena.dataset_field(
        HUB,
        CLINICAL_DATASET,
    )
)

age_field = "age_at_initial_pathologic_diagnosis"

if age_field not in clinical_fields:
    raise RuntimeError(
        f"Missing clinical field: {age_field}"
    )

_, age_values = xena.dataset_probe_values(
    HUB,
    CLINICAL_DATASET,
    paired,
    [age_field],
)

clinical_df = pd.DataFrame({
    "sample_id": paired,
    "age": pd.to_numeric(
        age_values[0],
        errors="coerce",
    ),
})

# ------------------------------------------------------------------
# Merge
# ------------------------------------------------------------------

df = mirna.merge(
    rna_df,
    on="sample_id",
    how="inner",
    validate="one_to_one",
    suffixes=("", "_rna"),
)

df = df.drop(
    columns=[
        c for c in ["age", "age_rna"]
        if c in df.columns
    ],
    errors="ignore",
)

df = df.merge(
    clinical_df,
    on="sample_id",
    how="inner",
    validate="one_to_one",
)

# The age queried directly from the TCGA clinical matrix
# is the authoritative clinical value.
if "age" not in df.columns:
    raise RuntimeError(
        "Age column missing after clinical merge."
    )

# ------------------------------------------------------------------
# Validation
# ------------------------------------------------------------------

mirna_cols = [
    c for c in df.columns
    if c.startswith("hsa-")
]

print()
print("=" * 80)
print("FINAL DATASET")
print("=" * 80)

print(
    "Patients:",
    len(df),
)

print(
    "miRNA features:",
    len(mirna_cols),
)

print(
    "RNA features:",
    len(RNA_GENES),
)

print(
    "Age observed:",
    int(df["age"].notna().sum()),
)

print(
    "Age missing:",
    int(df["age"].isna().sum()),
)

print()
print("RNA missing values:")

for gene in RNA_GENES:

    print(
        f"  {gene:8} "
        f"{int(df[gene].isna().sum())}"
    )

print()
print("miRNA missing values:", int(
    df[mirna_cols]
    .isna()
    .sum()
    .sum()
))

# ------------------------------------------------------------------
# Save
# ------------------------------------------------------------------

OUT.parent.mkdir(
    parents=True,
    exist_ok=True,
)

df.to_csv(
    OUT,
    index=False,
    compression="gzip",
)

print()
print("=" * 80)
print("SAVED")
print("=" * 80)

print(OUT)

print()
print("IMPORTANT:")
print(
    "This is a paired TCGA molecular dataset."
)
print(
    "No diagnostic model was trained."
)
print(
    "TCGA-OV remains tumor-only and is NOT being"
)
print(
    "used as a cancer-vs-control validation cohort."
)
