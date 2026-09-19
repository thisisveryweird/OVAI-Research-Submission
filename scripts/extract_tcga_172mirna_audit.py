from pathlib import Path
import pandas as pd
import numpy as np
import xenaPython as xena

ROOT = Path.home() / "ovarian-ai"

HUB = "https://tcga.xenahubs.net"

MIRNA_DATASET = (
    "TCGA.OV.sampleMap/miRNA_HiSeq_gene"
)

CLINICAL_DATASET = (
    "TCGA.OV.sampleMap/OV_clinicalMatrix"
)

MAPPING = (
    ROOT
    / "reports"
    / "tables"
    / "historical_mirna_to_mimat.csv"
)

OUT = (
    ROOT
    / "reports"
    / "tables"
    / "tcga_ov_172mirna_age_audit.csv"
)


print("=" * 80)
print("TCGA-OV 172-miRNA + AGE AUDIT")
print("=" * 80)

mapping = pd.read_csv(MAPPING)

mapping = mapping[
    mapping["matched"] == True
].copy()

if len(mapping) != 172:
    raise RuntimeError(
        f"Expected 172 matched miRNAs, found {len(mapping)}."
    )

mimat_ids = (
    mapping["mimat"]
    .astype(str)
    .tolist()
)

name_by_mimat = dict(
    zip(
        mapping["mimat"].astype(str),
        mapping["historical_name"].astype(str),
    )
)

print()
print("Mapped historical miRNAs:", len(mapping))

# -----------------------------------------------------------------
# TCGA miRNA samples
# -----------------------------------------------------------------

mirna_samples = xena.dataset_samples(
    HUB,
    MIRNA_DATASET,
    None,
)

print(
    "TCGA miRNA samples:",
    len(mirna_samples),
)

# -----------------------------------------------------------------
# Determine which requested MIMAT IDs are actually in Xena.
# -----------------------------------------------------------------

tcga_fields = set(
    str(x)
    for x in xena.dataset_field(
        HUB,
        MIRNA_DATASET,
    )
)

usable_mimat = [
    x for x in mimat_ids
    if x in tcga_fields
]

missing_tcga = [
    x for x in mimat_ids
    if x not in tcga_fields
]

print(
    "Requested MIMAT IDs:",
    len(mimat_ids),
)

print(
    "Present in TCGA dataset:",
    len(usable_mimat),
)

print(
    "Absent from TCGA dataset:",
    len(missing_tcga),
)

if missing_tcga:
    print()
    print("First missing MIMAT IDs:")
    for x in missing_tcga[:20]:
        print(" ", x)

# -----------------------------------------------------------------
# Query all usable miRNA values.
# -----------------------------------------------------------------

print()
print("Querying TCGA miRNA values...")

position, values = xena.dataset_probe_values(
    HUB,
    MIRNA_DATASET,
    mirna_samples,
    usable_mimat,
)

if len(values) != len(usable_mimat):
    raise RuntimeError(
        "Unexpected number of miRNA result rows."
    )

mirna_df = pd.DataFrame(
    np.asarray(values).T,
    columns=[
        name_by_mimat[x]
        for x in usable_mimat
    ],
)

mirna_df.insert(
    0,
    "sample_id",
    mirna_samples,
)

# -----------------------------------------------------------------
# Clinical samples
# -----------------------------------------------------------------

clinical_samples = xena.dataset_samples(
    HUB,
    CLINICAL_DATASET,
    None,
)

print(
    "TCGA clinical samples:",
    len(clinical_samples),
)

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
        f"Could not find {age_field}."
    )

clinical_probe_samples = [
    x for x in mirna_samples
    if x in set(clinical_samples)
]

print(
    "Samples with both miRNA + clinical records:",
    len(clinical_probe_samples),
)

_, age_values = xena.dataset_probe_values(
    HUB,
    CLINICAL_DATASET,
    clinical_probe_samples,
    [age_field],
)

age_df = pd.DataFrame({
    "sample_id": clinical_probe_samples,
    "age": pd.to_numeric(
        age_values[0],
        errors="coerce",
    ),
})

# -----------------------------------------------------------------
# Merge
# -----------------------------------------------------------------

df = mirna_df.merge(
    age_df,
    on="sample_id",
    how="inner",
    validate="one_to_one",
)

print()
print("=" * 80)
print("FINAL TCGA AUDIT DATASET")
print("=" * 80)

print(
    "Patients with miRNA + age:",
    len(df),
)

print(
    "MiRNA features:",
    len(usable_mimat),
)

print(
    "Age observed:",
    int(df["age"].notna().sum()),
)

print(
    "Age missing:",
    int(df["age"].isna().sum()),
)

# -----------------------------------------------------------------
# Feature completeness
# -----------------------------------------------------------------

missing_counts = (
    df[
        [
            name_by_mimat[x]
            for x in usable_mimat
        ]
    ]
    .isna()
    .sum()
    .sort_values(
        ascending=False
    )
)

print()
print("Features with most missing values:")

for name, count in missing_counts.head(10).items():
    print(
        f"  {name}: {int(count)}"
    )

# -----------------------------------------------------------------
# Save
# -----------------------------------------------------------------

df.to_csv(
    OUT,
    index=False,
)

print()
print("=" * 80)
print("SAVED")
print("=" * 80)
print(OUT)

print()
print("IMPORTANT:")
print("This is an AUDIT dataset.")
print("No model was trained.")
print("TCGA has NOT been used as a cancer-vs-noncancer test.")
print("Measurement-scale compatibility must be checked before")
print("any external prediction is attempted.")
