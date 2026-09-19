from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path("data/external/gse9899")

clinical_file = ROOT / "GSE9899_clinical_anns.csv"
matrix_file = ROOT / "GSE9899_series_matrix.txt.gz"

print("=" * 100)
print("GSE9899 CLINICAL / LABEL AUDIT")
print("=" * 100)

# ---------------------------------------------------------------------
# 1. Clinical annotations
# ---------------------------------------------------------------------

clinical = pd.read_csv(clinical_file)

print("\nCLINICAL FILE")
print("-" * 100)
print("Shape:", clinical.shape)

print("\nColumns:")
for i, col in enumerate(clinical.columns):
    print(f"[{i}] {col!r}")

print("\nFirst 10 rows:")
print(clinical.head(10).to_string())

print("\nData types:")
print(clinical.dtypes.to_string())

# Every categorical-looking column with a manageable number of values
print("\n" + "=" * 100)
print("CATEGORICAL / LOW-CARDINALITY VARIABLES")
print("=" * 100)

for col in clinical.columns:

    s = clinical[col]

    nunique = s.nunique(dropna=True)

    if (
        s.dtype == "object"
        or nunique <= 20
    ):
        print("\n" + "-" * 80)
        print(f"{col!r} | unique={nunique}")
        print("-" * 80)

        print(
            s.astype(str)
             .replace("nan", np.nan)
             .value_counts(dropna=False)
             .head(50)
             .to_string()
        )

# ---------------------------------------------------------------------
# 2. Search column names for key concepts
# ---------------------------------------------------------------------

print("\n" + "=" * 100)
print("POTENTIALLY IMPORTANT CLINICAL COLUMNS")
print("=" * 100)

keywords = [
    "diagn",
    "status",
    "type",
    "tumor",
    "tumour",
    "border",
    "lmp",
    "low",
    "malig",
    "cancer",
    "stage",
    "figo",
    "grade",
    "hist",
    "mucin",
    "serous",
    "endometr",
    "clear",
    "age",
]

for col in clinical.columns:
    text = str(col).lower()

    if any(k in text for k in keywords):
        print(col)

# ---------------------------------------------------------------------
# 3. Matrix header / sample IDs
# ---------------------------------------------------------------------

print("\n" + "=" * 100)
print("GEO SERIES MATRIX SAMPLE AUDIT")
print("=" * 100)

sample_ids = []

try:
    with open(matrix_file, "rt", encoding="utf-8", errors="replace") as f:

        for line in f:

            if line.startswith("!Sample_geo_accession"):
                print("\nFOUND GEO SAMPLE ACCESSIONS:")
                print(line[:5000])
                sample_ids = (
                    line.strip()
                    .split("\t")[1:]
                )
                break

            if line.startswith("!Sample_title"):
                print("\nFOUND SAMPLE TITLES:")
                print(line[:5000])

            if line.startswith("!Sample_characteristics_ch1"):
                print("\nFOUND SAMPLE CHARACTERISTICS:")
                print(line[:5000])

except Exception as e:
    print("Matrix header error:", repr(e))

# ---------------------------------------------------------------------
# 4. Match clinical sample IDs against GEO IDs
# ---------------------------------------------------------------------

print("\n" + "=" * 100)
print("SAMPLE ID MATCHING")
print("=" * 100)

if sample_ids:

    sample_ids = [
        x.strip().strip('"')
        for x in sample_ids
        if x.strip()
    ]

    print("GEO matrix samples:", len(sample_ids))

    # Try every clinical column to find exact overlap with GEO IDs.
    for col in clinical.columns:

        vals = set(
            clinical[col]
            .dropna()
            .astype(str)
            .str.strip()
        )

        overlap = set(sample_ids) & vals

        if overlap:
            print(
                f"\nColumn {col!r}: "
                f"{len(overlap)} exact GEO-ID matches"
            )
            print("Examples:", list(sorted(overlap))[:10])

# ---------------------------------------------------------------------
# 5. Save a basic audited copy
# ---------------------------------------------------------------------

out = ROOT / "GSE9899_clinical_audit.csv"
clinical.to_csv(out, index=False)

print("\nSaved:", out)

print("\n" + "=" * 100)
print("DONE")
print("=" * 100)
