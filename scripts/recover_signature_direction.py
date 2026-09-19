from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path("data/external/gse8841")

files = [
    Path("data/processed/gse8841_borderline_expression.csv.gz"),
    Path("data/processed/gse8841_frozen_borderline_genes.txt"),
]

for f in files:
    print(f, "exists =", f.exists())

expr_file = files[0]
genes_file = files[1]

df = pd.read_csv(expr_file)

print("\nColumns:", len(df.columns))
print("Shape:", df.shape)

genes = [
    x.strip()
    for x in genes_file.read_text().splitlines()
    if x.strip()
]

print("\nFrozen genes:")
print(", ".join(genes))

print("\nPossible label columns:")
for c in df.columns:
    if any(
        k in c.lower()
        for k in [
            "label",
            "diagn",
            "type",
            "status",
            "stage",
            "border",
            "cancer",
            "tumor",
        ]
    ):
        print(" ", c)

print("\nFirst rows:")
print(df.head().to_string())

