from pathlib import Path
import gzip
import re

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path("data/external/gse9899")
MATRIX = ROOT / "GSE9899_series_matrix.txt.gz"
ANNOT = ROOT / "GSE9899_sample_annotations.csv"

GENES_FILE = Path(
    "data/processed/gse8841_frozen_borderline_genes.txt"
)

genes = [
    x.strip()
    for x in GENES_FILE.read_text().splitlines()
    if x.strip()
]

# Frozen transferable genes from GSE12172:
# IFI30 was not present on GPL570.
genes = [g for g in genes if g != "IFI30"]

print("=" * 90)
print("GSE9899 — FROZEN 12-GENE SIGNATURE TEST")
print("=" * 90)
print("Genes:", ", ".join(genes))

ann = pd.read_csv(ANNOT)

# Keep only samples with known clinical diagnosis.
ann = ann[ann["Type"].isin(["LMP", "Malignant"])].copy()

print("\nClinical samples:", len(ann))
print("LMP:", (ann["Type"] == "LMP").sum())
print("Malignant:", (ann["Type"] == "Malignant").sum())

# ------------------------------------------------------------
# Read GEO matrix.
# First determine the sample order from !Sample_title.
# ------------------------------------------------------------

sample_ids = None
data_started = False
rows = []

with gzip.open(
    MATRIX,
    "rt",
    encoding="utf-8",
    errors="replace",
) as f:

    for line in f:

        if line.startswith("!Sample_title"):
            sample_ids = [
                x.strip().strip('"')
                for x in line.rstrip().split("\t")[1:]
                if x.strip()
            ]

        elif line.startswith("!series_matrix_table_begin"):
            data_started = True
            break

    if sample_ids is None:
        raise RuntimeError("Could not find !Sample_title.")

    # Read expression rows from here.
    for line in f:

        if line.startswith("!series_matrix_table_end"):
            break

        parts = line.rstrip().split("\t")

        if len(parts) < 2:
            continue

        probe = parts[0].strip().strip('"')

        # Skip metadata/header-like rows that are not expression data.
        if probe.startswith("!"):
            continue

        values = parts[1:]

        clean_values = []

        for x in values:
            x = x.strip().strip('"')

            try:
                clean_values.append(
                    np.nan if x in {"null", "NaN", ""} else float(x)
                )
            except ValueError:
                clean_values.append(np.nan)

        rows.append([probe] + clean_values)

expr = pd.DataFrame(
    rows,
    columns=["probe"] + sample_ids,
)

print("\nExpression matrix:", expr.shape)

# ------------------------------------------------------------
# Need GPL570 probe -> gene-symbol mapping.
# ------------------------------------------------------------

print("\nLooking for GPL570 annotation...")

platform_candidates = [
    Path("data/external/gse_platforms/GPL570.txt"),
    Path("data/external/gse_platforms/GPL570.annot.gz"),
    Path("data/external/gse_platforms/GPL570.txt.gz"),
]

platform_file = next(
    (p for p in platform_candidates if p.exists()),
    None,
)

if platform_file is None:
    raise FileNotFoundError(
        "GPL570 annotation not found. "
        "Expected it in data/external/gse_platforms/"
    )

print("Using:", platform_file)

# Parse GPL570 annotation.
if platform_file.suffix == ".gz":
    opener = gzip.open
else:
    opener = open

with opener(
    platform_file,
    "rt",
    encoding="utf-8",
    errors="replace",
) as f:

    platform_lines = f.readlines()

header_idx = None

for i, line in enumerate(platform_lines):
    if line.startswith("!platform_table_begin"):
        header_idx = i + 1
        break

if header_idx is None:
    raise RuntimeError("Could not find platform table.")

header = platform_lines[header_idx].rstrip().split("\t")

# Find gene-symbol column.
symbol_candidates = [
    i for i, h in enumerate(header)
    if any(
        x in h.upper()
        for x in [
            "GENE_SYMBOL",
            "GENE SYMBOL",
            "SYMBOL",
        ]
    )
]

if not symbol_candidates:
    raise RuntimeError(
        f"Could not locate gene-symbol column. Headers: {header}"
    )

symbol_idx = symbol_candidates[0]

annotation_rows = []

for line in platform_lines[header_idx + 1:]:

    if line.startswith("!platform_table_end"):
        break

    parts = line.rstrip().split("\t")

    if len(parts) <= max(0, symbol_idx):
        continue

    annotation_rows.append(
        (
            parts[0].strip().strip('"'),
            parts[symbol_idx].strip().strip('"'),
        )
    )

probe_to_gene = pd.DataFrame(
    annotation_rows,
    columns=["probe", "gene"],
)

# Clean symbols.
probe_to_gene["gene"] = (
    probe_to_gene["gene"]
    .astype(str)
    .str.strip()
    .str.replace('"', "", regex=False)
)

# Keep only our signature genes.
probe_to_gene = probe_to_gene[
    probe_to_gene["gene"].isin(genes)
].copy()

print(
    "\nSignature genes represented:",
    sorted(probe_to_gene["gene"].unique())
)

missing_genes = sorted(
    set(genes) - set(probe_to_gene["gene"])
)

print("Missing signature genes:", missing_genes)

# ------------------------------------------------------------
# Convert probe expression to gene expression.
# Multiple probes -> median.
# ------------------------------------------------------------

expr_sub = expr[
    expr["probe"].isin(probe_to_gene["probe"])
].copy()

expr_long = expr_sub.merge(
    probe_to_gene,
    on="probe",
    how="left",
)

gene_expr = expr_long.groupby("gene")[sample_ids].median()

print("\nGene expression shape:", gene_expr.shape)

# ------------------------------------------------------------
# Build signature score.
# Standardize each gene within this cohort ONLY for a ranking score.
# This does not alter gene selection.
# ------------------------------------------------------------

available = [
    g for g in genes
    if g in gene_expr.index
]

X = gene_expr.loc[available].T

# z-score each gene across GSE9899.
X = (X - X.mean()) / X.std(ddof=0)

score = X.mean(axis=1)

score_df = pd.DataFrame({
    "AOCSID": score.index,
    "signature_score": score.values,
})

merged = ann.merge(
    score_df,
    on="AOCSID",
    how="inner",
)

print("\nScored samples:", len(merged))

# LMP = positive class for this analysis.
merged["y"] = (merged["Type"] == "LMP").astype(int)

def auc_for(data, name):
    if data["y"].nunique() < 2:
        print(f"{name}: insufficient class variation")
        return

    auc = roc_auc_score(
        data["y"],
        data["signature_score"],
    )

    print(
        f"{name}: "
        f"n={len(data)}, "
        f"LMP={int(data['y'].sum())}, "
        f"malignant={int((data['y']==0).sum())}, "
        f"AUC={auc:.6f}"
    )

print("\n" + "=" * 90)
print("AUC RESULTS")
print("=" * 90)

auc_for(
    merged,
    "LMP vs all malignant",
)

auc_for(
    merged[
        merged["StageCode"].isin(["IA", "IB", "IC"]) |
        (merged["Type"] == "LMP")
    ],
    "LMP vs Stage-I malignant",
)

auc_for(
    merged[
        ~merged["StageCode"].isin(["IA", "IB", "IC"]) |
        (merged["Type"] == "LMP")
    ],
    "LMP vs non-Stage-I malignant",
)

# Save scores.
out = ROOT / "GSE9899_frozen_signature_scores.csv"
merged.to_csv(out, index=False)

print("\nSaved:", out)
print("DONE")
