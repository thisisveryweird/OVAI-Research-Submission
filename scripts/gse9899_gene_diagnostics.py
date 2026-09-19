from pathlib import Path
import gzip
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path("data/external/gse9899")

MATRIX = ROOT / "GSE9899_series_matrix.txt.gz"
ANNOT = ROOT / "GSE9899_sample_annotations.csv"

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

ann = pd.read_csv(ANNOT)
ann = ann[ann["Type"].isin(["LMP", "Malignant"])].copy()

# ------------------------------------------------------------
# Read expression matrix
# ------------------------------------------------------------

with gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace") as f:

    for line in f:
        if line.startswith("!Sample_title"):
            sample_ids = [
                x.strip().strip('"')
                for x in line.rstrip().split("\t")[1:]
                if x.strip()
            ]
            break

    for line in f:
        if line.startswith("!series_matrix_table_begin"):
            break

    rows = []

    for line in f:
        if line.startswith("!series_matrix_table_end"):
            break

        parts = line.rstrip().split("\t")

        if len(parts) < 2:
            continue

        probe = parts[0].strip().strip('"')

        if probe.startswith("!"):
            continue

        vals = []

        for x in parts[1:]:
            x = x.strip().strip('"')

            try:
                vals.append(float(x))
            except ValueError:
                vals.append(np.nan)

        rows.append([probe] + vals)

expr = pd.DataFrame(rows, columns=["probe"] + sample_ids)

# ------------------------------------------------------------
# GPL570 annotation
# ------------------------------------------------------------

platform = Path("data/external/gse_platforms/GPL570.txt")

lines = platform.read_text(
    encoding="utf-8",
    errors="replace"
).splitlines()

for i, line in enumerate(lines):
    if line.startswith("!platform_table_begin"):
        header_idx = i + 1
        break
else:
    raise RuntimeError("GPL570 table not found.")

header = lines[header_idx].split("\t")

symbol_candidates = [
    i for i, h in enumerate(header)
    if any(
        x in h.upper()
        for x in ["GENE_SYMBOL", "GENE SYMBOL", "SYMBOL"]
    )
]

if not symbol_candidates:
    raise RuntimeError("Gene-symbol column not found.")

symbol_idx = symbol_candidates[0]

mapping = []

for line in lines[header_idx + 1:]:
    if line.startswith("!platform_table_end"):
        break

    parts = line.split("\t")

    if len(parts) <= symbol_idx:
        continue

    probe = parts[0].strip().strip('"')
    gene = parts[symbol_idx].strip().strip('"')

    if gene in GENES:
        mapping.append((probe, gene))

mapping = pd.DataFrame(
    mapping,
    columns=["probe", "gene"]
)

# ------------------------------------------------------------
# Probe -> gene median
# ------------------------------------------------------------

sub = expr[
    expr["probe"].isin(mapping["probe"])
].merge(
    mapping,
    on="probe",
    how="left"
)

gene_expr = sub.groupby("gene")[sample_ids].median()

X = gene_expr.T

# Standardize within GSE9899 only, purely for comparability.
X = (X - X.mean()) / X.std(ddof=0)

# Attach diagnosis
data = ann.merge(
    pd.DataFrame({"AOCSID": X.index}),
    on="AOCSID",
    how="inner"
)

# Sample order should be identical after merge.
X = X.loc[data["AOCSID"]]

y = (data["Type"] == "LMP").astype(int).to_numpy()

print("=" * 100)
print("GSE9899 INDIVIDUAL FROZEN-GENE DIAGNOSTICS")
print("=" * 100)

print("Samples:", len(data))
print("LMP:", int(y.sum()))
print("Malignant:", int((y == 0).sum()))

results = []

for gene in GENES:

    if gene not in X.columns:
        print(gene, "MISSING")
        continue

    score = X[gene].to_numpy()

    auc = roc_auc_score(y, score)
    reversed_auc = 1 - auc
    best_direction_auc = max(auc, reversed_auc)

    results.append({
        "gene": gene,
        "auc": auc,
        "reversed_auc": reversed_auc,
        "best_direction_auc": best_direction_auc,
        "lmp_mean": score[y == 1].mean(),
        "malignant_mean": score[y == 0].mean(),
        "direction": "higher_LMP" if auc >= 0.5 else "higher_malignant",
    })

results_df = pd.DataFrame(results).sort_values(
    "best_direction_auc",
    ascending=False
)

print("\n" + "-" * 100)
print("ALL-SAMPLE RESULTS")
print("-" * 100)
print(results_df.to_string(index=False))

# ------------------------------------------------------------
# Stage-I malignant diagnostic
# ------------------------------------------------------------

stage1 = (
    (data["Type"] == "LMP") |
    (
        (data["Type"] == "Malignant") &
        data["StageCode"].isin(["IA", "IB", "IC"])
    )
)

X1 = X.loc[stage1.values]
data1 = data.loc[stage1].copy()
y1 = (data1["Type"] == "LMP").astype(int).to_numpy()

print("\n" + "-" * 100)
print("STAGE-I MALIGNANT VS LMP")
print("-" * 100)

print("Samples:", len(y1))
print("LMP:", int(y1.sum()))
print("Stage-I malignant:", int((y1 == 0).sum()))

stage_results = []

for gene in GENES:

    if gene not in X1.columns:
        continue

    score = X1[gene].to_numpy()

    auc = roc_auc_score(y1, score)

    stage_results.append({
        "gene": gene,
        "auc": auc,
        "reversed_auc": 1 - auc,
        "best_direction_auc": max(auc, 1 - auc),
    })

stage_df = pd.DataFrame(stage_results).sort_values(
    "best_direction_auc",
    ascending=False
)

print(stage_df.to_string(index=False))

# ------------------------------------------------------------
# Frozen combined score
# ------------------------------------------------------------

combined = X[GENES].mean(axis=1)

overall_auc = roc_auc_score(y, combined)

combined_stage = combined.loc[stage1.values]
stage_auc = roc_auc_score(y1, combined_stage)

print("\n" + "=" * 100)
print("COMBINED SCORE CHECK")
print("=" * 100)

print("All malignant vs LMP AUC:", overall_auc)
print("Stage-I malignant vs LMP AUC:", stage_auc)

out = ROOT / "GSE9899_gene_level_diagnostics.csv"
results_df.to_csv(out, index=False)

stage_out = ROOT / "GSE9899_stage1_gene_diagnostics.csv"
stage_df.to_csv(stage_out, index=False)

print("\nSaved:")
print(out)
print(stage_out)

print("\nDONE")
