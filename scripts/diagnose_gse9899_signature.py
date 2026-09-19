from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

ROOT = Path("data/external/gse9899")

df = pd.read_csv(
    ROOT / "GSE9899_frozen_signature_scores.csv"
)

genes = [
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

print("=" * 90)
print("GSE9899 FROZEN SIGNATURE DIAGNOSTIC")
print("=" * 90)

print("\nOverall score:")
print(
    df.groupby("Type")["signature_score"]
      .agg(["count", "mean", "median", "std"])
      .to_string()
)

print("\nStage-I comparison:")
stage1 = df[
    df["StageCode"].isin(["IA", "IB", "IC"]) |
    (df["Type"] == "LMP")
].copy()

print(
    stage1.groupby("Type")["signature_score"]
          .agg(["count", "mean", "median", "std"])
          .to_string()
)

print("\n" + "-" * 90)
print("GENE-LEVEL LMP VS MALIGNANT")
print("-" * 90)

# Reconstruct gene values is not present in the score file.
# Instead, inspect whether each gene is available in the saved file.
present = [g for g in genes if g in df.columns]

print("Gene columns available in score file:", present)

if not present:
    print(
        "\nThe score file only contains the combined score."
        "\nWe need to regenerate the individual gene values."
    )

print("\n" + "=" * 90)
print("DONE")
print("=" * 90)
