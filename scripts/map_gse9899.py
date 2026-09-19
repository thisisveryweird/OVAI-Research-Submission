from pathlib import Path
import gzip
import re
import pandas as pd

ROOT = Path("data/external/gse9899")
matrix = ROOT / "GSE9899_series_matrix.txt.gz"
clinical = pd.read_csv(ROOT / "GSE9899_clinical_anns.csv")

# Get sample titles from GEO matrix
with gzip.open(matrix, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("!Sample_title"):
            titles = [
                x.strip().strip('"')
                for x in line.rstrip().split("\t")[1:]
                if x.strip()
            ]
            break

mapped = pd.DataFrame({"AOCSID": titles}).merge(
    clinical,
    on="AOCSID",
    how="left",
    indicator=True
)

print("=" * 80)
print("GSE9899 SAMPLE MAPPING")
print("=" * 80)

print("Expression samples:", len(titles))
print("Clinical patients:", len(clinical))
print("Matched:", (mapped["_merge"] == "both").sum())
print("Unmatched expression:", (mapped["_merge"] == "left_only").sum())
print("Clinical without expression:", (mapped["_merge"] == "right_only").sum())

print("\nTYPE:")
print(mapped["Type"].value_counts(dropna=False).to_string())

print("\nTYPE × STAGE:")
print(
    pd.crosstab(
        mapped["Type"],
        mapped["StageCode"],
        dropna=False
    ).to_string()
)

print("\nTYPE × SUBTYPE:")
print(
    pd.crosstab(
        mapped["Type"],
        mapped["Subtype"],
        dropna=False
    ).to_string()
)

print("\nLMP SAMPLES:")
print(
    mapped.loc[
        mapped["Type"].eq("LMP"),
        ["AOCSID", "Subtype", "StageCode", "Consolidated.Grade"]
    ].to_string(index=False)
)

print("\nSTAGE I MALIGNANT:")
stage1 = ["IA", "IB", "IC"]
print(
    mapped.loc[
        mapped["Type"].eq("Malignant") &
        mapped["StageCode"].isin(stage1),
        ["AOCSID", "Subtype", "StageCode", "Consolidated.Grade"]
    ].to_string(index=False)
)

mapped = mapped.drop(columns="_merge")

out = ROOT / "GSE9899_sample_annotations.csv"
mapped.to_csv(out, index=False)

print("\nSaved:", out)
print("\nDONE")
