from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path("data/external/mendeley_ovarian/raw")

raw_file = next(ROOT.rglob("Supplementary data 1.xlsx"))
train_file = next(ROOT.rglob("Supplementary data 4.xlsx"))
test_file = next(ROOT.rglob("Supplementary data 5.xlsx"))

print("=" * 100)
print("MENDELEY OVARIAN DATASET — LABEL / COHORT AUDIT")
print("=" * 100)

def clean_numeric(series):
    """
    Convert numeric-looking strings such as:
      '15.36\\t'
      '>1000'
      '<0.1'
    to numeric values.

    For inequality values, keep the numeric boundary.
    """
    s = series.astype(str).str.strip()

    s = (
        s.str.replace("\t", "", regex=False)
         .str.replace(",", "", regex=False)
    )

    s = s.str.replace(r"^[<>]=?", "", regex=True)

    return pd.to_numeric(s, errors="coerce")


def audit(df, name):
    print("\n" + "=" * 100)
    print(name)
    print("=" * 100)

    print("SHAPE:", df.shape)

    if "TYPE" in df.columns:
        print("\nTYPE COUNTS:")
        print(df["TYPE"].value_counts(dropna=False).sort_index().to_string())

        print("\nTYPE PERCENTAGES:")
        pct = df["TYPE"].value_counts(
            normalize=True,
            dropna=False
        ).sort_index() * 100
        print(pct.round(2).to_string())

        print("\nCROSSTAB TYPE x MENOPAUSE:")
        if "Menopause" in df.columns:
            print(pd.crosstab(df["TYPE"], df["Menopause"], margins=True).to_string())

    print("\nMISSINGNESS:")
    miss = df.isna().sum().sort_values(ascending=False)
    miss_pct = (miss / len(df) * 100).round(2)

    missing_table = pd.DataFrame({
        "missing_n": miss,
        "missing_pct": miss_pct
    })

    print(missing_table.to_string())

    print("\nKEY VARIABLE SUMMARY:")

    key_vars = [
        "Age",
        "Menopause",
        "CA125",
        "HE4",
        "CA19-9",
        "CA72-4",
        "CEA",
        "AFP",
        "NEU",
        "LYM#",
        "PLT",
        "ALB",
    ]

    rows = []

    for col in key_vars:
        if col not in df.columns:
            continue

        x = clean_numeric(df[col])

        rows.append({
            "feature": col,
            "n": len(x),
            "numeric_n": x.notna().sum(),
            "missing_n": x.isna().sum(),
            "missing_pct": round(x.isna().mean() * 100, 2),
            "min": x.min(),
            "median": x.median(),
            "mean": x.mean(),
            "max": x.max(),
        })

    summary = pd.DataFrame(rows)

    if len(summary):
        print(summary.to_string(index=False))

    print("\nSUBJECT ID RANGE:")
    if "SUBJECT_ID" in df.columns:
        print("min:", df["SUBJECT_ID"].min())
        print("max:", df["SUBJECT_ID"].max())
        print("unique:", df["SUBJECT_ID"].nunique())

    if "TYPE" in df.columns:
        print("\nTYPE BY AGE:")
        age_table = (
            df.groupby("TYPE")["Age"]
            .agg(["count", "min", "median", "mean", "max"])
        )
        print(age_table.round(2).to_string())

    print("\n" + "-" * 100)


raw = pd.read_excel(raw_file)
train = pd.read_excel(train_file)
test = pd.read_excel(test_file)

audit(raw, "1. ORIGINAL RAW DATA — 349 PATIENTS")
audit(train, "2. RAW TRAINING DATA — EXPECTED 235 PATIENTS")
audit(test, "3. RAW TEST DATA — EXPECTED 114 PATIENTS")

print("\n" + "=" * 100)
print("TRAIN / TEST SUBJECT OVERLAP")
print("=" * 100)

train_ids = set(train["SUBJECT_ID"])
test_ids = set(test["SUBJECT_ID"])
raw_ids = set(raw["SUBJECT_ID"])

print("Train IDs:", len(train_ids))
print("Test IDs:", len(test_ids))
print("Overlap:", len(train_ids & test_ids))
print("Train missing from raw:", len(train_ids - raw_ids))
print("Test missing from raw:", len(test_ids - raw_ids))
print("Union:", len(train_ids | test_ids))
print("Raw IDs:", len(raw_ids))

print("\n" + "=" * 100)
print("TRAIN + TEST TYPE DISTRIBUTION")
print("=" * 100)

combined = pd.concat(
    [
        train.assign(SPLIT="TRAIN"),
        test.assign(SPLIT="TEST"),
    ],
    ignore_index=True
)

print(
    pd.crosstab(
        combined["SPLIT"],
        combined["TYPE"],
        margins=True
    ).to_string()
)

print("\nTYPE CROSS-TAB WITH SPLIT PERCENTAGES:")

ct = pd.crosstab(
    combined["SPLIT"],
    combined["TYPE"],
    normalize="index"
) * 100

print(ct.round(2).to_string())

print("\n" + "=" * 100)
print("DONE")
print("=" * 100)
