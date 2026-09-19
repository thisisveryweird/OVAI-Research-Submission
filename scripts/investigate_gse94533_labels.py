#!/usr/bin/env python3

from pathlib import Path
import gzip
import re
import urllib.request
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INTERIM = ROOT / "data" / "interim"
RAW = ROOT / "data" / "raw"

INTERIM.mkdir(parents=True, exist_ok=True)

SOFT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/"
    "GSE94nnn/GSE94533/soft/GSE94533_family.soft.gz"
)

SOFT = INTERIM / "GSE94533_family.soft.gz"

# ------------------------------------------------------------------
# Download official GEO metadata if necessary
# ------------------------------------------------------------------

if not SOFT.exists():
    print("Downloading official GEO SOFT metadata...")
    req = urllib.request.Request(
        SOFT_URL,
        headers={"User-Agent": "ovarian-ai-research/1.0"},
    )

    with urllib.request.urlopen(req, timeout=120) as r:
        SOFT.write_bytes(r.read())

# ------------------------------------------------------------------
# Parse sample metadata
# ------------------------------------------------------------------

records = []
current = None

with gzip.open(SOFT, "rt", encoding="utf-8", errors="replace") as f:

    for line in f:

        line = line.rstrip("\n")

        if line.startswith("^SAMPLE = "):

            if current is not None:
                records.append(current)

            current = {
                "geo_accession": line.split("=", 1)[1].strip(),
                "title": None,
                "characteristics": [],
            }

            continue

        if current is None:
            continue

        if line.startswith("!Sample_title = "):

            current["title"] = line.split("=", 1)[1].strip()

        elif line.startswith("!Sample_characteristics_ch1 = "):

            value = line.split("=", 1)[1].strip()
            current["characteristics"].append(value)

if current is not None:
    records.append(current)

rows = []

for rec in records:

    row = {
        "geo_accession": rec["geo_accession"],
        "title": rec["title"],
    }

    for characteristic in rec["characteristics"]:

        if ":" in characteristic:
            key, value = characteristic.split(":", 1)
            key = key.strip().lower().replace(" ", "_")
            value = value.strip()
        else:
            key = "other"
            value = characteristic.strip()

        # Keep all characteristics so we don't lose information.
        if key in row:
            row[f"{key}_2"] = value
        else:
            row[key] = value

    rows.append(row)

meta = pd.DataFrame(rows)

print()
print("=" * 100)
print("OFFICIAL GEO METADATA")
print("=" * 100)

print(f"\nNumber of GEO samples: {len(meta)}")

print("\nColumns:")
for c in meta.columns:
    print(" ", repr(c))

print("\nDiagnosis:")
if "diagnosis" in meta.columns:
    print(meta["diagnosis"].value_counts(dropna=False).to_string())
else:
    print("No diagnosis column found.")

# ------------------------------------------------------------------
# Print E-1056 and suspicious records
# ------------------------------------------------------------------

def norm(x):
    if pd.isna(x):
        return ""
    return re.sub(r"[^a-z0-9]", "", str(x).lower())

meta["_norm_title"] = meta["title"].map(norm)

print()
print("=" * 100)
print("E-1056")
print("=" * 100)

e1056 = meta[meta["_norm_title"] == "e1056"]

if len(e1056):
    print(e1056.T.to_string())
else:
    print("E-1056 was not found in the GEO title field.")

# ------------------------------------------------------------------
# Print every diagnosis group with sample IDs
# ------------------------------------------------------------------

print()
print("=" * 100)
print("SAMPLES BY GEO DIAGNOSIS")
print("=" * 100)

if "diagnosis" in meta.columns:

    for diagnosis, group in meta.groupby(
        "diagnosis",
        dropna=False
    ):

        print()
        print(f"### {diagnosis!r} ({len(group)})")

        ids = group["title"].tolist()

        for sample_id in ids:
            print(" ", sample_id)

# ------------------------------------------------------------------
# Inspect both processed workbooks' sample IDs
# ------------------------------------------------------------------

print()
print("=" * 100)
print("PROCESSED WORKBOOK SAMPLE IDS")
print("=" * 100)

all_expression_ids = []

for workbook in sorted(
    RAW.glob("GSE94533_Processed_file_Cohort*.xlsx")
):

    print()
    print(workbook.name)

    try:
        raw = pd.read_excel(
            workbook,
            sheet_name="mirna_tpm",
            header=None,
            nrows=4,
        )
    except Exception as exc:
        print("Could not read:", exc)
        continue

    print("Header rows:")

    print(raw.to_string(index=False, header=False))

    # Row 2 is the actual sample-header row in these workbooks.
    if len(raw) >= 3:

        header = raw.iloc[2]

        for value in header:

            normalized = norm(value)

            if normalized and normalized != "name":
                if normalized not in {
                    "sequence",
                    "namecount",
                    "identicalmirnas",
                }:
                    all_expression_ids.append(
                        (str(value), normalized, workbook.name)
                    )

# ------------------------------------------------------------------
# Compare expression IDs against GEO IDs
# ------------------------------------------------------------------

geo_ids = set(meta["_norm_title"])

expression_ids = {
    normalized
    for _, normalized, _ in all_expression_ids
}

common = expression_ids & geo_ids
expression_only = expression_ids - geo_ids
geo_only = geo_ids - expression_ids

print()
print("=" * 100)
print("ID CROSS-CHECK")
print("=" * 100)

print(f"\nUnique expression sample IDs: {len(expression_ids)}")
print(f"GEO sample IDs:               {len(geo_ids)}")
print(f"Common IDs:                   {len(common)}")
print(f"Expression-only IDs:          {len(expression_only)}")
print(f"GEO-only IDs:                 {len(geo_only)}")

if expression_only:

    print("\nExpression IDs not found in GEO:")
    for x in sorted(expression_only):
        print(" ", x)

if geo_only:

    print("\nGEO IDs not found in expression files:")
    for x in sorted(geo_only):
        print(" ", x)

# ------------------------------------------------------------------
# Save complete metadata
# ------------------------------------------------------------------

out = INTERIM / "gse94533_geo_metadata_full.csv"

meta.drop(
    columns=["_norm_title"]
).to_csv(
    out,
    index=False,
)

print()
print("=" * 100)
print("COMPLETE")
print("=" * 100)
print(f"\nSaved:")
print(out)

print(
    """
IMPORTANT:

Do not train the model yet.

We need to determine whether:
1. the current GEO diagnosis annotations differ from the labels used
   in the original Elias et al. model, or
2. the published 98/21/45/15 counts refer to a particular analytical
   subset/labeling that is preserved in the authors' supplementary data.

The next step is to compare these exact sample IDs with the original
published supplementary dataset rather than guessing.
"""
)
