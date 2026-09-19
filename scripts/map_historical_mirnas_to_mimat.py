from pathlib import Path
import gzip
import re
import urllib.request
import pandas as pd

ROOT = Path.home() / "ovarian-ai"

XLSX = (
    ROOT
    / "data"
    / "external"
    / "elife_28932"
    / "elife_poa_e28932_Supplementary_File_6.xlsx"
)

OUT = (
    ROOT
    / "reports"
    / "tables"
    / "historical_mirna_to_mimat.csv"
)

FASTA = ROOT / "data" / "external" / "mirbase_mature.fa"

URL = "https://www.mirbase.org/download/CURRENT/mature.fa"


def norm_name(x):
    s = str(x).strip().lower()
    s = s.replace("mirna", "mir")
    return re.sub(r"[^a-z0-9-]", "", s)


print("=" * 80)
print("MAPPING HISTORICAL miRNA NAMES TO miRBase MIMAT IDs")
print("=" * 80)

# ---------------------------------------------------------------
# Get current miRBase mature annotation.
# ---------------------------------------------------------------

FASTA.parent.mkdir(
    parents=True,
    exist_ok=True,
)

if not FASTA.exists():

    print()
    print("Downloading miRBase mature.fa.gz...")

    urllib.request.urlretrieve(
        URL,
        FASTA,
    )

# ---------------------------------------------------------------
# Parse human mature miRNAs.
#
# Typical header:
# >MIMAT0000418 hsa-miR-23b-3p ...
# ---------------------------------------------------------------

mapping = {}

with open(
    FASTA,
    "rt",
    encoding="utf-8",
) as f:

    for line in f:

        if not line.startswith(">"):
            continue

        parts = line[1:].strip().split()

        if len(parts) < 2:
            continue

        name = parts[0]
        accession = parts[1]

        if not name.lower().startswith("hsa-"):
            continue

        mapping[norm_name(name)] = accession

# ---------------------------------------------------------------
# Extract the 192 historical miRNAs.
# ---------------------------------------------------------------

raw = pd.read_excel(
    XLSX,
    sheet_name="Supplemental Dataset 2",
    header=None,
    engine="openpyxl",
)

historical = list(dict.fromkeys(
    str(x).strip()
    for x in raw.iloc[:, 0]
    if str(x).strip().startswith("hsa-")
))

print()
print("Historical miRNAs:", len(historical))
print("Human miRBase names:", len(mapping))

records = []

for name in historical:

    key = norm_name(name)
    mimat = mapping.get(key)

    records.append({
        "historical_name": name,
        "mimat": mimat,
        "matched": mimat is not None,
    })

result = pd.DataFrame(records)

result.to_csv(
    OUT,
    index=False,
)

matched = int(result["matched"].sum())
missing = len(result) - matched

print()
print("=" * 80)
print("MAPPING RESULT")
print("=" * 80)

print("Historical miRNAs:", len(result))
print("Matched to MIMAT:", matched)
print("Not matched:", missing)

if missing:

    print()
    print("UNMATCHED miRNAs:")

    for x in result.loc[
        ~result["matched"],
        "historical_name",
    ]:
        print(" ", x)

print()
print("EXAMPLE MATCHES:")

print(
    result[
        result["matched"]
    ].head(20).to_string(index=False)
)

print()
print("Saved:")
print(OUT)
