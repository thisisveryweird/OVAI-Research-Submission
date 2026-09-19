from pathlib import Path
import gzip
import csv
import re
import pandas as pd

ROOT = Path.home() / "ovarian-ai"

MATRIX = (
    ROOT
    / "data"
    / "external"
    / "gse106817"
    / "GSE106817_series_matrix.txt.gz"
)

MAPPING = (
    ROOT
    / "reports"
    / "tables"
    / "historical_mirna_to_mimat.csv"
)


def parse_metadata_line(line):
    parts = next(
        csv.reader(
            [line.rstrip("\n")],
            delimiter="\t",
            quotechar='"',
        )
    )
    return parts[0], parts[1:]


print("=" * 80)
print("GSE106817 EXTERNAL VALIDATION AUDIT")
print("=" * 80)

metadata = {}

with gzip.open(
    MATRIX,
    "rt",
    encoding="utf-8",
    errors="replace",
) as f:

    in_table = False
    feature_ids = []

    for line in f:

        if line.startswith("!series_matrix_table_begin"):
            in_table = True
            continue

        if line.startswith("!series_matrix_table_end"):
            break

        if not in_table and line.startswith("!Sample_"):
            key, values = parse_metadata_line(line)
            metadata[key] = values
            continue

        if in_table and not line.startswith("#"):
            first = line.split("\t", 1)[0].strip().strip('"')

            if first:
                feature_ids.append(first)


# ------------------------------------------------------------------
# Basic metadata
# ------------------------------------------------------------------

accessions = metadata.get(
    "!Sample_geo_accession",
    []
)

titles = metadata.get(
    "!Sample_title",
    []
)

characteristics = metadata.get(
    "!Sample_characteristics_ch1",
    []
)

print()
print("Samples:", len(accessions))
print("Sample titles:", len(titles))
print("Feature rows:", len(feature_ids))


# ------------------------------------------------------------------
# Sample title group audit
# ------------------------------------------------------------------

print()
print("=" * 80)
print("SAMPLE GROUP COUNTS")
print("=" * 80)

patterns = {
    "Ovarian cancer": r"^Ovarian Cancer\b",
    "Borderline ovarian": r"^Borderline\b.*Ovarian|^Ovarian.*Borderline",
    "Benign ovarian": r"^Benign\b.*Ovarian|^Ovarian.*Benign",
    "Non-cancer": r"^Non[- ]cancer\b",
    "Other solid cancer": r"Solid Cancer|Cancer\b",
}

for label, pattern in patterns.items():

    count = sum(
        bool(re.search(
            pattern,
            title,
            flags=re.I,
        ))
        for title in titles
    )

    print(f"{label:25} {count}")


# ------------------------------------------------------------------
# Print unique useful title prefixes
# ------------------------------------------------------------------

print()
print("=" * 80)
print("TITLE PREFIX COUNTS")
print("=" * 80)

prefixes = {}

for title in titles:

    # Remove sample ID at the end.
    cleaned = re.sub(
        r"\s+[A-Z]{1,5}\d+$",
        "",
        title.strip(),
    )

    # First three words.
    prefix = " ".join(
        cleaned.split()[:3]
    )

    prefixes[prefix] = (
        prefixes.get(prefix, 0) + 1
    )

for prefix, count in sorted(
    prefixes.items(),
    key=lambda x: -x[1],
)[:30]:

    print(
        f"{count:5}  {prefix}"
    )


# ------------------------------------------------------------------
# Age audit
# ------------------------------------------------------------------

print()
print("=" * 80)
print("AGE METADATA AUDIT")
print("=" * 80)

age_hits = []

for key, values in metadata.items():

    joined = " ".join(
        str(x)
        for x in values
    )

    if re.search(
        r"\bage\b|\byears?\b",
        joined,
        flags=re.I,
    ):
        age_hits.append(key)

print(
    "Metadata fields containing age/years:",
    age_hits if age_hits else "NONE",
)


# ------------------------------------------------------------------
# MIMAT overlap
# ------------------------------------------------------------------

mapping = pd.read_csv(
    MAPPING
)

mapped = set(
    mapping.loc[
        mapping["matched"] == True,
        "mimat",
    ]
    .astype(str)
)

feature_set = set(
    feature_ids
)

present = sorted(
    mapped & feature_set
)

missing = sorted(
    mapped - feature_set
)

print()
print("=" * 80)
print("MIMAT FEATURE OVERLAP")
print("=" * 80)

print(
    "Historical MIMAT IDs:",
    len(mapped),
)

print(
    "Present in GSE106817:",
    len(present),
)

print(
    "Missing from GSE106817:",
    len(missing),
)

if missing:
    print()
    print("Missing MIMAT IDs:")
    for x in missing:
        print(" ", x)

print()
print("=" * 80)
print("AUDIT COMPLETE")
print("=" * 80)
