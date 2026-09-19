from pathlib import Path
import pandas as pd

ROOT = Path.home() / "ovarian-ai"

BENCHMARK = ROOT / "data/processed/gse94533_14mirna_clinical_benchmark.csv.gz"
PMML = ROOT / "reports/tables/original_historical_clean_predictions.csv"
MAPPING = ROOT / "reports/tables/original_historical_clean_sample_mapping.csv"

OUT = ROOT / "data/processed/gse94533_14mirna_clinical_benchmark_fixed.csv.gz"

print("=" * 80)
print("FIXING VERIFIED CLINICAL BENCHMARK")
print("=" * 80)

benchmark = pd.read_csv(BENCHMARK)
pmml = pd.read_csv(PMML)
mapping = pd.read_csv(MAPPING)

verified_ids = pmml["sample_id"].astype(str).tolist()
verified_set = set(verified_ids)

print()
print("Current benchmark:", len(benchmark))
print("Verified PMML cohort:", len(verified_ids))

# ---------------------------------------------------------------
# Remove the incorrectly included excluded patient.
# ---------------------------------------------------------------

benchmark = benchmark[
    benchmark["sample_id"].astype(str).isin(verified_set)
].copy()

# ---------------------------------------------------------------
# Check whether e1056 is available anywhere in the mapping/
# benchmark sources.
# ---------------------------------------------------------------

has_e1056 = (
    benchmark["sample_id"]
    .astype(str)
    .str.lower()
    .eq("e1056")
    .any()
)

print()
print("e1056 already present:", has_e1056)
print("81195 present after filtering:",
      "81195" in set(benchmark["sample_id"].astype(str)))

# ---------------------------------------------------------------
# If e1056 is missing, recover it from Dataset 2 using the
# existing verified PMML prediction file's feature values.
#
# The historical clean reproduction already produced the exact
# PMML feature matrix in its forward-trace/prediction ecosystem.
# We use the mapping file to locate e1056's Dataset-2 column.
# ---------------------------------------------------------------

if not has_e1056:

    raise RuntimeError(
        "e1056 is missing after filtering. "
        "Do not train models yet. "
        "We need to recover e1056's 14 miRNA values from Dataset 2."
    )

# ---------------------------------------------------------------
# Final exact cohort check.
# ---------------------------------------------------------------

final_ids = set(
    benchmark["sample_id"].astype(str)
)

if final_ids != verified_set:
    print()
    print("Missing:", sorted(verified_set - final_ids))
    print("Extra:", sorted(final_ids - verified_set))
    raise RuntimeError(
        "Clinical benchmark does not exactly match verified cohort."
    )

if "81195" in final_ids:
    raise RuntimeError("81195 is still present.")

if "e1056" not in final_ids:
    raise RuntimeError("e1056 is still missing.")

benchmark.to_csv(
    OUT,
    index=False,
    compression="gzip",
)

print()
print("=" * 80)
print("SUCCESS")
print("=" * 80)
print("Patients:", len(benchmark))
print("81195 present:", "81195" in final_ids)
print("e1056 present:", "e1056" in final_ids)
print("Saved:", OUT)
