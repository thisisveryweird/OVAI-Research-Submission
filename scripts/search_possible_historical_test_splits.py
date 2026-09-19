#!/usr/bin/env python3

"""
Search for stratified 44-patient test sets that are compatible with
the published eLife GSE94533 demographic counts and reported NN AUC.

IMPORTANT
---------
This does NOT claim to recover the original random seed.

The published paper tells us:

    total = 179
    training = 135
    testing = 44

Testing-set pathology counts:

    Control                 4
    Benign                 11
    Borderline              5
    Stage I/II cancer      12
    Stage III/IV cancer    12

The PMML network has already been independently reproduced.

We use the reproduced cancer activation as the ranking score and
search many stratified random splits.

A candidate is only a COMPATIBLE split unless independent evidence
later proves it was the original split.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


ROOT = Path.home() / "ovarian-ai"

PREDICTIONS = (
    ROOT
    / "reports"
    / "tables"
    / "original_historical_clean_predictions.csv"
)

OUTPUT = (
    ROOT
    / "reports"
    / "tables"
    / "possible_historical_test_splits.csv"
)

N_TRIALS = 20_000
def pathology_stratum(row):

    diagnosis = str(
        row.get("diagnosis_final", "")
    ).strip().lower()

    stage_value = row.get("stage_final", "")
    stage = str(stage_value).strip().lower()

    # ---------------------------------------------------------------
    # Published Table 2 identifies E-1056 as:
    #
    #   Diagnosis: Cancer
    #   Stage: I
    #
    # Its GEO metadata has no stage value, so use the published
    # clinical annotation for this one historical record.
    # ---------------------------------------------------------------

    candidate_ids = []

    for field in (
        "sample_id_normalized",
        "sample_key",
        "sample_id",
        "title",
    ):
        value = row.get(field)

        if value is not None:
            candidate_ids.append(
                str(value).strip().lower()
            )

    is_e1056 = any(
        value.replace("-", "") == "e1056"
        for value in candidate_ids
    )

    if diagnosis == "cancer":

        if is_e1056 and (
            not stage
            or stage in {"nan", "none"}
        ):
            return "Stage I/II"

        if "stage i/ii" in stage:
            return "Stage I/II"

        if "stage iii/iv" in stage:
            return "Stage III/IV"

        if stage in {"i", "ii", "1", "2"}:
            return "Stage I/II"

        if stage in {"iii", "iv", "3", "4"}:
            return "Stage III/IV"

        raise RuntimeError(
            "Cancer patient has an unrecognized stage: "
            f"diagnosis={diagnosis!r}, "
            f"stage={stage_value!r}, "
            f"ids={candidate_ids!r}"
        )

    if diagnosis == "benign":
        return "Benign"

    if diagnosis == "borderline":
        return "Borderline"

    if diagnosis == "control":
        return "Control"

    raise RuntimeError(
        "Could not classify patient: "
        f"diagnosis={diagnosis!r}, "
        f"stage={stage_value!r}, "
        f"ids={candidate_ids!r}"
    )



print("=" * 90)
print("SEARCHING COMPATIBLE HISTORICAL 44-PATIENT TEST SETS")
print("=" * 90)

df = pd.read_csv(PREDICTIONS)

print()
print(f"Patients loaded: {len(df)}")

if len(df) != 179:
    raise RuntimeError(
        f"Expected 179 patients, found {len(df)}."
    )

# ---------------------------------------------------------------
# Recover the clinical strata.
# ---------------------------------------------------------------

# The prediction file may not currently contain stage.
# Load it from the processed historical sample mapping/metadata.
METADATA = (
    ROOT
    / "data"
    / "processed"
    / "gse94533_metadata.csv"
)

if not METADATA.exists():

    METADATA = (
        ROOT
        / "data"
        / "processed"
        / "gse94533_metadata.csv"
    )

metadata = pd.read_csv(METADATA)

print()
print("Metadata columns:")
print(metadata.columns.tolist())

# Normalize identifiers.

def norm_id(value):

    s = str(value).strip().lower()

    if s.endswith(".0") and s[:-2].isdigit():
        s = s[:-2]

    return s.replace("-", "")


df["sample_key"] = df["sample_id"].map(norm_id)

metadata["sample_key"] = (
    metadata["sample_id_normalized"]
    .map(norm_id)
)

merged = df.merge(
    metadata,
    on="sample_key",
    how="left",
    suffixes=("", "_metadata"),
)

if len(merged) != 179:
    raise RuntimeError(
        "Metadata join changed patient count."
    )

# Prefer metadata stage/diagnosis.

# ---------------------------------------------------------------
# Resolve diagnosis/stage after the metadata join.
#
# `diagnosis` exists in both inputs, so pandas creates:
#     diagnosis
#     diagnosis_metadata
#
# `stage` only exists in the metadata file, so pandas keeps:
#     stage
#
# Therefore we must NOT assume `stage_metadata` exists.
# ---------------------------------------------------------------

if "diagnosis_metadata" in merged.columns:
    merged["diagnosis_final"] = (
        merged["diagnosis_metadata"]
        .fillna(merged.get("diagnosis"))
    )
elif "diagnosis" in merged.columns:
    merged["diagnosis_final"] = merged["diagnosis"]
else:
    raise RuntimeError(
        "Could not find a diagnosis column after metadata join."
    )

if "stage_metadata" in merged.columns:
    merged["stage_final"] = merged["stage_metadata"]
elif "stage" in merged.columns:
    merged["stage_final"] = merged["stage"]
else:
    raise RuntimeError(
        "Could not find a stage column after metadata join."
    )

# Construct strata explicitly.

def get_stratum(row):

    diagnosis = str(
        row.get("diagnosis_final", "")
    ).strip().lower()

    stage_value = row.get("stage_final")
    stage = str(stage_value).strip().lower()

    # ---------------------------------------------------------------
    # E-1056 is the only historical cancer patient whose GEO metadata
    # lacks a stage annotation.
    #
    # Published Table 2 identifies E-1056 as Stage I.
    # Therefore it belongs to the Stage I/II stratum.
    # ---------------------------------------------------------------

    sample_values = []

    for field in (
        "sample_id",
        "sample_key",
        "sample_id_normalized",
        "title",
    ):
        value = row.get(field)

        if value is not None:
            sample_values.append(
                str(value).strip().lower()
            )

    is_e1056 = any(
        value.replace("-", "") == "e1056"
        for value in sample_values
    )

    if diagnosis == "control":
        return "Control"

    if diagnosis == "benign":
        return "Benign"

    if diagnosis == "borderline":
        return "Borderline"

    if diagnosis == "cancer":

        if is_e1056 and (
            not stage
            or stage in {"nan", "none"}
        ):
            return "Stage I/II"

        if "stage i/ii" in stage:
            return "Stage I/II"

        if "stage iii/iv" in stage:
            return "Stage III/IV"

        if stage in {"i", "ii", "1", "2"}:
            return "Stage I/II"

        if stage in {"iii", "iv", "3", "4"}:
            return "Stage III/IV"

        raise RuntimeError(
            "Cancer patient has an unrecognized or missing stage: "
            f"sample_values={sample_values!r}, "
            f"stage={stage_value!r}"
        )

    raise RuntimeError(
        "Could not classify patient: "
        f"sample_values={sample_values!r}, "
        f"diagnosis={diagnosis!r}, "
        f"stage={stage_value!r}"
    )



merged["stratum"] = merged.apply(
    get_stratum,
    axis=1,
)

print()
print("Available strata:")
print(
    merged["stratum"]
    .value_counts(dropna=False)
    .to_string()
)

if merged["stratum"].isna().any():

    bad = merged.loc[
        merged["stratum"].isna(),
        [
            "sample_id",
            "diagnosis_final",
            "stage_final",
        ],
    ]

    print()
    print("UNCLASSIFIED PATIENTS:")
    print(bad.to_string(index=False))

    raise RuntimeError(
        "Could not assign every patient to one of the five historical strata."
    )

TARGET_COUNTS = {
    "Control": 4,
    "Benign": 11,
    "Borderline": 5,
    "Stage I/II": 12,
    "Stage III/IV": 12,
}

for stratum, n in TARGET_COUNTS.items():

    available = int(
        (merged["stratum"] == stratum).sum()
    )

    print(
        f"{stratum:15s}: available={available:3d} "
        f"required={n:3d}"
    )

    if available < n:
        raise RuntimeError(
            f"Not enough patients in stratum {stratum}."
        )

# ---------------------------------------------------------------
# Search stratified random splits.
# ---------------------------------------------------------------

rng = np.random.default_rng(94533)

records = []

for trial in range(N_TRIALS):

    if trial % 1000 == 0:
        print(
            f"Progress: {trial:,}/{N_TRIALS:,} trials",
            flush=True,
        )

    selected = []

    for stratum, n in TARGET_COUNTS.items():

        candidates = merged.index[
            merged["stratum"] == stratum
        ].to_numpy()

        chosen = rng.choice(
            candidates,
            size=n,
            replace=False,
        )

        selected.extend(chosen.tolist())

    test = merged.loc[selected]

    auc = roc_auc_score(
        test["true_label"],
        test["cancer_activation"],
    )

    records.append(
        {
            "trial": trial,
            "auc": auc,
            "n_test": len(test),
            "controls": int(
                (test["stratum"] == "Control").sum()
            ),
            "benign": int(
                (test["stratum"] == "Benign").sum()
            ),
            "borderline": int(
                (test["stratum"] == "Borderline").sum()
            ),
            "stage_I_II": int(
                (test["stratum"] == "Stage I/II").sum()
            ),
            "stage_III_IV": int(
                (test["stratum"] == "Stage III/IV").sum()
            ),
            "sample_ids": "|".join(
                sorted(
                    test["sample_id"].astype(str)
                )
            ),
        }
    )


results = pd.DataFrame(records)

results = results.sort_values(
    by="auc",
    key=lambda s: (s - 0.90).abs(),
)

OUTPUT.parent.mkdir(
    parents=True,
    exist_ok=True,
)

results.to_csv(
    OUTPUT,
    index=False,
)

print()
print("=" * 90)
print("SEARCH COMPLETE")
print("=" * 90)

print()
print(f"Trials searched: {N_TRIALS}")

print()
print("Closest candidate splits to published AUC = 0.90:")

print(
    results[
        [
            "trial",
            "auc",
            "n_test",
            "controls",
            "benign",
            "borderline",
            "stage_I_II",
            "stage_III_IV",
        ]
    ]
    .head(20)
    .to_string(index=False)
)

print()
print("AUC distribution:")
print(
    results["auc"]
    .describe()
    .to_string()
)

print()
print("Saved:")
print(OUTPUT)

print()
print("=" * 90)
print("IMPORTANT")
print("=" * 90)

print(
    "A candidate with AUC near 0.90 is NOT proof that it was the "
    "original 44-patient test set."
)

print(
    "The purpose of this experiment is to determine whether the "
    "published AUC and demographic constraints substantially narrow "
    "the possible splits."
)
