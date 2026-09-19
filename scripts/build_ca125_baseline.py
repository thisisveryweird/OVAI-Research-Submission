#!/usr/bin/env python3

"""
CA-125 BASELINE

Evaluates CA-125 alone on the exact 179-patient historical cohort.

This is a NON-ML clinical baseline.

Positive prediction:
    higher CA-125 -> greater probability of cancer

Patients:
    179
    97 cancer
    82 non-cancer

CA-125 values are extracted from published eLife Table 2.

Missing CA-125 values are handled using median imputation.

Evaluation:
    5-fold patient-level OOF predictions

The purpose is to establish a clinically meaningful baseline
before adding miRNA features.
"""

from __future__ import annotations

import json
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold


ROOT = Path.home() / "ovarian-ai"

DOCX = (
    ROOT
    / "data"
    / "external"
    / "elife_28932"
    / "elife_poa_e28932_Supplementary_File_1.docx"
)

HISTORICAL_PREDICTIONS = (
    ROOT
    / "reports"
    / "tables"
    / "original_historical_clean_predictions.csv"
)

OUTPUT_DIR = (
    ROOT
    / "reports"
    / "tables"
)

PREDICTIONS_OUTPUT = (
    OUTPUT_DIR
    / "ca125_baseline_oof_predictions.csv"
)

RESULTS_OUTPUT = (
    OUTPUT_DIR
    / "ca125_baseline_results.json"
)


def norm_id(value):

    s = str(value).strip().lower()

    if re.fullmatch(r"\d+\.0+", s):
        s = s.split(".", 1)[0]

    return re.sub(
        r"[^a-z0-9]",
        "",
        s,
    )


def cell_text(cell):

    parts = []

    for node in cell.iter():

        if node.tag.endswith("}t"):
            parts.append(node.text or "")

    return "".join(parts).strip()


def extract_tables(path):

    with zipfile.ZipFile(path) as z:
        xml = z.read(
            "word/document.xml"
        )

    root = ET.fromstring(xml)

    tables = []

    for tbl in root.iter():

        if not tbl.tag.endswith("}tbl"):
            continue

        rows = []

        for tr in tbl:

            if not tr.tag.endswith("}tr"):
                continue

            cells = []

            for tc in tr:

                if tc.tag.endswith("}tc"):
                    cells.append(
                        cell_text(tc)
                    )

            if cells:
                rows.append(cells)

        if rows:
            tables.append(rows)

    return tables


# =====================================================================
# LOAD EXACT HISTORICAL COHORT
# =====================================================================

print("=" * 90)
print("CA-125 BASELINE")
print("=" * 90)

pred = pd.read_csv(
    HISTORICAL_PREDICTIONS
)

if len(pred) != 179:
    raise RuntimeError(
        f"Expected 179 patients, found {len(pred)}."
    )

pred["sample_key"] = (
    pred["sample_id"].map(norm_id)
)

label_map = (
    pred
    .set_index("sample_key")["true_label"]
    .astype(int)
)

print()
print("Historical cohort:", len(pred))
print(
    "Cancer:",
    int((pred.true_label == 1).sum())
)
print(
    "Non-cancer:",
    int((pred.true_label == 0).sum())
)


# =====================================================================
# READ TABLE 2
# =====================================================================

print()
print("Reading published Table 2...")

tables = extract_tables(DOCX)

if len(tables) < 2:
    raise RuntimeError(
        "Could not find two supplementary tables."
    )

table2 = tables[1]

header = table2[0]

try:
    id_col = header.index("ID")
except ValueError:
    raise RuntimeError(
        "Could not locate ID column."
    )

try:
    ca_col = header.index("CA-125")
except ValueError:
    raise RuntimeError(
        "Could not locate CA-125 column."
    )


# =====================================================================
# PARSE PATIENT ROWS
# =====================================================================

rows = []

for row in table2[1:]:

    if len(row) <= max(
        id_col,
        ca_col,
    ):
        continue

    sid = row[id_col].strip()

    if not sid:
        continue

    if sid.lower() == "id":
        continue

    ca125 = row[ca_col].strip()

    rows.append(
        {
            "table2_id": sid,
            "ca125_raw": ca125,
        }
    )


ca = pd.DataFrame(rows)

print()
print(
    "Table 2 patient rows:",
    len(ca),
)


# =====================================================================
# IMPORTANT:
# TABLE 2 IDs ARE ACTUAL SAMPLE IDS
# =====================================================================

# ---------------------------------------------------------------
# Table 2 uses TWO kinds of identifiers:
#
#   numeric 1..120  = historical patient positions
#   e1001...        = actual e-sample IDs
#
# The historical PMML reproduction established that original
# numeric position 56 was excluded because it corresponds to 81195.
#
# Therefore:
#
#   Table 2 ID 1..55  -> prediction rows 1..55
#   Table 2 ID 57..120 -> prediction rows 56..119
#
# e-samples retain their explicit IDs.
# ---------------------------------------------------------------

prediction_by_key = (
    pred
    .set_index("sample_key")
)

def table2_to_actual_sample_id(value):

    key = norm_id(value)

    # Explicit e-sample.
    if key.startswith("e"):
        return key

    # Numeric Table-2 position.
    try:
        position = int(float(value))
    except (TypeError, ValueError):
        raise RuntimeError(
            f"Could not interpret Table 2 ID: {value!r}"
        )

    if position == 56:
        raise RuntimeError(
            "Table 2 unexpectedly contains excluded position 56."
        )

    if position < 1 or position > 120:
        raise RuntimeError(
            f"Unexpected numeric Table 2 ID: {position}"
        )

    # Convert original historical position to retained prediction
    # row after excluding original position 56.
    retained_index = (
        position - 1
        if position <= 55
        else position - 2
    )

    return pred.iloc[
        retained_index
    ]["sample_key"]


ca["actual_sample_key"] = (
    ca["table2_id"].map(
        table2_to_actual_sample_id
    )
)

ca["true_label"] = (
    ca["actual_sample_key"].map(
        label_map
    )
)


if ca["true_label"].isna().any():

    missing = ca.loc[
        ca["true_label"].isna(),
        "table2_id",
    ].tolist()

    raise RuntimeError(
        "Missing labels for Table 2 IDs:\n"
        + "\n".join(
            map(str, missing)
        )
    )


# =====================================================================
# PARSE CA-125
# =====================================================================

def parse_ca125(value):

    s = str(value).strip()

    if (
        not s
        or s in {"‡", "*", "nan"}
    ):
        return np.nan

    try:
        return float(s)
    except ValueError:
        return np.nan


ca["ca125"] = (
    ca["ca125_raw"].map(
        parse_ca125
    )
)


# =====================================================================
# ALIGN EXACTLY TO HISTORICAL COHORT
# =====================================================================

historical = (
    pred[
        [
            "sample_id",
            "sample_key",
            "true_label",
        ]
    ]
    .merge(
        ca[
            [
                "actual_sample_key",
                "ca125",
            ]
        ],
        left_on="sample_key",
        right_on="actual_sample_key",
        how="left",
        validate="one_to_one",
    )
)


if len(historical) != 179:
    raise RuntimeError(
        f"Expected 179 aligned patients, got {len(historical)}."
    )


# =====================================================================
# VALIDATE
# =====================================================================

y = historical[
    "true_label"
].to_numpy(dtype=int)

x_raw = historical[
    ["ca125"]
].to_numpy(dtype=float)

print()
print("=" * 90)
print("CA-125 DATA")
print("=" * 90)

print()
print(
    "Patients:",
    len(historical),
)

print(
    "Cancer:",
    int((y == 1).sum())
)

print(
    "Non-cancer:",
    int((y == 0).sum())
)

print(
    "Observed CA-125:",
    int(np.isfinite(x_raw).sum())
)

print(
    "Missing CA-125:",
    int(np.isnan(x_raw).sum())
)


# =====================================================================
# OOF MEDIAN IMPUTATION
# =====================================================================

cv = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=2026,
)

oof = np.full(
    len(y),
    np.nan,
)


for fold, (
    train_idx,
    test_idx,
) in enumerate(
    cv.split(
        x_raw,
        y,
    ),
    start=1,
):

    imputer = SimpleImputer(
        strategy="median"
    )

    x_train = imputer.fit_transform(
        x_raw[train_idx]
    )

    x_test = imputer.transform(
        x_raw[test_idx]
    )

    # CA-125 itself is the score.
    # Higher CA-125 => greater cancer likelihood.
    oof[test_idx] = (
        x_test[:, 0]
    )


# =====================================================================
# METRICS
# =====================================================================

auc = roc_auc_score(
    y,
    oof,
)

pr_auc = average_precision_score(
    y,
    oof,
)

threshold = 35.0

prediction = (
    oof > threshold
).astype(int)

accuracy = accuracy_score(
    y,
    prediction,
)

tn, fp, fn, tp = confusion_matrix(
    y,
    prediction,
).ravel()

sensitivity = (
    tp / (tp + fn)
)

specificity = (
    tn / (tn + fp)
)


# =====================================================================
# SAVE
# =====================================================================

historical[
    "oof_ca125"
] = oof

historical[
    "prediction_ca125_gt35"
] = prediction

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

historical.to_csv(
    PREDICTIONS_OUTPUT,
    index=False,
)


results = {
    "model": "CA-125 only",
    "cohort": {
        "patients": 179,
        "cancer": int((y == 1).sum()),
        "non_cancer": int((y == 0).sum()),
    },
    "missing_ca125": int(
        np.isnan(x_raw).sum()
    ),
    "evaluation": {
        "method": "5-fold patient-level OOF",
        "random_state": 2026,
    },
    "threshold": 35.0,
    "metrics": {
        "roc_auc": float(auc),
        "pr_auc": float(pr_auc),
        "accuracy": float(accuracy),
        "sensitivity": float(sensitivity),
        "specificity": float(specificity),
    },
    "confusion_matrix": {
        "TN": int(tn),
        "FP": int(fp),
        "FN": int(fn),
        "TP": int(tp),
    },
}


with RESULTS_OUTPUT.open(
    "w"
) as f:

    json.dump(
        results,
        f,
        indent=2,
    )


# =====================================================================
# REPORT
# =====================================================================

print()
print("=" * 90)
print("CA-125 BASELINE RESULTS")
print("=" * 90)

print()
print(f"ROC AUC:      {auc:.6f}")
print(f"PR AUC:       {pr_auc:.6f}")
print(f"Accuracy:     {accuracy:.6f}")
print(f"Sensitivity:  {sensitivity:.6f}")
print(f"Specificity:  {specificity:.6f}")

print()
print("Confusion matrix:")
print(f"  TP={tp} FN={fn}")
print(f"  FP={fp} TN={tn}")

print()
print("Saved:")
print(PREDICTIONS_OUTPUT)
print(RESULTS_OUTPUT)
