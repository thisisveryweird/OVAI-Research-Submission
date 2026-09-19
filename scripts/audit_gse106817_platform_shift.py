#!/usr/bin/env python3
"""Audit GSE106817 distribution shift without fitting on external labels."""

from __future__ import annotations

import csv
import gzip
import io
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
from sklearn.preprocessing import QuantileTransformer, StandardScaler


PROJECT = Path(__file__).resolve().parents[1]
TRAIN = PROJECT / "data" / "processed" / "gse94533_14mirna_clinical.csv.gz"
RAW = PROJECT / "data" / "external" / "gse106817_series_matrix.txt.gz"
PLATFORM = PROJECT / "data" / "external" / "GPL21263.txt"
RESULTS = PROJECT / "reports" / "tables" / "gse106817_external_results.json"
OUT_JSON = PROJECT / "reports" / "gse106817_platform_shift_audit.json"
OUT_REPORT = PROJECT / "reports" / "gse106817_platform_shift_audit.md"

MIRNAS = [
    "hsa-miR-23b-3p", "hsa-miR-29a-3p", "hsa-miR-32-5p", "hsa-miR-92a-3p",
    "hsa-miR-150-5p", "hsa-miR-200a-3p", "hsa-miR-200c-3p", "hsa-miR-203a",
    "hsa-miR-320c", "hsa-miR-320d", "hsa-miR-335-5p", "hsa-miR-450b-5p",
    "hsa-miR-1246", "hsa-miR-1307-5p",
]


def load_external() -> tuple[pd.DataFrame, dict[str, int]]:
    with gzip.open(RAW, "rt", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_begin"))
    metadata = {}
    for line in lines[:start]:
        if line.startswith("!Sample_"):
            parts = next(csv.reader([line.rstrip("\n")], delimiter="\t", quotechar='"'))
            metadata[parts[0]] = parts[1:]
    table = pd.read_csv(
        io.StringIO("".join(lines[start + 1:])),
        sep="\t",
        comment="!",
        quotechar='"',
    ).rename(columns=lambda value: value)
    table = table.rename(columns={table.columns[0]: "feature_id"})
    table = table[table.feature_id.notna() & (table.feature_id != "!series_matrix_table_end")]
    names = {}
    inside = False
    for line in PLATFORM.read_text(errors="replace").splitlines():
        if line.startswith("!platform_table_begin"):
            inside = True
            continue
        if inside and line.startswith("!platform_table_end"):
            break
        if inside:
            parts = re.split(r"\s+", line.strip())
            if len(parts) >= 3:
                names[parts[0]] = parts[2]
    feature_by_name = {name: feature for feature, name in names.items()}
    if "hsa-miR-203a-3p" in feature_by_name:
        feature_by_name["hsa-miR-203a"] = feature_by_name["hsa-miR-203a-3p"]
    missing = [name for name in MIRNAS if name not in feature_by_name]
    if missing:
        raise RuntimeError(f"Missing platform features: {missing}")
    sample_ids = metadata["!Sample_geo_accession"]
    titles = metadata["!Sample_title"]
    labels = []
    keep = []
    groups = {"cancer": 0, "benign": 0, "borderline": 0, "control": 0, "other": 0}
    for sample_id, title in zip(sample_ids, titles):
        title_lower = title.lower()
        if re.match(r"^Ovarian Cancer\b", title, re.I):
            keep.append(sample_id)
            labels.append(1)
            groups["cancer"] += 1
        elif re.match(r"^Benign Ovarian Disease\b", title, re.I):
            keep.append(sample_id)
            labels.append(0)
            groups["benign"] += 1
        elif re.match(r"^Borderline Ovarian Tumor\b", title, re.I):
            groups["borderline"] += 1
        elif re.match(r"^Non-Cancer\b", title, re.I):
            keep.append(sample_id)
            labels.append(0)
            groups["control"] += 1
        else:
            groups["other"] += 1
    values = table.set_index("feature_id").loc[[feature_by_name[name] for name in MIRNAS], keep].T
    values.columns = MIRNAS
    return values.apply(pd.to_numeric, errors="coerce"), groups


def main() -> None:
    train = pd.read_csv(TRAIN)
    external, groups = load_external()
    train_values = np.log1p(train[MIRNAS].clip(lower=0))
    external_values = np.log1p(external.clip(lower=0))
    scaler = StandardScaler().fit(train_values)
    quantile = QuantileTransformer(
        output_distribution="normal",
        n_quantiles=min(100, len(train_values)),
        random_state=42,
    ).fit(train_values)
    external_scaled = scaler.transform(external_values)
    external_quantiled = quantile.transform(external_values)
    rows = []
    for feature in MIRNAS:
        train_feature = train_values[feature].dropna().to_numpy()
        external_feature = external_values[feature].dropna().to_numpy()
        rows.append({
            "feature": feature,
            "train_median_log1p": float(np.median(train_feature)),
            "external_median_log1p": float(np.median(external_feature)),
            "median_shift_log1p": float(np.median(external_feature) - np.median(train_feature)),
            "train_q10_log1p": float(np.quantile(train_feature, 0.10)),
            "external_q10_log1p": float(np.quantile(external_feature, 0.10)),
            "train_q90_log1p": float(np.quantile(train_feature, 0.90)),
            "external_q90_log1p": float(np.quantile(external_feature, 0.90)),
            "external_median_source_z": float(np.median(external_scaled[:, MIRNAS.index(feature)])),
            "external_median_source_quantile": float(np.median(external_quantiled[:, MIRNAS.index(feature)])),
            "ks_statistic": float(ks_2samp(train_feature, external_feature).statistic),
            "ks_pvalue": float(ks_2samp(train_feature, external_feature).pvalue),
        })
    distribution = pd.DataFrame(rows)
    existing = json.loads(RESULTS.read_text()) if RESULTS.exists() else {}
    output = {
        "training_cohort": "GSE94533",
        "external_cohort": "GSE106817",
        "external_groups": groups,
        "training_patients": int(len(train)),
        "external_patients_evaluated": int(len(external)),
        "normalization_policy": "Only source-fitted transforms are assessed. No external labels are used to fit a transform, select features, tune hyperparameters, or choose a threshold.",
        "feature_shift": rows,
        "existing_frozen_results": existing,
        "interpretation": "Source-fitted scaling and quantile transformation change numeric scale but do not establish assay equivalence, patient independence, or a clinically valid operating threshold.",
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(output, indent=2) + "\n")
    largest = distribution.sort_values("ks_statistic", ascending=False).head(5)
    lines = [
        "# GSE106817 Platform-Shift Audit",
        "",
        "This audit uses GSE94533 as the training cohort and keeps GSE106817 untouched. It does not fit a transform on external labels and does not add GSE106817 to training.",
        "",
        f"GSE106817 samples evaluated: {len(external)}. Full title-derived groups: {groups}.",
        "",
        "## What Normalization Can and Cannot Do",
        "",
        "The existing external evaluation already fits standard scaling and a source-fitted quantile transform on GSE94533 only. This audit quantifies feature distribution differences and records those frozen results. These transforms can alter numeric scale; they cannot prove that GPL21263 measurements are assay-equivalent to the GSE94533 processed serum miRNA measurements, recover missing batch covariates, or make an external cohort part of the training pool.",
        "",
        "## Largest Distribution Differences",
        "",
    ]
    for row in largest.itertuples(index=False):
        lines.append(f"- `{row.feature}`: KS={row.ks_statistic:.3f}, median log1p shift={row.median_shift_log1p:.3f}, source-fitted external median z={row.external_median_source_z:.3f}.")
    lines.extend([
        "",
        "## Frozen External Result",
        "",
        "The existing standard-scaling result has ROC-AUC 0.8232 but predicts every evaluated external sample as non-cancer at the source-derived threshold: sensitivity 0%, specificity 100%, balanced accuracy 0.5. The source-fitted quantile result has ROC-AUC 0.6893 and the same threshold failure.",
        "",
        "## Decision",
        "",
        "GSE106817 remains an untouched external stress test and platform-shift study. It is not added to training, and normalization alone is not accepted as a fix. A valid next experiment requires prespecified assay harmonization, a positive-control operating point defined without target labels, overlap auditing, and replication on another independent serum cohort.",
        "",
        "Machine-readable output: `reports/gse106817_platform_shift_audit.json`.",
    ])
    OUT_REPORT.write_text("\n".join(lines) + "\n")
    print(f"Saved {OUT_REPORT}")
    print(f"Saved {OUT_JSON}")


if __name__ == "__main__":
    main()
