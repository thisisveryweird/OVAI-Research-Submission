#!/usr/bin/env python3
"""Run bidirectional and platform-shift experiments for GSE94533/GSE106817."""

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
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import QuantileTransformer, StandardScaler


PROJECT = Path(__file__).resolve().parents[1]
GSE94533 = PROJECT / "data" / "processed" / "gse94533_14mirna_clinical.csv.gz"
GSE106817 = PROJECT / "data" / "external" / "gse106817_series_matrix.txt.gz"
GPL21263 = PROJECT / "data" / "external" / "GPL21263.txt"
OUT_RESULTS = PROJECT / "reports" / "gse106817_phase5_results.json"
OUT_PREDICTIONS = PROJECT / "reports" / "gse106817_phase5_predictions.csv"
OUT_DISTRIBUTIONS = PROJECT / "reports" / "gse106817_phase5_distribution.csv"
OUT_REPORT = PROJECT / "reports" / "gse106817_phase5.md"

SEED = 20260919
C_GRID = [0.03, 0.1, 0.3, 1.0]
THRESHOLDS = np.linspace(0.10, 0.90, 81)
MIRNAS = [
    "hsa-miR-23b-3p", "hsa-miR-29a-3p", "hsa-miR-32-5p", "hsa-miR-92a-3p",
    "hsa-miR-150-5p", "hsa-miR-200a-3p", "hsa-miR-200c-3p", "hsa-miR-203a",
    "hsa-miR-320c", "hsa-miR-320d", "hsa-miR-335-5p", "hsa-miR-450b-5p",
    "hsa-miR-1246", "hsa-miR-1307-5p",
]
MODES = [
    "raw_values",
    "log1p_values",
    "within_dataset_standard",
    "training_derived_standard",
    "within_dataset_quantile",
    "training_derived_quantile",
]


def load_gse94533() -> tuple[pd.DataFrame, np.ndarray]:
    frame = pd.read_csv(GSE94533)
    return frame[MIRNAS].apply(pd.to_numeric, errors="coerce"), frame["binary_label"].astype(int).to_numpy()


def load_gse106817() -> tuple[pd.DataFrame, np.ndarray, dict[str, int]]:
    with gzip.open(GSE106817, "rt", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_begin"))
    metadata = {}
    for line in lines[:start]:
        if line.startswith("!Sample_"):
            parts = next(csv.reader([line.rstrip("\n")], delimiter="\t", quotechar='"'))
            metadata[parts[0]] = parts[1:]
    table = pd.read_csv(io.StringIO("".join(lines[start + 1:])), sep="\t", comment="!", quotechar='"')
    table = table.rename(columns={table.columns[0]: "feature_id"})
    table = table[table.feature_id.notna() & (table.feature_id != "!series_matrix_table_end")]
    names = {}
    inside = False
    for line in GPL21263.read_text(errors="replace").splitlines():
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
        raise RuntimeError(f"Missing GSE106817 features: {missing}")
    sample_ids = metadata["!Sample_geo_accession"]
    titles = metadata["!Sample_title"]
    keep = []
    labels = []
    groups = {"cancer": 0, "benign": 0, "borderline": 0, "control": 0, "other": 0}
    for sample_id, title in zip(sample_ids, titles):
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
    matrix = table.set_index("feature_id").loc[[feature_by_name[name] for name in MIRNAS], keep].T
    matrix.columns = MIRNAS
    return matrix.apply(pd.to_numeric, errors="coerce"), np.asarray(labels, dtype=int), groups


def transform_pair(train: pd.DataFrame, test: pd.DataFrame, mode: str) -> tuple[np.ndarray, np.ndarray]:
    train_values = train.to_numpy(dtype=float)
    test_values = test.to_numpy(dtype=float)
    if mode.endswith("_values"):
        if mode == "log1p_values":
            train_values = np.log1p(np.clip(train_values, 0, None))
            test_values = np.log1p(np.clip(test_values, 0, None))
        return train_values, test_values
    if mode.startswith("within_dataset"):
        train_values = np.log1p(np.clip(train_values, 0, None))
        test_values = np.log1p(np.clip(test_values, 0, None))
        transformer = StandardScaler() if mode.endswith("standard") else QuantileTransformer(
            output_distribution="normal", n_quantiles=min(100, len(train_values)), random_state=SEED
        )
        return transformer.fit_transform(train_values), transformer.fit_transform(test_values)
    train_values = np.log1p(np.clip(train_values, 0, None))
    test_values = np.log1p(np.clip(test_values, 0, None))
    transformer = StandardScaler() if mode.endswith("standard") else QuantileTransformer(
        output_distribution="normal", n_quantiles=min(100, len(train_values)), random_state=SEED
    )
    return transformer.fit_transform(train_values), transformer.transform(test_values)


def choose_c_and_threshold(X: pd.DataFrame, y: np.ndarray, mode: str) -> tuple[float, float, np.ndarray]:
    outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    cv_scores = {c: [] for c in C_GRID}
    oof_by_c = {c: np.full(len(y), np.nan) for c in C_GRID}
    for train_idx, validation_idx in outer.split(X, y):
        for c in C_GRID:
            train_values, validation_values = transform_pair(X.iloc[train_idx], X.iloc[validation_idx], mode)
            model = LogisticRegression(C=c, class_weight="balanced", solver="liblinear", max_iter=5000)
            model.fit(train_values, y[train_idx])
            probability = model.predict_proba(validation_values)[:, 1]
            oof_by_c[c][validation_idx] = probability
            cv_scores[c].append(roc_auc_score(y[validation_idx], probability))
    best_c = max(C_GRID, key=lambda c: (np.mean(cv_scores[c]), -abs(np.log10(c) - np.log10(0.3))))
    source_probability = oof_by_c[best_c]
    threshold = max(
        THRESHOLDS,
        key=lambda candidate: balanced_accuracy_score(y, source_probability >= candidate),
    )
    return float(best_c), float(threshold), source_probability


def metrics(y: np.ndarray, probability: np.ndarray, threshold: float) -> dict[str, object]:
    prediction = (probability >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, prediction, labels=[0, 1]).ravel()
    return {
        "n": int(len(y)),
        "positive": int(y.sum()),
        "negative": int((y == 0).sum()),
        "roc_auc": float(roc_auc_score(y, probability)),
        "pr_auc": float(average_precision_score(y, probability)),
        "accuracy": float(accuracy_score(y, prediction)),
        "balanced_accuracy": float(balanced_accuracy_score(y, prediction)),
        "sensitivity": float(tp / (tp + fn)) if tp + fn else None,
        "specificity": float(tn / (tn + fp)) if tn + fp else None,
        "brier_score": float(brier_score_loss(y, probability)),
        "threshold": float(threshold),
        "confusion_matrix": {"tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp)},
    }


def run_direction(train_name: str, train_X: pd.DataFrame, train_y: np.ndarray, test_name: str, test_X: pd.DataFrame, test_y: np.ndarray) -> tuple[list[dict[str, object]], list[pd.DataFrame]]:
    results = []
    prediction_rows = []
    for mode in MODES:
        best_c, threshold, source_probability = choose_c_and_threshold(train_X, train_y, mode)
        train_values, test_values = transform_pair(train_X, test_X, mode)
        model = LogisticRegression(C=best_c, class_weight="balanced", solver="liblinear", max_iter=5000)
        model.fit(train_values, train_y)
        probability = model.predict_proba(test_values)[:, 1]
        result = {
            "experiment": "A" if train_name == "GSE94533" else "B",
            "train_dataset": train_name,
            "test_dataset": test_name,
            "normalization": mode,
            "best_C_source_only": best_c,
            "source_threshold": threshold,
            **metrics(test_y, probability, threshold),
        }
        results.append(result)
        prediction_rows.append(pd.DataFrame({
            "experiment": result["experiment"],
            "train_dataset": train_name,
            "test_dataset": test_name,
            "normalization": mode,
            "sample_index": np.arange(len(test_y)),
            "true_label": test_y,
            "probability": probability,
            "threshold": threshold,
            "prediction": (probability >= threshold).astype(int),
        }))
    return results, prediction_rows


def main() -> None:
    gse94533_X, gse94533_y = load_gse94533()
    gse106817_X, gse106817_y, external_groups = load_gse106817()
    results_a, predictions_a = run_direction("GSE94533", gse94533_X, gse94533_y, "GSE106817", gse106817_X, gse106817_y)
    results_b, predictions_b = run_direction("GSE106817", gse106817_X, gse106817_y, "GSE94533", gse94533_X, gse94533_y)
    distributions = []
    for feature in MIRNAS:
        source = np.log1p(np.clip(gse94533_X[feature].to_numpy(), 0, None))
        target = np.log1p(np.clip(gse106817_X[feature].to_numpy(), 0, None))
        ks = ks_2samp(source, target)
        distributions.append({
            "feature": feature,
            "gse94533_median_log1p": float(np.median(source)),
            "gse106817_median_log1p": float(np.median(target)),
            "median_shift_log1p_gse106817_minus_gse94533": float(np.median(target) - np.median(source)),
            "ks_statistic": float(ks.statistic),
            "ks_pvalue": float(ks.pvalue),
        })
    result_payload = {
        "shared_features": MIRNAS,
        "experiments": results_a + results_b,
        "experiment_C": {
            "status": "not_distinctly_executable",
            "compatible_datasets_available": ["GSE94533", "GSE106817"],
            "reason": "Only two compatible serum-miRNA datasets are available. Pooling both and holding one out is identical to Experiments A and B; there is no third dataset to rotate.",
        },
        "external_groups": external_groups,
        "preprocessing_rule": "All model selection, transforms, source thresholds, and hyperparameters use training data only. External labels are never used for preprocessing or selection.",
    }
    OUT_RESULTS.parent.mkdir(parents=True, exist_ok=True)
    OUT_RESULTS.write_text(json.dumps(result_payload, indent=2) + "\n")
    pd.concat(predictions_a + predictions_b, ignore_index=True).to_csv(OUT_PREDICTIONS, index=False)
    pd.DataFrame(distributions).to_csv(OUT_DISTRIBUTIONS, index=False)
    lines = [
        "# Phase 5: GSE106817 Bidirectional Validation",
        "",
        "This evaluation uses only shared historical miRNAs and never uses external test labels for preprocessing, model selection, normalization, or threshold selection.",
        "",
        "## Experiment A: GSE94533 -> GSE106817",
        "",
    ]
    for row in results_a:
        lines.append(f"- `{row['normalization']}`: ROC-AUC {row['roc_auc']:.4f}, balanced accuracy {row['balanced_accuracy']:.4f}, sensitivity {row['sensitivity']:.4f}, specificity {row['specificity']:.4f}, threshold {row['source_threshold']:.2f}.")
    lines.extend(["", "## Experiment B: GSE106817 -> GSE94533", ""])
    for row in results_b:
        lines.append(f"- `{row['normalization']}`: ROC-AUC {row['roc_auc']:.4f}, balanced accuracy {row['balanced_accuracy']:.4f}, sensitivity {row['sensitivity']:.4f}, specificity {row['specificity']:.4f}, threshold {row['source_threshold']:.2f}.")
    lines.extend([
        "",
        "## Experiment C: Pooled Leave-One-Dataset-Out",
        "",
        "Only two compatible serum-miRNA datasets are available. A pooled train/held-out rotation with these two datasets is exactly Experiments A and B, so no additional pooled result is reported as if it were an independent third rotation.",
        "",
        "## Experiment D: Platform Shift",
        "",
        "The distribution table compares log1p shared-miRNA values using medians and two-sample KS statistics. The raw, within-dataset, training-derived, and quantile modes are evaluated above. Within-dataset normalization uses unlabeled target feature distributions and is therefore explicitly labeled transductive; it does not use target labels.",
        "",
        "The external cohort remains untouched. A high ranking score with a failed source-derived threshold is reported as threshold/platform failure, not hidden or converted into diagnostic success.",
        "",
        "Machine-readable outputs: `reports/gse106817_phase5_results.json`, `reports/gse106817_phase5_predictions.csv`, and `reports/gse106817_phase5_distribution.csv`.",
    ])
    OUT_REPORT.write_text("\n".join(lines) + "\n")
    print(f"Saved {OUT_REPORT}")
    print(f"Saved {OUT_RESULTS}")
    print(f"Saved {OUT_PREDICTIONS}")
    print(f"Saved {OUT_DISTRIBUTIONS}")


if __name__ == "__main__":
    main()
