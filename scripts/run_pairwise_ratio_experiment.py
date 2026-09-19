#!/usr/bin/env python3
"""Evaluate a closed-loop pairwise miRNA representation across GEO cohorts.

Each sample is represented only by the ten ordered comparisons among the
frozen five-miRNA signature. No external labels or cohort-level statistics are
used to construct the features.
"""

from __future__ import annotations

import argparse
import logging
from itertools import combinations
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

try:
    from run_rank_transform_experiment import (
        EXTERNAL_ACCESSIONS,
        Cohort,
        load_external,
        load_local,
    )
except ModuleNotFoundError:  # Supports package-style imports in tests.
    from scripts.run_rank_transform_experiment import (
        EXTERNAL_ACCESSIONS,
        Cohort,
        load_external,
        load_local,
    )


PROJECT = Path(__file__).resolve().parents[1]
REPORT = PROJECT / "reports" / "pairwise_ratio_experiment.md"
SEED = 20260919
SIGNATURE = (
    "hsa-miR-200b-3p",
    "hsa-miR-296-3p",
    "hsa-miR-6826-5p",
    "hsa-miR-32-3p",
    "hsa-miR-1255a",
)
PAIRS = tuple(combinations(SIGNATURE, 2))
LOGGER = logging.getLogger("ovai.pairwise_ratio")


def pairwise_features(values: pd.DataFrame) -> pd.DataFrame:
    """Return ten binary A>B features using only each row's five values."""

    missing = [feature for feature in SIGNATURE if feature not in values.columns]
    if missing:
        raise ValueError(f"Missing pairwise signature features: {missing}")
    numeric = values.loc[:, SIGNATURE].apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any():
        missing_rows = int(numeric.isna().any(axis=1).sum())
        raise ValueError(f"Pairwise input contains missing/non-numeric target values in {missing_rows} rows")
    output = {
        f"{left}>{right}": (numeric[left].to_numpy() > numeric[right].to_numpy()).astype(np.int8)
        for left, right in PAIRS
    }
    return pd.DataFrame(output, index=values.index)


def classifier() -> Pipeline:
    return Pipeline([
        ("scale", StandardScaler()),
        ("model", LogisticRegression(C=0.3, class_weight="balanced", solver="liblinear", max_iter=5000, random_state=SEED)),
    ])


def calculate_metrics(labels: np.ndarray, probabilities: np.ndarray) -> dict[str, float]:
    predictions = (probabilities >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(labels, predictions, labels=[0, 1]).ravel()
    return {
        "roc_auc": float(roc_auc_score(labels, probabilities)),
        "balanced_accuracy": float(balanced_accuracy_score(labels, predictions)),
        "sensitivity": float(tp / (tp + fn)) if tp + fn else np.nan,
        "specificity": float(tn / (tn + fp)) if tn + fp else np.nan,
    }


def cross_validate(values: pd.DataFrame, labels: np.ndarray) -> dict[str, float]:
    features = pairwise_features(values)
    probabilities = np.full(len(labels), np.nan)
    splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    for train_index, validation_index in splitter.split(features, labels):
        fitted = classifier().fit(features.iloc[train_index], labels[train_index])
        probabilities[validation_index] = fitted.predict_proba(features.iloc[validation_index])[:, 1]
    return calculate_metrics(labels, probabilities)


def evaluate_external(model: Pipeline, cohort: Cohort) -> dict[str, float]:
    return calculate_metrics(
        cohort.labels,
        model.predict_proba(pairwise_features(cohort.values))[:, 1],
    )


def format_metric(value: object) -> str:
    return f"{value:.4f}" if isinstance(value, (float, int)) and np.isfinite(value) else "-"


def write_report(rows: list[dict[str, object]]) -> None:
    lines = [
        "# Pairwise Ratio Experiment",
        "",
        f"Frozen signature: `{', '.join(SIGNATURE)}`.",
        f"Pairwise feature count: {len(PAIRS)}. A feature is 1 when the left miRNA expression exceeds the right miRNA expression, otherwise 0.",
        "External cohorts are evaluated with a model fitted only on the local GSE94533 cohort.",
        "",
        "| Cohort | Evaluation | ROC-AUC | Balanced Accuracy | Sensitivity | Specificity | Detail |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        detail = str(row.get("detail", "")).replace("|", "\\|")
        lines.append(
            f"| {row['cohort']} | {row['evaluation']} | {format_metric(row.get('roc_auc'))} | "
            f"{format_metric(row.get('balanced_accuracy'))} | {format_metric(row.get('sensitivity'))} | "
            f"{format_metric(row.get('specificity'))} | {detail} |"
        )
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines) + "\n")


def run() -> None:
    local_values, local_labels = load_local()
    local_pairwise = pairwise_features(local_values)
    if local_pairwise.shape[1] != 10:
        raise RuntimeError(f"Expected 10 pairwise features, found {local_pairwise.shape[1]}")
    local_cv = cross_validate(local_values, local_labels)
    frozen_model = classifier().fit(local_pairwise, local_labels)
    rows: list[dict[str, object]] = [{
        "cohort": "GSE94533",
        "evaluation": "5-fold stratified CV",
        **local_cv,
    }]
    for accession in EXTERNAL_ACCESSIONS[:3]:
        try:
            cohort = load_external(accession)
            external_features = pairwise_features(cohort.values)
            LOGGER.info("%s: %d samples, %d pairwise features", accession, len(external_features), external_features.shape[1])
            rows.append({
                "cohort": accession,
                "evaluation": "frozen local model",
                **evaluate_external(frozen_model, cohort),
            })
        except (FileNotFoundError, ValueError) as exc:
            LOGGER.warning("Skipping %s: %s", accession, exc)
            rows.append({
                "cohort": accession,
                "evaluation": "not evaluable",
                "detail": str(exc),
            })
    table = pd.DataFrame(rows)
    print(table.to_string(index=False, float_format=lambda value: f"{value:.4f}"))
    write_report(rows)
    LOGGER.info("Wrote %s", REPORT)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())