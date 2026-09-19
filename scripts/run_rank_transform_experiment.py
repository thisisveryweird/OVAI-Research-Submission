#!/usr/bin/env python3
"""Evaluate a within-patient miRNA rank transform across cached GEO cohorts.

Ranks are computed independently within each sample against that sample's
complete available ``hsa-`` expression background. No transform parameters,
labels, thresholds, or external observations are learned from the target
cohort.
"""

from __future__ import annotations

import argparse
import logging
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import GEOparse
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

try:
    from bulk_test_external_datasets import label_from_metadata
except ModuleNotFoundError:  # Supports package-style imports in tests.
    from scripts.bulk_test_external_datasets import label_from_metadata

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=pd.errors.DtypeWarning)

PROJECT = Path(__file__).resolve().parents[1]
LOCAL_MATRIX = PROJECT / "data" / "processed" / "gse94533_mirna_tpm.csv.gz"
LOCAL_METADATA = PROJECT / "data" / "processed" / "gse94533_metadata.csv"
GEO_CACHE = PROJECT / "data" / "raw" / "geo"
REPORT = PROJECT / "reports" / "rank_transform_experiment.md"
SEED = 20260919
SIGNATURE = (
    "hsa-miR-200b-3p",
    "hsa-miR-296-3p",
    "hsa-miR-6826-5p",
    "hsa-miR-32-3p",
    "hsa-miR-1255a",
)
EXTERNAL_ACCESSIONS = ("GSE106817", "GSE113486", "GSE113740", "GSE138092")
LOGGER = logging.getLogger("ovai.rank_transform")


@dataclass(frozen=True)
class Cohort:
    name: str
    values: pd.DataFrame
    labels: np.ndarray
    background: pd.DataFrame | None = None
    probe_map: Mapping[str, tuple[str, ...]] | None = None


class WithinPatientBackgroundRankTransformer:
    """Convert target features to per-sample percentile ranks.

    ``fit`` is deliberately a no-op. Ranking depends only on each row's
    expression values, so calling it separately inside every CV fold makes the
    no-leakage contract explicit without learning cross-sample statistics.
    """

    def __init__(self, target_features: tuple[str, ...]) -> None:
        self.target_features = target_features

    def fit(self, background: pd.DataFrame, y: np.ndarray | None = None) -> "WithinPatientBackgroundRankTransformer":
        missing = [feature for feature in self.target_features if feature not in background.columns]
        if missing:
            raise ValueError(f"Missing target features from background: {missing}")
        return self

    def transform(self, background: pd.DataFrame) -> pd.DataFrame:
        missing = [feature for feature in self.target_features if feature not in background.columns]
        if missing:
            raise ValueError(f"Missing target features from background: {missing}")
        numeric = background.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        output = np.full((len(background), len(self.target_features)), np.nan, dtype=float)
        target_indices = [background.columns.get_loc(feature) for feature in self.target_features]
        for row_index, row in enumerate(numeric):
            finite = np.isfinite(row)
            if not finite.any():
                continue
            ranks = rankdata(row[finite], method="average") / finite.sum()
            positions = np.flatnonzero(finite)
            position_by_column = dict(zip(positions.tolist(), ranks.tolist()))
            for target_index, column_index in enumerate(target_indices):
                if column_index in position_by_column:
                    output[row_index, target_index] = position_by_column[column_index]
        return pd.DataFrame(output, index=background.index, columns=self.target_features)

    def transform_external(
        self,
        background: pd.DataFrame,
        probe_map: Mapping[str, tuple[str, ...]],
    ) -> pd.DataFrame:
        """Rank canonical target probes against each external full probe row."""

        output = np.full((len(background), len(self.target_features)), np.nan, dtype=float)
        for row_index, row in enumerate(background.to_numpy(dtype=float)):
            finite = np.isfinite(row)
            if not finite.any():
                continue
            ranks = rankdata(row[finite], method="average") / finite.sum()
            ranked = pd.Series(np.nan, index=background.columns, dtype=float)
            ranked.iloc[np.flatnonzero(finite)] = ranks
            for feature_index, feature in enumerate(self.target_features):
                probe_values = pd.to_numeric(ranked.reindex(probe_map.get(feature, ())), errors="coerce").dropna()
                if not probe_values.empty:
                    output[row_index, feature_index] = float(probe_values.mean())
        return pd.DataFrame(output, index=background.index, columns=self.target_features)


def metrics(labels: np.ndarray, probabilities: np.ndarray) -> dict[str, float]:
    predictions = (probabilities >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(labels, predictions, labels=[0, 1]).ravel()
    return {
        "roc_auc": float(roc_auc_score(labels, probabilities)),
        "balanced_accuracy": float(balanced_accuracy_score(labels, predictions)),
        "sensitivity": float(tp / (tp + fn)) if tp + fn else np.nan,
        "specificity": float(tn / (tn + fp)) if tn + fp else np.nan,
    }


def model() -> Pipeline:
    return Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        ("model", LogisticRegression(C=0.3, class_weight="balanced", solver="liblinear", max_iter=5000, random_state=SEED)),
    ])


def load_local() -> tuple[pd.DataFrame, np.ndarray]:
    matrix = pd.read_csv(LOCAL_MATRIX)
    metadata = pd.read_csv(LOCAL_METADATA)
    frame = metadata[["sample_id_normalized", "diagnosis", "stage"]].merge(
        matrix, on="sample_id_normalized", how="inner", validate="one_to_one"
    )
    keep = frame["diagnosis"].eq("Cancer") & frame["stage"].eq("Stage I/II")
    keep |= frame["diagnosis"].isin(["Benign", "Control"])
    frame = frame.loc[keep].reset_index(drop=True)
    labels = (frame["diagnosis"].eq("Cancer") & frame["stage"].eq("Stage I/II")).astype(int).to_numpy()
    background = frame[[column for column in matrix.columns if column.startswith("hsa-")]].apply(pd.to_numeric, errors="coerce")
    return background, labels


def normalize(value: Any) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def feature_probe_map(gse: Any) -> dict[str, list[str]]:
    mapping = {feature: [] for feature in SIGNATURE}
    aliases = {normalize(feature): feature for feature in SIGNATURE}
    for feature in SIGNATURE:
        if feature == "hsa-miR-203a":
            aliases[normalize("hsa-miR-203a-3p")] = feature
    for platform in (getattr(gse, "gpls", {}) or {}).values():
        table = getattr(platform, "table", None)
        if not isinstance(table, pd.DataFrame) or table.empty or table.shape[1] == 0:
            continue
        id_column = next((column for column in table.columns if str(column).upper() in {"ID", "ID_REF", "PROBEID"}), table.columns[0])
        for _, row in table.iterrows():
            matched = [aliases[normalize(value)] for value in row.tolist() if normalize(value) in aliases]
            for feature in matched:
                mapping[feature].append(str(row[id_column]))
    return mapping


def numeric_values(table: pd.DataFrame) -> pd.Series:
    columns = [column for column in table.columns if str(column).upper() in {"VALUE", "VALUES", "EXPRESSION"}]
    if not columns:
        numeric_counts = table.apply(pd.to_numeric, errors="coerce").notna().sum()
        columns = [numeric_counts.idxmax()] if numeric_counts.max() else []
    if not columns:
        raise ValueError("No numeric GEO expression column")
    return pd.to_numeric(table[columns[0]], errors="coerce")


def load_external(accession: str) -> Cohort:
    path = GEO_CACHE / accession / f"{accession}_family.soft.gz"
    if not path.exists():
        raise FileNotFoundError(f"Missing cached cohort: {path}")
    gse = GEOparse.get_GEO(filepath=str(path), geotype="GSE", silent=True, how="full")
    probes = feature_probe_map(gse)
    rows: list[dict[str, float]] = []
    background_rows: list[dict[str, float]] = []
    labels: list[int | None] = []
    for sample in (getattr(gse, "gsms", {}) or {}).values():
        table = getattr(sample, "table", None)
        if not isinstance(table, pd.DataFrame) or table.empty or table.shape[1] == 0:
            continue
        identifiers = table.iloc[:, 0].astype(str)
        values = numeric_values(table)
        by_probe = pd.Series(values.to_numpy(), index=identifiers).groupby(level=0).mean()
        background_rows.append({str(key): float(value) for key, value in by_probe.items() if np.isfinite(value)})
        row = {}
        for feature in SIGNATURE:
            matched = pd.to_numeric(by_probe.reindex(probes[feature]), errors="coerce").dropna()
            row[feature] = float(matched.mean()) if not matched.empty else np.nan
        rows.append(row)
        metadata = getattr(sample, "metadata", {}) or {}
        text = " ".join(str(value) for key, value in metadata.items() if key in {"title", "source_name_ch1", "characteristics_ch1", "description"})
        labels.append(label_from_metadata(text))
    values = pd.DataFrame(rows, columns=SIGNATURE)
    background = pd.DataFrame(background_rows).replace([np.inf, -np.inf], np.nan)
    valid = np.array([label is not None for label in labels]) & values.notna().any(axis=1).to_numpy()
    labels_array = np.asarray([labels[index] for index, keep in enumerate(valid) if keep], dtype=int)
    values = values.loc[valid].reset_index(drop=True)
    background = background.loc[valid].reset_index(drop=True)
    if len(labels_array) == 0 or len(np.unique(labels_array)) != 2:
        raise ValueError(f"{accession} does not contain two evaluable label classes")
    return Cohort(
        accession,
        values,
        labels_array,
        background=background,
        probe_map={feature: tuple(str(probe) for probe in probes[feature]) for feature in SIGNATURE},
    )


def local_cross_validation(background: pd.DataFrame, labels: np.ndarray, transformer: WithinPatientBackgroundRankTransformer) -> tuple[dict[str, float], dict[str, float]]:
    raw_probabilities = np.full(len(labels), np.nan)
    rank_probabilities = np.full(len(labels), np.nan)
    splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    for train_index, validation_index in splitter.split(background, labels):
        train_raw = background.iloc[train_index]
        validation_raw = background.iloc[validation_index]
        train_rank = transformer.fit(train_raw).transform(train_raw)
        validation_rank = transformer.fit(train_raw).transform(validation_raw)
        raw_model = model().fit(train_raw[list(SIGNATURE)], labels[train_index])
        rank_model = model().fit(train_rank, labels[train_index])
        raw_probabilities[validation_index] = raw_model.predict_proba(validation_raw[list(SIGNATURE)])[:, 1]
        rank_probabilities[validation_index] = rank_model.predict_proba(validation_rank)[:, 1]
    return metrics(labels, raw_probabilities), metrics(labels, rank_probabilities)


def fit_and_evaluate(background: pd.DataFrame, labels: np.ndarray, transformer: WithinPatientBackgroundRankTransformer, train_rank_model: Pipeline, train_raw_model: Pipeline) -> tuple[dict[str, float], dict[str, float]]:
    rank_values = transformer.fit(background).transform(background)
    raw_values = background[list(SIGNATURE)]
    rank_model = train_rank_model.fit(transformer.transform(background), labels)
    raw_model = train_raw_model.fit(raw_values, labels)
    return raw_model, rank_model


def run() -> None:
    local_background, local_labels = load_local()
    transformer = WithinPatientBackgroundRankTransformer(SIGNATURE)
    local_raw, local_rank = local_cross_validation(local_background, local_labels, transformer)
    local_rank_values = transformer.fit(local_background).transform(local_background)
    raw_model = model().fit(local_background[list(SIGNATURE)], local_labels)
    rank_model = model().fit(local_rank_values, local_labels)
    rows: list[dict[str, object]] = []
    for name, raw_result, rank_result in [("GSE94533 (5-fold CV)", local_raw, local_rank)]:
        for space, result in [("raw", raw_result), ("rank", rank_result)]:
            rows.append({"cohort": name, "space": space, **result})
    for accession in EXTERNAL_ACCESSIONS:
        try:
            cohort = load_external(accession)
        except ValueError as exc:
            LOGGER.warning("Skipping %s: %s", accession, exc)
            rows.append({
                "cohort": accession,
                "space": "raw/rank",
                "roc_auc": np.nan,
                "balanced_accuracy": np.nan,
                "sensitivity": np.nan,
                "specificity": np.nan,
                "detail": f"Not evaluable: {exc}",
            })
            continue
        if cohort.background is None or cohort.probe_map is None:
            raise RuntimeError(f"{accession} has no external background matrix")
        external_rank = transformer.transform_external(cohort.background, cohort.probe_map)
        raw_result = metrics(cohort.labels, raw_model.predict_proba(cohort.values[list(SIGNATURE)])[:, 1])
        rank_result = metrics(cohort.labels, rank_model.predict_proba(external_rank)[:, 1])
        rows.extend([
            {"cohort": accession, "space": "raw", **raw_result},
            {"cohort": accession, "space": "rank", **rank_result},
        ])
    result_frame = pd.DataFrame(rows)
    print(result_frame.to_string(index=False, float_format=lambda value: f"{value:.4f}"))
    lines = [
        "# Within-Patient Rank-Transform Experiment", "",
        f"Frozen signature: `{', '.join(SIGNATURE)}`.",
        "Ranks are computed independently within each sample against that sample's full available hsa-miRNA background. Age and CA-125 are excluded.", "",
        "| Cohort | Space | ROC-AUC | Balanced Accuracy | Sensitivity | Specificity | Detail |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        def display_metric(value: object) -> str:
            return f"{value:.4f}" if isinstance(value, (float, int)) and np.isfinite(value) else "-"
        detail = str(row.get("detail", "")).replace("|", "\\|")
        lines.append(
            f"| {row['cohort']} | {row['space']} | {display_metric(row['roc_auc'])} | "
            f"{display_metric(row['balanced_accuracy'])} | {display_metric(row['sensitivity'])} | "
            f"{display_metric(row['specificity'])} | {detail} |"
        )
    REPORT.write_text("\n".join(lines) + "\n")
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