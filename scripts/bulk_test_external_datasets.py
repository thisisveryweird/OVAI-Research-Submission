#!/usr/bin/env python3
"""Download and externally test GEO ovarian-miRNA candidate series.

The runner is intentionally conservative: external labels are used only for
final evaluation, never for feature selection, preprocessing, threshold
tuning, or model fitting. GEO metadata and matrices are untrusted inputs and
each accession is isolated so one malformed series cannot stop the batch.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import re
import shutil
import urllib.request
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, confusion_matrix, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=pd.errors.DtypeWarning)
warnings.filterwarnings(
    "ignore",
    message=r"Inconsistent values: penalty=l1 with l1_ratio=0\.0.*",
    category=UserWarning,
)

try:
    import GEOparse
except ImportError as exc:  # pragma: no cover - exercised in deployment
    GEOparse = None
    GEOPARSE_IMPORT_ERROR = exc


PROJECT = Path(__file__).resolve().parents[1]
REPORT_INPUT = PROJECT / "reports" / "potential_external_datasets.md"
REPORT_OUTPUT = PROJECT / "reports" / "bulk_external_test_results.md"
MATRIX_PATH = PROJECT / "data" / "processed" / "gse94533_mirna_tpm.csv.gz"
METADATA_PATH = PROJECT / "data" / "processed" / "gse94533_metadata.csv"
DEFAULT_CACHE = PROJECT / "data" / "raw" / "geo"
DEFAULT_WORKERS = 3
SIGNATURE_SIZE = 5
SEED = 20260919
LOGGER = logging.getLogger("ovai.bulk_geo_test")


@dataclass
class Result:
    accession: str
    status: str
    recommendation: str
    sample_count: int | None = None
    labeled_count: int | None = None
    roc_auc: float | None = None
    balanced_accuracy: float | None = None
    sensitivity: float | None = None
    specificity: float | None = None
    missing_features: tuple[str, ...] = ()
    detail: str = ""


@dataclass(frozen=True)
class Candidate:
    accession: str
    sample_count: int | None


def read_accessions(path: Path) -> list[str]:
    """Extract and deduplicate GSE accessions from the scout report."""

    if not path.exists():
        raise FileNotFoundError(f"Candidate report not found: {path}")
    accessions = re.findall(r"^\|\s*(GSE\d+)\s*\|", path.read_text(), re.MULTILINE | re.IGNORECASE)
    unique = list(dict.fromkeys(accession.upper() for accession in accessions))
    if not unique:
        raise RuntimeError(f"No GSE accessions found in {path}")
    LOGGER.info("Parsed %d candidate accessions from %s", len(unique), path)
    return unique


def read_candidates(path: Path) -> list[Candidate]:
    """Read accession/sample-count pairs from the scout report."""

    candidates: list[Candidate] = []
    for line in path.read_text().splitlines():
        if not line.startswith("|") or "Accession ID (GSE)" in line or line.startswith("|---"):
            continue
        columns = [column.strip() for column in line.strip("|").split("|")]
        if len(columns) < 3 or not re.fullmatch(r"GSE\d+", columns[0], re.IGNORECASE):
            continue
        sample_count = int(columns[2]) if columns[2].isdigit() else None
        candidates.append(Candidate(columns[0].upper(), sample_count))
    if not candidates:
        raise RuntimeError(f"No candidate rows found in {path}")
    return list({candidate.accession: candidate for candidate in candidates}.values())


def text_value(value: Any) -> str:
    if isinstance(value, (list, tuple, set)):
        return " ".join(text_value(item) for item in value)
    if isinstance(value, Mapping):
        return " ".join(f"{key} {text_value(item)}" for key, item in value.items())
    return str(value).strip()


def normalize_feature(value: Any) -> str:
    return re.sub(r"[^a-z0-9]", "", text_value(value).lower())


def aliases(feature: str) -> set[str]:
    normalized = normalize_feature(feature)
    names = {normalized}
    if feature == "hsa-miR-203a":
        names.add(normalize_feature("hsa-miR-203a-3p"))
    return names


def sample_metadata(sample: Any) -> str:
    metadata = getattr(sample, "metadata", {}) or {}
    values: list[str] = []
    for key, value in metadata.items():
        if key in {"title", "source_name_ch1", "characteristics_ch1", "description"}:
            values.append(text_value(value))
    return " ".join(values)


def label_from_metadata(metadata: str) -> int | None:
    """Return the prespecified early-stage-compatible binary label."""

    value = metadata.lower()
    if "borderline" in value:
        return None
    ovarian_cancer = bool(re.search(r"ovarian\s+(cancer|carcinoma)|ovarian\s+tumou?r", value))
    if ovarian_cancer:
        return 1
    if re.search(r"\b(benign|healthy|control|non[- ]?cancer|normal)\b", value):
        return 0
    return None


def numeric_series(table: pd.DataFrame) -> pd.Series:
    if table is None or table.empty or table.shape[1] == 0:
        raise ValueError("GEO sample table is empty")
    value_columns = [column for column in table.columns if str(column).upper() in {"VALUE", "VALUES", "EXPRESSION"}]
    if not value_columns:
        numeric = table.apply(pd.to_numeric, errors="coerce").notna().sum()
        value_columns = [str(numeric.idxmax())] if numeric.max() else []
    if not value_columns:
        raise ValueError("No numeric expression column found in GEO sample table")
    return pd.to_numeric(table[value_columns[0]], errors="coerce")


def metadata_expression_values(sample: Any, features: tuple[str, ...]) -> dict[str, float]:
    """Extract feature/value pairs from non-tabular GSM metadata when present."""

    aliases_to_feature = {
        alias: feature for feature in features for alias in aliases(feature)
    }
    values: dict[str, float] = {}

    def visit(value: Any, key_hint: str = "") -> None:
        if isinstance(value, Mapping):
            for key, child in value.items():
                visit(child, normalize_feature(key))
            return
        if isinstance(value, (list, tuple, np.ndarray, pd.Series)):
            for child in value:
                visit(child, key_hint)
            return
        if not key_hint:
            return
        matched = next(
            (feature for alias, feature in aliases_to_feature.items() if alias in key_hint),
            None,
        )
        if matched is None:
            return
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            return
        if np.isfinite(numeric):
            values[matched] = numeric

    visit(getattr(sample, "metadata", {}) or {})
    return values


def platform_feature_map(gse: Any, features: Iterable[str]) -> dict[str, list[Any]]:
    """Map canonical miRNAs to platform probe IDs using GPL annotations."""

    targets = {feature: aliases(feature) for feature in features}
    mapping: dict[str, list[Any]] = {feature: [] for feature in targets}
    for platform in (getattr(gse, "gpls", {}) or {}).values():
        table = getattr(platform, "table", None)
        if not isinstance(table, pd.DataFrame) or table.empty:
            continue
        id_column = next((column for column in table.columns if str(column).upper() in {"ID", "ID_REF", "PROBEID"}), table.columns[0])
        for _, row in table.iterrows():
            row_text = {normalize_feature(value) for value in row.tolist() if pd.notna(value)}
            for feature, target_aliases in targets.items():
                if row_text & target_aliases:
                    mapping[feature].append(str(row[id_column]))
    return mapping


def expression_frame(gse: Any, features: tuple[str, ...]) -> tuple[pd.DataFrame, list[str]]:
    mapping = platform_feature_map(gse, features)
    sample_rows: list[dict[str, float]] = []
    labels: list[int | None] = []
    sample_ids: list[str] = []
    for sample_id, sample in (getattr(gse, "gsms", {}) or {}).items():
        table = getattr(sample, "table", None)
        row: dict[str, float]
        if not isinstance(table, pd.DataFrame) or table.empty or table.shape[1] == 0:
            LOGGER.warning(
                "[%s] GSM table is empty or has no columns; trying metadata expression fallback",
                sample_id,
            )
            row = {feature: np.nan for feature in features}
            row.update(metadata_expression_values(sample, features))
        else:
            ids = table.iloc[:, 0].astype(str)
            normalized_ids = {normalize_feature(value): value for value in ids}
            for feature in features:
                if not mapping[feature]:
                    mapping[feature] = [normalized_ids[name] for name in aliases(feature) if name in normalized_ids]
            values = numeric_series(table)
            by_probe = pd.Series(values.to_numpy(), index=ids).groupby(level=0).mean()
            row = {}
            for feature in features:
                probe_values = pd.to_numeric(by_probe.reindex(mapping[feature]), errors="coerce").dropna()
                row[feature] = float(probe_values.mean()) if not probe_values.empty else np.nan
        sample_rows.append(row)
        labels.append(label_from_metadata(sample_metadata(sample)))
        sample_ids.append(str(sample_id))
    frame = pd.DataFrame(sample_rows, index=sample_ids)
    return frame, labels


def load_local_training() -> tuple[pd.DataFrame, np.ndarray, tuple[str, ...]]:
    if not MATRIX_PATH.exists() or not METADATA_PATH.exists():
        raise FileNotFoundError(
            "GSE94533 processed inputs are required: "
            f"{MATRIX_PATH} and {METADATA_PATH}"
        )
    matrix = pd.read_csv(MATRIX_PATH)
    metadata = pd.read_csv(METADATA_PATH)
    required = {"sample_id_normalized", "diagnosis", "stage"}
    missing = required - set(metadata.columns)
    if missing:
        raise ValueError(f"Local metadata missing columns: {sorted(missing)}")
    features = tuple(column for column in matrix.columns if column.startswith("hsa-"))
    frame = metadata[["sample_id_normalized", "diagnosis", "stage"]].merge(
        matrix, on="sample_id_normalized", how="inner", validate="one_to_one"
    )
    keep = frame["diagnosis"].eq("Cancer") & frame["stage"].eq("Stage I/II")
    keep |= frame["diagnosis"].isin(["Benign", "Control"])
    frame = frame.loc[keep].reset_index(drop=True)
    labels = (frame["diagnosis"].eq("Cancer") & frame["stage"].eq("Stage I/II")).astype(int).to_numpy()
    if len(np.unique(labels)) != 2 or len(features) < SIGNATURE_SIZE:
        raise ValueError("Local GSE94533 training cohort does not contain the expected endpoint/features")
    return frame, labels, features


def fit_model(frame: pd.DataFrame, labels: np.ndarray, features: tuple[str, ...]) -> tuple[Pipeline, tuple[str, ...]]:
    values = np.log1p(frame[list(features)].clip(lower=0).to_numpy(dtype=float))
    preprocessor = Pipeline([("impute", SimpleImputer(strategy="median")), ("scale", StandardScaler())])
    training_values = preprocessor.fit_transform(values)
    selector = LogisticRegression(penalty="l1", solver="liblinear", C=0.1, class_weight="balanced", max_iter=5000, random_state=SEED)
    selector.fit(training_values, labels)
    ranking = np.argsort(-np.abs(selector.coef_[0]), kind="stable")[:SIGNATURE_SIZE]
    signature = tuple(features[index] for index in ranking)
    pipeline = Pipeline([("impute", SimpleImputer(strategy="median")), ("scale", StandardScaler()), ("model", LogisticRegression(solver="liblinear", class_weight="balanced", max_iter=5000, random_state=SEED))])
    search = GridSearchCV(pipeline, {"model__C": [0.03, 0.1, 0.3, 1.0]}, scoring="balanced_accuracy", cv=StratifiedKFold(3, shuffle=True, random_state=SEED), n_jobs=-1, error_score="raise")
    search.fit(np.log1p(frame[list(signature)].clip(lower=0)), labels)
    LOGGER.info("Frozen local signature: %s", ", ".join(signature))
    return search.best_estimator_, signature


def evaluate(model: Pipeline, frame: pd.DataFrame, labels: np.ndarray, signature: tuple[str, ...]) -> dict[str, float]:
    values = np.log1p(frame[list(signature)].clip(lower=0))
    probability = model.predict_proba(values)[:, 1]
    prediction = (probability >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(labels, prediction, labels=[0, 1]).ravel()
    return {
        "roc_auc": float(roc_auc_score(labels, probability)),
        "balanced_accuracy": float(balanced_accuracy_score(labels, prediction)),
        "sensitivity": float(tp / (tp + fn)) if tp + fn else np.nan,
        "specificity": float(tn / (tn + fp)) if tn + fp else np.nan,
    }


def geo_family_url(accession: str) -> str:
    range_subdir = re.sub(r"\d{1,3}$", "nnn", accession)
    return (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{range_subdir}/{accession}/soft/{accession}_family.soft.gz"
    )


def download_one(accession: str, cache: Path) -> Path:
    """Download one GSE family file atomically, reusing completed cache files."""

    destination = cache / accession / f"{accession}_family.soft.gz"
    if destination.exists() and destination.stat().st_size > 0:
        LOGGER.info("[%s] using cached GEO file", accession)
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    temporary.unlink(missing_ok=True)
    LOGGER.info("[%s] downloading GEO family file", accession)
    request = urllib.request.Request(
        geo_family_url(accession),
        headers={"User-Agent": "ovai-bulk-geo-test/1.0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=180) as response, temporary.open("wb") as output:
            shutil.copyfileobj(response, output)
        if temporary.stat().st_size == 0:
            raise RuntimeError("downloaded file is empty")
        with gzip.open(temporary, "rb") as compressed:
            compressed.read(1)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


def prefetch(accessions: list[str], cache: Path, workers: int) -> dict[str, Path | Exception]:
    """Prefetch GSE files concurrently before any model evaluation starts."""

    downloaded: dict[str, Path | Exception] = {}
    with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="geo-download") as pool:
        futures = {pool.submit(download_one, accession, cache): accession for accession in accessions}
        for future in as_completed(futures):
            accession = futures[future]
            try:
                downloaded[accession] = future.result()
                LOGGER.info("[%s] download ready", accession)
            except Exception as exc:
                downloaded[accession] = exc
                LOGGER.error("[%s] download failed: %s", accession, exc)
    return downloaded


def test_one(accession: str, model: Pipeline, signature: tuple[str, ...], source: Path | Exception) -> Result:
    try:
        if GEOparse is None:
            raise RuntimeError(f"GEOparse is unavailable: {GEOPARSE_IMPORT_ERROR}")
        if isinstance(source, Exception):
            raise RuntimeError(f"download failed: {source}") from source
        LOGGER.info("[%s] parsing cached GEO family file", accession)
        gse = GEOparse.get_GEO(filepath=str(source), geotype="GSE", silent=True, how="full")
        values, raw_labels = expression_frame(gse, signature)
        missing = tuple(feature for feature in signature if values[feature].notna().sum() == 0)
        if missing:
            return Result(accession, "Skipped: Missing Features", "Skipped: Missing Features", len(values), detail="Missing: " + ", ".join(missing), missing_features=missing)
        keep = np.array([label is not None for label in raw_labels]) & values[list(signature)].notna().any(axis=1).to_numpy()
        labels = np.array([raw_labels[index] for index, value in enumerate(keep) if value], dtype=int)
        values = values.loc[keep].copy()
        if len(labels) != len(values) or len(np.unique(labels)) != 2:
            return Result(accession, "Failed: Missing Labels", "Failed: Missing Labels", len(values), detail="Need both ovarian cancer and non-cancer labeled samples")
        metrics = evaluate(model, values, labels, signature)
        success = metrics["sensitivity"] > 0.50 and metrics["balanced_accuracy"] > 0.60
        recommendation = "Approved for Integration" if success else "Failed: Platform Shift / 0% Sensitivity"
        status = "SUCCESSFULLY COMPATIBLE" if success else "FAILED"
        return Result(accession, status, recommendation, len(values), len(labels), **metrics)
    except Exception as exc:
        LOGGER.exception("[%s] failed", accession)
        return Result(accession, "Failed: Download or Parse Error", "Failed: Download or Parse Error", detail=str(exc))


def write_report(path: Path, results: list[Result], signature: tuple[str, ...]) -> None:
    priority = {"SUCCESSFULLY COMPATIBLE": 0, "FAILED": 1, "Skipped: Missing Features": 2}
    ordered = sorted(results, key=lambda result: (priority.get(result.status, 3), result.accession))
    lines = [
        "# Bulk External GEO Dataset Test Results", "",
        "Evaluation uses the complete local GSE94533 Stage I/II cancer versus benign/control cohort. The model is a regularized logistic regression with local median imputation and scaling; external labels are never used for fitting or tuning.", "",
        f"- Candidates tested: {len(results)}", f"- Frozen five-miRNA signature: `{', '.join(signature)}`", "- Success rule: sensitivity > 50% and balanced accuracy > 60%.", "",
        "| GSE | Status | Recommendation | N evaluated | ROC-AUC | Balanced Accuracy | Sensitivity | Specificity | Detail |",
        "|---|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for result in ordered:
        def metric(value: float | None) -> str:
            return f"{value:.3f}" if value is not None and np.isfinite(value) else "-"
        detail = result.detail.replace("|", "\\|").replace("\n", " ")
        lines.append(f"| {result.accession} | {result.status} | {result.recommendation} | {result.labeled_count or '-'} | {metric(result.roc_auc)} | {metric(result.balanced_accuracy)} | {metric(result.sensitivity)} | {metric(result.specificity)} | {detail} |")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-report", type=Path, default=REPORT_INPUT)
    parser.add_argument("--output-report", type=Path, default=REPORT_OUTPUT)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--top-n", type=int, default=10, help="Test the top N candidates by reported sample count.")
    parser.add_argument("--limit", type=int, help="Deprecated alias for --top-n.")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS, help="Concurrent GEO download workers (default: 3).")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    if args.top_n < 1:
        raise ValueError("--top-n must be positive")
    if args.limit is not None:
        if args.limit < 1:
            raise ValueError("--limit must be positive")
        args.top_n = args.limit
    if args.workers < 1 or args.workers > 3:
        raise ValueError("--workers must be between 1 and 3 to respect NCBI request limits")
    candidates = sorted(
        read_candidates(args.input_report),
        key=lambda candidate: (candidate.sample_count is not None, candidate.sample_count or -1),
        reverse=True,
    )
    selected = candidates[:args.top_n]
    accessions = [candidate.accession for candidate in selected]
    LOGGER.info("Selected top %d candidates by reported sample count: %s", len(accessions), ", ".join(accessions))
    frame, labels, features = load_local_training()
    model, signature = fit_model(frame, labels, features)
    downloaded = prefetch(accessions, args.cache_dir, args.workers)
    results: list[Result] = []
    for index, accession in enumerate(accessions, start=1):
        LOGGER.info("[%d/%d] testing %s", index, len(accessions), accession)
        result = test_one(accession, model, signature, downloaded[accession])
        results.append(result)
    write_report(args.output_report, results, signature)
    LOGGER.info("Wrote %s", args.output_report)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        LOGGER.error("Bulk GEO test failed: %s", exc)
        raise SystemExit(2) from exc