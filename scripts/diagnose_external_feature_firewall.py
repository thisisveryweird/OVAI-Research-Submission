#!/usr/bin/env python3
"""Audit cached GEO cohorts against the frozen five-miRNA firewall.

This diagnostic does not fit a model and never uses external labels. It only
records whether the frozen features and a broad expression matrix are present,
which distinguishes direct feature salvage from a required rank-transform
experiment.
"""

from __future__ import annotations

import argparse
import logging
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import GEOparse
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=pd.errors.DtypeWarning)

PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_CACHE = PROJECT / "data" / "raw" / "geo"
DEFAULT_QUEUE = PROJECT / "reports" / "potential_external_datasets.md"
DEFAULT_OUTPUT = PROJECT / "reports" / "external_feature_firewall_diagnostic.md"
SIGNATURE = (
    "hsa-miR-200b-3p",
    "hsa-miR-296-3p",
    "hsa-miR-6826-5p",
    "hsa-miR-32-3p",
    "hsa-miR-1255a",
)
LOGGER = logging.getLogger("ovai.feature_firewall")


@dataclass(frozen=True)
class Audit:
    accession: str
    cached: bool
    samples: int = 0
    nonempty_samples: int = 0
    feature_count: int = 0
    present: tuple[str, ...] = ()
    missing: tuple[str, ...] = SIGNATURE
    decision: str = "Not cached"
    detail: str = ""


def normalize(value: Any) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def read_queue(path: Path) -> list[str]:
    pattern = re.compile(r"^\|\s*(GSE\d+)\s*\|", re.IGNORECASE)
    accessions = [
        match.group(1).upper()
        for line in path.read_text().splitlines()
        if (match := pattern.match(line))
    ]
    return list(dict.fromkeys(accessions))


def cached_path(cache: Path, accession: str) -> Path:
    return cache / accession / f"{accession}_family.soft.gz"


def expression_identifiers(gse: Any) -> tuple[set[str], int, int, int]:
    identifiers: set[str] = set()
    nonempty = 0
    row_counts: list[int] = []
    for sample in (getattr(gse, "gsms", {}) or {}).values():
        table = getattr(sample, "table", None)
        if not isinstance(table, pd.DataFrame) or table.empty or table.shape[1] == 0:
            continue
        nonempty += 1
        row_counts.append(table.shape[0])
        identifiers.update(normalize(value) for value in table.iloc[:, 0].astype(str))
    median_rows = int(np.median(row_counts)) if row_counts else 0
    return identifiers, nonempty, median_rows, max(row_counts, default=0)


def platform_identifiers(gse: Any) -> tuple[set[str], int]:
    identifiers: set[str] = set()
    feature_count = 0
    for platform in (getattr(gse, "gpls", {}) or {}).values():
        table = getattr(platform, "table", None)
        if not isinstance(table, pd.DataFrame) or table.empty or table.shape[1] == 0:
            continue
        feature_count = max(feature_count, len(table))
        for value in table.astype(str).to_numpy().ravel():
            identifiers.add(normalize(value))
    return identifiers, feature_count


def audit_one(accession: str, cache: Path) -> Audit:
    path = cached_path(cache, accession)
    if not path.exists() or path.stat().st_size == 0:
        return Audit(accession=accession, cached=False)
    try:
        LOGGER.info("Auditing %s", accession)
        gse = GEOparse.get_GEO(filepath=str(path), geotype="GSE", silent=True, how="full")
        sample_ids, nonempty, median_rows, max_rows = expression_identifiers(gse)
        platform_ids, platform_count = platform_identifiers(gse)
        available = sample_ids | platform_ids
        present = tuple(feature for feature in SIGNATURE if any(normalize(feature) in identifier for identifier in available))
        missing = tuple(feature for feature in SIGNATURE if feature not in present)
        feature_count = max(platform_count, median_rows, max_rows)
        if not nonempty:
            decision = "No expression matrix; cannot salvage raw features"
        elif not missing and feature_count >= 100:
            decision = "Broad panel with exact signature; prioritize within-patient rank transform"
        elif not missing:
            decision = "Exact signature present; direct feature salvage is plausible"
        elif feature_count >= 100:
            decision = "Broad panel but incomplete signature; rank transform may be salvageable"
        else:
            decision = "Incompatible Feature Space"
        return Audit(
            accession=accession,
            cached=True,
            samples=len(getattr(gse, "gsms", {}) or {}),
            nonempty_samples=nonempty,
            feature_count=feature_count,
            present=present,
            missing=missing,
            decision=decision,
            detail=f"median expression rows/sample={median_rows}; platform rows={platform_count}",
        )
    except Exception as exc:
        LOGGER.exception("Could not audit %s", accession)
        return Audit(accession=accession, cached=True, decision="Parse error", detail=str(exc))


def write_report(path: Path, audits: list[Audit]) -> None:
    lines = [
        "# External Feature Validation Firewall Diagnostic",
        "",
        f"Frozen signature: `{', '.join(SIGNATURE)}`.",
        "This is a structural cache audit only. It does not establish assay equivalence or clinical validity.",
        "",
        "| GSE | Cached | Samples | Non-empty GSMs | Matrix features | Exact signature | Missing features | Decision | Detail |",
        "|---|---:|---:|---:|---:|---|---|---|---|",
    ]
    for audit in audits:
        present = ", ".join(audit.present) or "None"
        missing = ", ".join(audit.missing) or "None"
        detail = audit.detail.replace("|", "\\|").replace("\n", " ")
        lines.append(
            f"| {audit.accession} | {'yes' if audit.cached else 'no'} | {audit.samples or '-'} | "
            f"{audit.nonempty_samples or '-'} | {audit.feature_count or '-'} | {present} | {missing} | "
            f"{audit.decision} | {detail} |"
        )
    lines.extend([
        "",
        "## Interpretation",
        "",
        "- A broad panel with all five frozen features is not a restricted-subset failure; it is a candidate for a prespecified within-patient rank-transform experiment because rank order can reduce platform-scale dependence.",
        "- A cohort without an expression table cannot be salvaged from its family SOFT file alone; obtain its series matrix or processed supplementary matrix first.",
        "- Missing frozen features block direct Phase 8 evaluation. A broad alternative panel may still support a separate rank-based model, but it must not be silently substituted into the frozen signature.",
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--queue", type=Path, default=DEFAULT_QUEUE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    audits = [audit_one(accession, args.cache_dir) for accession in read_queue(args.queue)]
    write_report(args.output, audits)
    LOGGER.info("Wrote %s", args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())