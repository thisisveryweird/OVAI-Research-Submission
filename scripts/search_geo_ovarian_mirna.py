#!/usr/bin/env python3
"""Scout GEO for candidate blood-based ovarian-cancer miRNA cohorts.

This script performs metadata discovery only. It does not download expression
matrices, infer patient overlap, or certify a study as a valid external
validation cohort. Every retained record requires manual follow-up before use.
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from Bio import Entrez
from Bio.Entrez.Parser import DictionaryElement, ListElement
from urllib.error import HTTPError, URLError


PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = PROJECT / "reports" / "potential_external_datasets.md"
DEFAULT_TOOL = "ovai-geo-ovarian-mirna-scout"
DEFAULT_RETMAX = 100
NO_TISSUE_TERMS = (
    "tissue",
    "tumor sample",
    "ovarian tissue",
    "ffpe",
    "cell line",
    "xenograft",
)
QUERY = (
    '("ovarian cancer" OR "ovarian carcinoma") AND '
    '(microrna OR mirna OR "small rna") AND '
    '(serum OR plasma OR blood OR "liquid biopsy")'
)

LOGGER = logging.getLogger("ovai.geo_scout")


@dataclass(frozen=True)
class GeoCandidate:
    accession: str
    platform: str
    sample_count: int | None
    title: str
    abstract_summary: str
    sample_attributes: tuple[str, ...]
    tissue_filter_text: str


class RateLimiter:
    """Enforce a conservative NCBI request interval."""

    def __init__(self, interval_seconds: float) -> None:
        self.interval_seconds = interval_seconds
        self._last_request = 0.0

    def wait(self) -> None:
        elapsed = time.monotonic() - self._last_request
        if elapsed < self.interval_seconds:
            time.sleep(self.interval_seconds - elapsed)
        self._last_request = time.monotonic()


def configure_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )


def configure_entrez(email: str | None, api_key: str | None, tool: str) -> None:
    """Configure Entrez contact information and refuse anonymous misuse."""

    resolved_email = email or os.environ.get("NCBI_EMAIL")
    if not resolved_email:
        raise ValueError(
            "NCBI contact email is required. Pass --email or set NCBI_EMAIL."
        )
    Entrez.email = resolved_email
    Entrez.tool = tool
    Entrez.api_key = api_key or os.environ.get("NCBI_API_KEY")
    LOGGER.info(
        "Configured Entrez with tool=%s api_key=%s",
        tool,
        "yes" if Entrez.api_key else "no; using <=3 requests/second",
    )


def as_text(value: Any) -> str:
    """Convert Entrez parser values to normalized text."""

    if value is None:
        return ""
    if isinstance(value, (list, tuple, ListElement)):
        return " ".join(as_text(item) for item in value)
    if isinstance(value, (dict, DictionaryElement)):
        return " ".join(f"{key} {as_text(item)}" for key, item in value.items())
    return str(value).strip()


def first_field(record: Mapping[str, Any], names: Sequence[str]) -> str:
    """Return the first non-empty field under any known GEO schema name."""

    lowered = {str(key).lower(): value for key, value in record.items()}
    for name in names:
        value = lowered.get(name.lower())
        if value not in (None, "", [], ()):
            return as_text(value)
    return ""


def request_entrez(
    operation: str,
    limiter: RateLimiter,
    retries: int = 4,
    **kwargs: Any,
) -> Any:
    """Run one Entrez request with throttling and transient-error retries."""

    delay = 2.0
    for attempt in range(1, retries + 1):
        limiter.wait()
        try:
            LOGGER.debug("Entrez %s attempt=%d kwargs=%s", operation, attempt, kwargs)
            with getattr(Entrez, operation)(**kwargs) as handle:
                return Entrez.read(handle)
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            if attempt == retries:
                LOGGER.error("Entrez %s failed after %d attempts: %s", operation, attempt, exc)
                raise
            LOGGER.warning("Entrez %s transient failure (%s); retrying in %.1fs", operation, exc, delay)
            time.sleep(delay)
            delay *= 2
        except Exception as exc:
            LOGGER.error("Entrez %s returned an unexpected error: %s", operation, exc)
            raise
    raise RuntimeError(f"Unreachable Entrez retry state for {operation}")


def search_geo(query: str, limiter: RateLimiter, retmax: int) -> list[str]:
    """Search GEO series records and return Entrez UIDs."""

    result = request_entrez(
        "esearch",
        limiter,
        db="gds",
        term=f"{query} AND gse[Entry Type]",
        retmode="xml",
        retmax=retmax,
        sort="relevance",
    )
    ids = [str(item) for item in result.get("IdList", [])]
    LOGGER.info("Entrez returned %d GEO series UIDs", len(ids))
    return ids


def fetch_summaries(uids: Sequence[str], limiter: RateLimiter) -> list[Mapping[str, Any]]:
    """Fetch GEO summaries in one request, preserving UID-level records."""

    if not uids:
        return []
    result = request_entrez(
        "esummary",
        limiter,
        db="gds",
        id=", ".join(uids),
        retmode="xml",
    )
    documents = result.get("DocumentSummarySet", {}).get("DocumentSummary", [])
    return [dict(document) for document in documents]


def sample_attributes(record: Mapping[str, Any]) -> tuple[str, ...]:
    """Extract sample-attribute text from common ESummary field shapes."""

    raw = first_field(
        record,
        ("samples", "sample", "sampleattributes", "sample_attributes", "characteristics"),
    )
    values = tuple(part.strip() for part in re.split(r"\s*\|\s*|\n+", raw) if part.strip())
    return values


def parse_platform(record: Mapping[str, Any]) -> str:
    """Extract one or more GPL accessions from GEO summary fields."""

    raw = first_field(record, ("platformacc", "platform", "gpl", "platforms"))
    platforms = sorted(set(re.findall(r"GPL\d+", raw, flags=re.IGNORECASE)))
    return ", ".join(platforms) if platforms else (raw or "Unknown")


def parse_accession(record: Mapping[str, Any]) -> str:
    raw = first_field(record, ("accession", "accession_version", "gse"))
    match = re.search(r"GSE\d+", raw, flags=re.IGNORECASE)
    return match.group(0).upper() if match else raw or "Unknown"


def parse_sample_count(record: Mapping[str, Any], attributes: Sequence[str]) -> int | None:
    raw = first_field(record, ("n_samples", "samples", "sample_count", "n"))
    matches = re.findall(r"\b\d+\b", raw)
    if matches:
        try:
            return max(int(value) for value in matches)
        except ValueError:
            pass
    accession_like = re.findall(r"(?:GSM|SRS|SRR)\d+", " ".join(attributes), flags=re.IGNORECASE)
    return len(set(accession_like)) or None


def normalized_filter_text(record: Mapping[str, Any], attributes: Sequence[str]) -> str:
    fields = [
        first_field(record, ("title",)),
        first_field(record, ("summary", "description", "abstract")),
        first_field(record, ("design", "type", "source")),
        " ".join(attributes),
    ]
    return " ".join(fields).lower()


def is_tissue_contaminated(text: str) -> str | None:
    for term in NO_TISSUE_TERMS:
        if term in text:
            return term
    return None


def candidate_from_summary(record: Mapping[str, Any]) -> tuple[GeoCandidate | None, str | None]:
    attributes = sample_attributes(record)
    filter_text = normalized_filter_text(record, attributes)
    matched_term = is_tissue_contaminated(filter_text)
    if matched_term:
        return None, matched_term
    candidate = GeoCandidate(
        accession=parse_accession(record),
        platform=parse_platform(record),
        sample_count=parse_sample_count(record, attributes),
        title=first_field(record, ("title",)) or "Untitled GEO series",
        abstract_summary=first_field(record, ("summary", "description", "abstract")) or "No abstract summary available.",
        sample_attributes=attributes,
        tissue_filter_text=filter_text,
    )
    return candidate, None


def compact_summary(text: str, limit: int = 420) -> str:
    compact = re.sub(r"\s+", " ", text).strip()
    return compact if len(compact) <= limit else compact[: limit - 1].rstrip() + "…"


def markdown_escape(text: str) -> str:
    return text.replace("|", "\\|").replace("\n", " ").strip()


def write_report(
    path: Path,
    query: str,
    candidates: Sequence[GeoCandidate],
    filtered: Mapping[str, int],
    searched_count: int,
) -> None:
    """Write a readable candidate report with audit metadata."""

    lines = [
        "# Potential External Blood-Based Ovarian miRNA Datasets",
        "",
        "This is an automated GEO metadata scout, not a validation decision. Candidate studies require manual inspection of specimen type, platform compatibility, phenotype definitions, patient overlap, and raw/processed data provenance before download or modeling.",
        "",
        "## Search Audit",
        "",
        f"- Query: `{query}`",
        f"- GEO series records inspected: {searched_count}",
        f"- Candidates retained after tissue filtering: {len(candidates)}",
        f"- Records dropped by tissue keyword: {filtered.get('tissue_keyword', 0)}",
        f"- Records with malformed/unknown metadata: {filtered.get('metadata_error', 0)}",
        "- Rate policy: Entrez requests throttled to at most approximately 3 requests/second without an API key.",
        "",
        "## Retained Candidates",
        "",
        "| Accession ID (GSE) | Platform (GPL) | Sample Count (N) | Title | Study Abstract Summary |",
        "|---|---|---:|---|---|",
    ]
    if candidates:
        for candidate in candidates:
            count = str(candidate.sample_count) if candidate.sample_count is not None else "Unknown"
            lines.append(
                f"| {markdown_escape(candidate.accession)} | {markdown_escape(candidate.platform)} | "
                f"{count} | {markdown_escape(candidate.title)} | "
                f"{markdown_escape(compact_summary(candidate.abstract_summary))} |"
            )
    else:
        lines.append("| None retained | - | - | - | No records passed the metadata filter. |")
    lines.extend([
        "",
        "## Scientific Safeguards",
        "",
        "- Tissue-keyword filtering is conservative and can produce false positives or false negatives; inspect original GEO metadata before use.",
        "- A serum/plasma/blood keyword match does not prove that every sample is blood-derived.",
        "- Sample counts are estimates from GEO summary fields and sample accession-like metadata.",
        "- This scout does not download expression matrices and does not infer patient-level independence.",
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--email", help="NCBI contact email; otherwise NCBI_EMAIL is used.")
    parser.add_argument("--api-key", help="NCBI API key; otherwise NCBI_API_KEY is used.")
    parser.add_argument("--retmax", type=int, default=DEFAULT_RETMAX, help="Maximum GEO series to inspect.")
    parser.add_argument("--output", type=Path, default=DEFAULT_REPORT, help="Markdown report path.")
    parser.add_argument("--tool", default=DEFAULT_TOOL, help="Entrez tool identifier.")
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    configure_logging(args.verbose)
    if args.retmax < 1 or args.retmax > 1000:
        raise ValueError("--retmax must be between 1 and 1000.")
    configure_entrez(args.email, args.api_key, args.tool)
    interval = 0.11 if (args.api_key or os.environ.get("NCBI_API_KEY")) else 0.34
    limiter = RateLimiter(interval)
    LOGGER.info("Searching GEO with query: %s", QUERY)
    uids = search_geo(QUERY, limiter, args.retmax)
    summaries = fetch_summaries(uids, limiter)
    candidates: list[GeoCandidate] = []
    filtered = {"tissue_keyword": 0, "metadata_error": 0}
    for summary in summaries:
        try:
            candidate, dropped_term = candidate_from_summary(summary)
            if dropped_term:
                LOGGER.info("Dropping %s due to tissue keyword: %s", parse_accession(summary), dropped_term)
                filtered["tissue_keyword"] += 1
            elif candidate and candidate.accession != "Unknown":
                candidates.append(candidate)
                LOGGER.info("Retaining candidate %s platform=%s N=%s", candidate.accession, candidate.platform, candidate.sample_count)
            else:
                filtered["metadata_error"] += 1
        except Exception:
            LOGGER.exception("Could not parse one GEO summary")
            filtered["metadata_error"] += 1
    deduplicated = {candidate.accession: candidate for candidate in candidates}
    write_report(args.output, QUERY, sorted(deduplicated.values(), key=lambda item: item.accession), filtered, len(summaries))
    LOGGER.info("Wrote %s with %d retained candidates", args.output, len(deduplicated))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        LOGGER.error("GEO scout failed: %s", exc)
        raise SystemExit(2) from exc
