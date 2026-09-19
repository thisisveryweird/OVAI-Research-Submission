#!/usr/bin/env python3
"""Screen uncached GEO candidates against the repository inclusion laws.

This is a metadata triage step, not a substitute for inspecting raw matrices.
A candidate is approved for automated caching only when the report metadata
contains explicit sequencing, circulating/small-RNA, and ovarian cancer versus
control evidence. Array evidence always rejects a candidate.
"""

from __future__ import annotations

import argparse
import logging
import re
from dataclasses import dataclass
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_QUEUE = PROJECT / "reports" / "potential_external_datasets.md"
DEFAULT_CACHE = PROJECT / "data" / "raw" / "geo"
DEFAULT_OUTPUT = PROJECT / "reports" / "remaining_queue_screen.md"
LOGGER = logging.getLogger("ovai.remaining_queue_screen")

SEQUENCING_TERMS = (
    "sequencing",
    "rna-seq",
    "rna seq",
    "small rna-seq",
    "small rna seq",
    "mirna-seq",
    "mirna seq",
    "next-generation",
    "next generation",
    "illumina",
)
ARRAY_TERMS = (
    "microarray",
    "microarrays",
    "array",
    "agilent",
    "affymetrix",
    "toray",
    "3d-gene",
    "beadchip",
    "mas5",
    "plier",
    "gcrma",
)
CIRCULATING_TERMS = (
    "serum",
    "plasma",
    "blood",
    "circulating",
    "extracellular vesicle",
    "exosome",
    "liquid biopsy",
)
SMALL_RNA_TERMS = ("microrna", "mirna", "small rna", "small-rna")
CONTROL_TERMS = (
    "benign",
    "healthy",
    "control",
    "non-cancer",
    "noncancer",
    "normal",
)
STAGE_TERMS = (
    "stage i",
    "stage ii",
    "stage i/ii",
    "stage 1",
    "stage 2",
    "early stage",
    "early-stage",
)


@dataclass(frozen=True)
class Candidate:
    accession: str
    sample_count: int | None
    platform: str
    title: str
    abstract: str

    @property
    def text(self) -> str:
        return f"{self.platform} {self.title} {self.abstract}".lower()


@dataclass(frozen=True)
class Decision:
    candidate: Candidate
    profiling_pass: bool
    feature_pass: bool
    phenotype_pass: bool
    early_stage_documented: bool
    approved: bool
    reasons: tuple[str, ...]


def parse_candidates(path: Path) -> list[Candidate]:
    """Parse the scout Markdown table without relying on abstract pipes."""

    candidates: list[Candidate] = []
    pattern = re.compile(r"^\|\s*(GSE\d+)\s*\|\s*(.*?)\s*\|\s*(.*?)\s*\|\s*(.*?)\s*\|\s*(.*?)\s*\|$")
    for line in path.read_text().splitlines():
        match = pattern.match(line)
        if not match or match.group(1).lower() == "accession id (gse)":
            continue
        accession, platform, sample, title, abstract = match.groups()
        sample_count = int(sample) if sample.isdigit() else None
        candidates.append(Candidate(accession.upper(), sample_count, platform, title, abstract))
    if not candidates:
        raise RuntimeError(f"No candidate rows found in {path}")
    return list({candidate.accession: candidate for candidate in candidates}.values())


def cached_accessions(cache: Path) -> set[str]:
    return {
        path.parent.name.upper()
        for path in cache.glob("GSE*/GSE*_family.soft.gz")
        if path.stat().st_size > 0
    }


def contains_any(text: str, terms: tuple[str, ...]) -> bool:
    return any(term in text for term in terms)


def evaluate(candidate: Candidate) -> Decision:
    text = candidate.text
    has_array = contains_any(text, ARRAY_TERMS)
    has_sequencing = contains_any(text, SEQUENCING_TERMS)
    profiling_pass = has_sequencing and not has_array
    has_circulating = contains_any(text, CIRCULATING_TERMS)
    has_small_rna = contains_any(text, SMALL_RNA_TERMS)
    feature_pass = has_circulating and has_small_rna
    has_ovarian = bool(re.search(r"ovarian\s+(cancer|carcinoma|tumou?r|clear cell)", text))
    phenotype_pass = has_ovarian and contains_any(text, CONTROL_TERMS)
    early_stage = contains_any(text, STAGE_TERMS)
    reasons: list[str] = []
    if not profiling_pass:
        reasons.append("profiling law failed: missing explicit sequencing or array platform evidence")
    if has_array:
        reasons.append("array evidence overrides generic sequencing/profiling terms")
    if not feature_pass:
        reasons.append("feature law failed: missing circulating/blood small-RNA evidence")
    if not phenotype_pass:
        reasons.append("phenotype law failed: missing ovarian cancer plus benign/healthy/control evidence")
    if not early_stage:
        reasons.append("early-stage documentation not found")
    return Decision(candidate, profiling_pass, feature_pass, phenotype_pass, early_stage, profiling_pass and feature_pass and phenotype_pass, tuple(reasons))


def write_report(path: Path, decisions: list[Decision], source_count: int) -> None:
    passing = [decision for decision in decisions if decision.approved]
    lines = [
        "# Remaining GEO Queue Screen",
        "",
        "Metadata-only triage of uncached candidates from the potential external dataset report. Approval means suitable for automated download/cache triage, not scientific validation or model integration.",
        "",
        f"- Candidates in source report: {source_count}",
        f"- Uncached candidates screened: {len(decisions)}",
        f"- Candidates passing all laws: {len(passing)}",
        "",
        "## Approved Candidates",
        "",
        "| GEO Accession | Sample Size | Platform Type | Early Stage Documented | Approved for Automated Download/Cache |",
        "|---|---:|---|---|---|",
    ]
    if passing:
        for decision in sorted(passing, key=lambda item: (item.candidate.sample_count or 0), reverse=True):
            candidate = decision.candidate
            lines.append(f"| {candidate.accession} | {candidate.sample_count or 'Unknown'} | {candidate.platform} | {'Yes' if decision.early_stage_documented else 'No'} | Yes |")
    else:
        lines.append("| None | - | No remaining candidate passed all three inclusion laws | - | No |")
    lines.extend([
        "",
        "## Full Screening Decisions",
        "",
        "| GSE | N | Platform | Profiling | Feature | Phenotype | Early Stage | Decision | Reasons |",
        "|---|---:|---|---|---|---|---|---|---|",
    ])
    for decision in sorted(decisions, key=lambda item: item.candidate.accession):
        candidate = decision.candidate
        reasons = "; ".join(decision.reasons).replace("|", "\\|") or "All laws passed"
        lines.append(f"| {candidate.accession} | {candidate.sample_count or 'Unknown'} | {candidate.platform} | {'PASS' if decision.profiling_pass else 'FAIL'} | {'PASS' if decision.feature_pass else 'FAIL'} | {'PASS' if decision.phenotype_pass else 'FAIL'} | {'Yes' if decision.early_stage_documented else 'No'} | {'APPROVED' if decision.approved else 'REJECTED'} | {reasons} |")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--queue", type=Path, default=DEFAULT_QUEUE)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    candidates = parse_candidates(args.queue)
    cached = cached_accessions(args.cache_dir)
    remaining = [candidate for candidate in candidates if candidate.accession not in cached]
    decisions = [evaluate(candidate) for candidate in remaining]
    for decision in decisions:
        LOGGER.info("%s: %s", decision.candidate.accession, "APPROVED" if decision.approved else "REJECTED")
    write_report(args.output, decisions, len(candidates))
    passing = [decision for decision in decisions if decision.approved]
    print("\nApproved candidates")
    print("GEO\tN\tPlatform\tEarly-stage\tDownload/cache")
    if passing:
        for decision in passing:
            candidate = decision.candidate
            print(f"{candidate.accession}\t{candidate.sample_count or 'Unknown'}\t{candidate.platform}\t{'Yes' if decision.early_stage_documented else 'No'}\tYes")
    else:
        print("None\t-\tNo remaining candidate passed all inclusion laws\t-\tNo")
    LOGGER.info("Screened %d uncached candidates; %d approved", len(decisions), len(passing))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
