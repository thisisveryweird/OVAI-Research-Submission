#!/usr/bin/env python3
"""Download public ovarian-cancer datasets used by this research prototype.

Internet access is expected on the user's Linux machine. This script does not
silently overwrite existing files.
"""
from __future__ import annotations

import argparse
import gzip
import shutil
import sys
from pathlib import Path
from urllib.request import Request, urlopen

PROJECT = Path(__file__).resolve().parents[1]
RAW = PROJECT / "data" / "raw"
EXTERNAL = PROJECT / "data" / "external"

FILES = {
    "gse94533_cohort1": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94533/suppl/"
        "GSE94533_Processed_file_Cohort1.xlsx.gz",
        RAW / "GSE94533_Processed_file_Cohort1.xlsx.gz",
    ),
    "gse94533_cohort2": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94533/suppl/"
        "GSE94533_Processed_file_Cohort2.xlsx.gz",
        RAW / "GSE94533_Processed_file_Cohort2.xlsx.gz",
    ),
    "methylbert_data_s1": (
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC11384945/bin/mmc2.xlsx",
        EXTERNAL / "MethylBERT_EOC_DataS1_mmc2.xlsx",
    ),
}


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0:
        print(f"[skip] {destination} already exists")
        return

    print(f"[download] {url}")
    request = Request(url, headers={"User-Agent": "ovarian-ai-research/0.1"})
    with urlopen(request, timeout=60) as response, destination.open("wb") as out:
        shutil.copyfileobj(response, out)
    print(f"[saved] {destination} ({destination.stat().st_size:,} bytes)")


def gunzip(source: Path, destination: Path) -> None:
    if destination.exists() and destination.stat().st_size > 0:
        return
    print(f"[extract] {source.name}")
    with gzip.open(source, "rb") as src, destination.open("wb") as dst:
        shutil.copyfileobj(src, dst)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset",
        choices=[*FILES.keys(), "all"],
        default="all",
        help="Dataset to download.",
    )
    args = parser.parse_args()

    names = list(FILES) if args.dataset == "all" else [args.dataset]
    failures = []
    for name in names:
        url, path = FILES[name]
        try:
            download(url, path)
            if path.suffix == ".gz":
                gunzip(path, path.with_suffix(""))
        except Exception as exc:
            print(f"[ERROR] {name}: {exc}", file=sys.stderr)
            failures.append(name)

    if failures:
        print("\nSome downloads failed:", ", ".join(failures), file=sys.stderr)
        print("Check internet access and rerun this script.", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
