#!/usr/bin/env python3

from __future__ import annotations

import os
import re
import shutil
import urllib.request
import zipfile
from pathlib import Path

import pandas as pd


ROOT = Path.home() / "ovarian-ai"

EXTERNAL = ROOT / "data" / "external"
ELIFE_DIR = EXTERNAL / "elife_28932"

ZIP_PATH = EXTERNAL / "elife-28932-supp-v1.zip"

# The eLife article identifies this as the complete supplementary archive.
SUP_URL = (
    "https://cdn.elifesciences.org/articles/"
    "28932/elife-28932-supp-v1.zip"
)

# Original source-code archive for the miRNA-seq neural network.
CODE_URL = (
    "https://cdn.elifesciences.org/articles/"
    "28932/elife-28932-code1-v2.zip"
)

CODE_ZIP = EXTERNAL / "elife-28932-code1-v2.zip"


def download(url: str, destination: Path) -> None:

    destination.parent.mkdir(parents=True, exist_ok=True)

    if destination.exists() and destination.stat().st_size > 0:
        print(f"[OK] Already downloaded: {destination}")
        return

    print()
    print(f"[DOWNLOAD] {url}")
    print(f"           -> {destination}")

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "ovarian-ai-research-prototype/1.0"
        },
    )

    with urllib.request.urlopen(
        request,
        timeout=180,
    ) as response:

        data = response.read()

    destination.write_bytes(data)

    print(
        f"[OK] Downloaded {destination.stat().st_size:,} bytes"
    )


def extract_zip(zip_path: Path, destination: Path) -> None:

    if destination.exists():

        # Don't automatically delete anything.
        # We simply reuse an existing extracted archive.
        existing = list(destination.rglob("*"))

        if existing:
            print(
                f"[OK] Existing extraction found: {destination}"
            )
            return

    destination.mkdir(parents=True, exist_ok=True)

    print()
    print(f"[EXTRACT] {zip_path}")
    print(f"          -> {destination}")

    with zipfile.ZipFile(zip_path, "r") as z:

        z.extractall(destination)

    print("[OK] Extraction complete.")


def find_files(root: Path, extensions: set[str]) -> list[Path]:

    results = []

    for path in root.rglob("*"):

        if path.is_file() and path.suffix.lower() in extensions:
            results.append(path)

    return sorted(results)


def normalize_id(value: object) -> str:

    if pd.isna(value):
        return ""

    text = str(value).strip().lower()

    # Correct handling of Excel-style numeric IDs:
    # 50057.0 -> 50057
    if re.fullmatch(r"\d+\.0+", text):
        text = text.split(".", 1)[0]

    # E-1001 -> e1001
    text = re.sub(r"[^a-z0-9]", "", text)

    return text


def inspect_excel(path: Path) -> None:

    print()
    print("=" * 110)
    print(f"EXCEL FILE: {path.name}")
    print("=" * 110)

    try:
        xls = pd.ExcelFile(path)
    except Exception as exc:
        print(f"Could not open workbook: {exc}")
        return

    print()
    print("SHEETS:")

    for i, sheet in enumerate(xls.sheet_names, 1):
        print(f"  {i:3d}. {sheet}")

    # Search every sheet for:
    #   - diagnostic words
    #   - sample IDs
    #   - clinical fields
    #
    # We do this in a deliberately bounded way first so the terminal
    # does not get flooded with millions of expression values.

    keywords = [
        "cancer",
        "benign",
        "borderline",
        "control",
        "diagnosis",
        "histology",
        "stage",
        "sample",
        "patient",
        "case",
        "validation",
        "training",
        "testing",
        "prediction",
    ]

    target_ids = {
        "e1001",
        "e1056",
        "e1057",
        "e1058",
        "e1063",
        "50217",
        "50057",
        "70088",
    }

    for sheet in xls.sheet_names:

        print()
        print("-" * 110)
        print(f"SHEET: {sheet}")
        print("-" * 110)

        try:
            # Read only first 25 rows without assuming a header.
            preview = pd.read_excel(
                path,
                sheet_name=sheet,
                header=None,
                nrows=25,
                engine="openpyxl",
            )
        except Exception as exc:
            print(f"Could not read sheet: {exc}")
            continue

        print(f"Preview shape: {preview.shape}")

        print()
        print("FIRST 15 ROWS / FIRST 30 COLUMNS:")

        print(
            preview.iloc[:, :30].to_string(
                index=True,
                header=False,
            )
        )

        # Search the first 100 rows for obvious metadata.
        try:
            search_df = pd.read_excel(
                path,
                sheet_name=sheet,
                header=None,
                nrows=100,
                engine="openpyxl",
            )

            flattened = []

            for value in search_df.to_numpy().ravel():

                if pd.isna(value):
                    continue

                text = str(value).strip()

                if text:
                    flattened.append(text)

            joined_lower = " | ".join(
                flattened
            ).lower()

            found_keywords = [
                keyword
                for keyword in keywords
                if keyword in joined_lower
            ]

            if found_keywords:

                print()
                print("KEYWORDS FOUND:")
                for keyword in found_keywords:
                    print(f"  - {keyword}")

            # Search specifically for sample IDs.
            found_ids = set()

            for value in flattened:

                normalized = normalize_id(value)

                if normalized in target_ids:
                    found_ids.add(normalized)

            if found_ids:

                print()
                print("TARGET SAMPLE IDS FOUND:")
                for sample_id in sorted(found_ids):
                    print(f"  - {sample_id}")

        except Exception as exc:
            print(
                f"Metadata search failed for {sheet}: {exc}"
            )


def inspect_code(root: Path) -> None:

    print()
    print("=" * 110)
    print("ORIGINAL eLIFE SOURCE CODE")
    print("=" * 110)

    files = []

    for path in root.rglob("*"):

        if (
            path.is_file()
            and path.suffix.lower() in {
                ".r",
                ".txt",
                ".m",
                ".pmml",
                ".xml",
                ".csv",
                ".xlsx",
            }
        ):

            files.append(path)

    if not files:
        print("No source-code files found.")
        return

    for path in sorted(files):

        print()
        print("-" * 110)
        print(f"FILE: {path}")
        print("-" * 110)

        try:
            text = path.read_text(
                encoding="utf-8",
                errors="replace",
            )
        except Exception:
            print("Could not read as text.")
            continue

        lower = text.lower()

        interesting_terms = [
            "cancer",
            "benign",
            "borderline",
            "control",
            "training",
            "testing",
            "validation",
            "mirna",
            "e-1001",
            "e-1056",
        ]

        found = [
            term
            for term in interesting_terms
            if term in lower
        ]

        if found:
            print(
                "Interesting terms:",
                ", ".join(found),
            )

        # Print lines containing diagnostic/classification information.
        lines = text.splitlines()

        relevant_lines = []

        for line_no, line in enumerate(lines, 1):

            low = line.lower()

            if any(
                term in low
                for term in [
                    "cancer",
                    "benign",
                    "borderline",
                    "control",
                    "class",
                    "label",
                    "training",
                    "testing",
                    "validation",
                    "sample",
                    "mirna",
                ]
            ):

                clean = line.strip()

                if clean:
                    relevant_lines.append(
                        (line_no, clean[:300])
                    )

        for line_no, line in relevant_lines[:120]:

            print(
                f"{line_no:5d}: {line}"
            )


def main() -> int:

    print("=" * 110)
    print("ORIGINAL GSE94533 / eLIFE DATASET INVESTIGATION")
    print("=" * 110)

    EXTERNAL.mkdir(
        parents=True,
        exist_ok=True,
    )

    # --------------------------------------------------------------
    # Download original eLife supplementary archive
    # --------------------------------------------------------------

    try:
        download(
            SUP_URL,
            ZIP_PATH,
        )
    except Exception as exc:

        print()
        print("FAILED TO DOWNLOAD eLife supplementary archive.")
        print(exc)

        print()
        print(
            "The Linux VM needs internet access. "
            "Do not modify the existing GSE94533 files."
        )

        return 1

    # --------------------------------------------------------------
    # Extract archive
    # --------------------------------------------------------------

    try:
        extract_zip(
            ZIP_PATH,
            ELIFE_DIR,
        )
    except Exception as exc:

        print()
        print("FAILED TO EXTRACT supplementary archive.")
        print(exc)
        return 1

    # --------------------------------------------------------------
    # Find Excel supplementary files
    # --------------------------------------------------------------

    excel_files = find_files(
        ELIFE_DIR,
        {".xlsx", ".xls"},
    )

    print()
    print("=" * 110)
    print("SUPPLEMENTARY EXCEL FILES")
    print("=" * 110)

    if not excel_files:

        print("No Excel files found.")
    else:

        for path in excel_files:
            print(
                f"{path} "
                f"({path.stat().st_size:,} bytes)"
            )

    # Inspect the main supplementary dataset first.
    for path in excel_files:

        if (
            "supp6" in path.name.lower()
            or "supp" in path.name.lower()
        ):

            inspect_excel(path)

    # --------------------------------------------------------------
    # Download and inspect original neural-network source
    # --------------------------------------------------------------

    try:
        download(
            CODE_URL,
            CODE_ZIP,
        )
    except Exception as exc:

        print()
        print(
            "WARNING: Could not download original source-code archive."
        )
        print(exc)

    CODE_DIR = EXTERNAL / "elife-28932-code1"

    if CODE_ZIP.exists():

        try:
            extract_zip(
                CODE_ZIP,
                CODE_DIR,
            )
            inspect_code(CODE_DIR)
        except Exception as exc:

            print()
            print(
                "WARNING: Could not inspect source code."
            )
            print(exc)

    # --------------------------------------------------------------
    # Final instructions
    # --------------------------------------------------------------

    print()
    print("=" * 110)
    print("INVESTIGATION COMPLETE")
    print("=" * 110)

    print()
    print("Original supplementary archive:")
    print(f"  {ZIP_PATH}")

    print()
    print("Extracted archive:")
    print(f"  {ELIFE_DIR}")

    print()
    print(
        "IMPORTANT: DO NOT TRAIN THE NEW MODEL YET."
    )

    print(
        """
The goal here is to recover the exact sample-level information used
by the original publication.

We specifically want to establish:

    sample ID
        ↓
    original diagnosis
        ↓
    original training/validation assignment
        ↓
    original miRNA expression

Once those are established, we can reproduce the historical miRNA
baseline and then deliberately improve upon it.

That gives your research project a defensible progression:

    Published model
          ↓
    faithful reproduction
          ↓
    clinical/symptom baseline
          ↓
    additional genomic modalities
          ↓
    multimodal model
          ↓
    independent validation
    """
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
