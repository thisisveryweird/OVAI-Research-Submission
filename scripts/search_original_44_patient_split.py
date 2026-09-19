#!/usr/bin/env python3

from __future__ import annotations

import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd


ROOT = Path.home() / "ovarian-ai"

EXTERNAL = ROOT / "data" / "external" / "elife_28932"

TERMS = [
    "training set",
    "testing set",
    "training",
    "testing",
    "test set",
    "train set",
    "stratified random",
    "random sampling",
    "randomization",
    "135",
    "44",
    "validation set",
    "validation",
    "patient 56",
    "sample 56",
]

FILES = sorted(
    p for p in EXTERNAL.rglob("*")
    if p.is_file()
)


# ================================================================
# HELPERS
# ================================================================

def normalize_text(text: str) -> str:

    text = text.replace("\xa0", " ")
    text = re.sub(r"\s+", " ", text)

    return text.strip()


def show_context(text: str, term: str, radius: int = 350):

    lower = text.lower()

    start = 0

    while True:

        pos = lower.find(term.lower(), start)

        if pos == -1:
            break

        left = max(0, pos - radius)
        right = min(len(text), pos + len(term) + radius)

        print()
        print("-" * 100)
        print(f"TERM: {term}")
        print("-" * 100)
        print(text[left:right])

        start = pos + len(term)


# ================================================================
# TEXT / XML / CODE FILES
# ================================================================

def search_text_file(path: Path):

    try:
        text = path.read_text(errors="ignore")
    except Exception:
        return

    normalized = normalize_text(text)

    for term in TERMS:

        if term.lower() in normalized.lower():

            print()
            print("=" * 100)
            print(path)
            print("=" * 100)

            show_context(
                normalized,
                term,
            )


# ================================================================
# DOCX
# ================================================================

def extract_docx_text(path: Path) -> str:

    try:
        with zipfile.ZipFile(path) as z:
            xml = z.read("word/document.xml")
    except Exception:
        return ""

    try:
        root = ET.fromstring(xml)
    except Exception:
        return ""

    parts = []

    for node in root.iter():

        if node.tag.endswith("}t"):
            parts.append(node.text or "")

    return normalize_text(
        " ".join(parts)
    )


# ================================================================
# XLSX
# ================================================================

def inspect_xlsx(path: Path):

    try:
        excel = pd.ExcelFile(
            path,
            engine="openpyxl",
        )
    except Exception as exc:

        print()
        print(f"Could not open {path}: {exc}")
        return

    print()
    print("=" * 100)
    print(f"XLSX: {path}")
    print("=" * 100)

    print()
    print("Sheets:")

    for sheet in excel.sheet_names:
        print(f"  {sheet}")

    for sheet in excel.sheet_names:

        try:

            raw = pd.read_excel(
                path,
                sheet_name=sheet,
                header=None,
                engine="openpyxl",
            )

        except Exception:
            continue

        # Search every cell as text.
        for row_idx in range(len(raw)):

            values = [
                str(x)
                for x in raw.iloc[row_idx].tolist()
            ]

            row_text = normalize_text(
                " | ".join(values)
            )

            lower = row_text.lower()

            if any(
                term.lower() in lower
                for term in TERMS
            ):

                print()
                print(
                    f"[{sheet}] row {row_idx + 1}"
                )
                print(row_text)


# ================================================================
# DIRECTORY SEARCH
# ================================================================

print("=" * 100)
print("SEARCHING ORIGINAL eLIFE MATERIALS FOR THE 135/44 SPLIT")
print("=" * 100)

print()
print("External directory:")
print(EXTERNAL)

print()
print("Files found:", len(FILES))

for path in FILES:

    suffix = path.suffix.lower()

    if suffix in {
        ".txt",
        ".xml",
        ".csv",
        ".r",
        ".py",
        ".sas",
        ".log",
        ".md",
    }:

        search_text_file(path)

    elif suffix == ".docx":

        text = extract_docx_text(path)

        for term in TERMS:

            if term.lower() in text.lower():

                print()
                print("=" * 100)
                print(path)
                print("=" * 100)

                show_context(
                    text,
                    term,
                )

    elif suffix == ".xlsx":

        inspect_xlsx(path)


# ================================================================
# SEARCH CODE SOURCE TREE TOO
# ================================================================

print()
print("=" * 100)
print("SEARCHING ENTIRE PROJECT FOR SPLIT INFORMATION")
print("=" * 100)

for path in ROOT.rglob("*"):

    if not path.is_file():
        continue

    if path.suffix.lower() not in {
        ".py",
        ".r",
        ".xml",
        ".txt",
        ".csv",
        ".md",
        ".json",
    }:
        continue

    if "reports" in path.parts:
        continue

    try:
        text = path.read_text(errors="ignore")
    except Exception:
        continue

    lower = text.lower()

    hits = [
        term
        for term in TERMS
        if term.lower() in lower
    ]

    if hits:

        print()
        print("-" * 100)
        print(path)
        print("Hits:", ", ".join(hits))

        # Show only relevant lines.
        for line_no, line in enumerate(
            text.splitlines(),
            start=1,
        ):

            if any(
                term.lower() in line.lower()
                for term in TERMS
            ):

                print(
                    f"{line_no}: {line[:500]}"
                )


print()
print("=" * 100)
print("SEARCH COMPLETE")
print("=" * 100)

print()
print("Important:")
print(
    "Finding '135' and '44' alone is not enough. "
    "We need evidence identifying WHICH 44 patients formed the "
    "original testing set."
)
