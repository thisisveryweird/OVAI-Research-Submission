#!/usr/bin/env python3

from __future__ import annotations

import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


ROOT = Path.home() / "ovarian-ai"
EXT = ROOT / "data" / "external" / "elife_28932"

DOCX_FILES = [
    EXT / "elife_poa_e28932_Supplementary_File_1.docx",
    EXT / "elife_poa_e28932_Supplementary_File_2.docx",
    EXT / "elife_poa_e28932_Supplementary_File_3.docx",
    EXT / "elife_poa_e28932_Supplementary_File_4.docx",
    EXT / "elife_poa_e28932_Supplementary_File_5.docx",
    EXT / "elife_poa_e28932_Transparent_Reporting_File.docx",
]

XML_FILES = list(
    EXT.glob("*.xml")
)

KEYWORDS = [
    "training set",
    "testing set",
    "test set",
    "train set",
    "testing",
    "training",
    "stratified",
    "random",
    "randomly",
    "randomization",
    "135",
    "44",
]


def extract_docx_text(path: Path) -> str:

    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml")

    root = ET.fromstring(xml)

    parts = []

    for node in root.iter():
        if node.tag.endswith("}t"):
            parts.append(node.text or "")

    text = " ".join(parts)

    return re.sub(r"\s+", " ", text)


def print_context(text: str, keyword: str, radius: int = 250):

    lower = text.lower()
    start = 0
    seen = 0

    while True:

        pos = lower.find(keyword.lower(), start)

        if pos == -1:
            break

        left = max(0, pos - radius)
        right = min(
            len(text),
            pos + len(keyword) + radius,
        )

        context = text[left:right]

        # Ignore isolated numeric hits.
        if keyword in {"135", "44"}:

            nearby = context.lower()

            phrases = [
                "training set",
                "testing set",
                "test set",
                "train set",
                "testing",
                "training",
                "sample",
                "patient",
            ]

            if not any(p in nearby for p in phrases):
                start = pos + len(keyword)
                continue

        print(f"\n  [{keyword}] ...{context}...")

        seen += 1

        # Avoid printing dozens of duplicate contexts.
        if seen >= 5:
            break

        start = pos + len(keyword)


print("=" * 90)
print("CLEAN SEARCH FOR ORIGINAL 135/44 TRAIN-TEST SPLIT EVIDENCE")
print("=" * 90)

print("\nSearching DOCX supplementary materials...\n")

for path in DOCX_FILES:

    if not path.exists():
        continue

    text = extract_docx_text(path)

    hits = []

    for keyword in KEYWORDS:

        if keyword.lower() in text.lower():
            hits.append(keyword)

    if not hits:
        continue

    print("=" * 90)
    print(path.name)
    print("=" * 90)

    # Prioritize meaningful textual phrases.
    priority = [
        "training set",
        "testing set",
        "test set",
        "train set",
        "stratified",
        "random",
        "randomly",
        "randomization",
    ]

    printed = set()

    for keyword in priority:

        if keyword.lower() in text.lower():

            print_context(
                text,
                keyword,
            )

            printed.add(keyword)


print("\n")
print("=" * 90)
print("SEARCHING XML SOURCE FILES")
print("=" * 90)

for path in XML_FILES:

    try:
        text = path.read_text(errors="ignore")
    except Exception:
        continue

    # Only report actual textual split terminology.
    pattern = re.compile(
        r".{0,250}"
        r"(training set|testing set|test set|train set|"
        r"stratif[a-z]*|random[a-z]*)"
        r".{0,250}",
        re.IGNORECASE,
    )

    matches = pattern.findall(text)

    if not matches:
        continue

    print(f"\n{path.name}")

    # Re-run with full context.
    count = 0

    for match in pattern.finditer(text):

        snippet = re.sub(
            r"\s+",
            " ",
            match.group(0),
        )

        print(f"  ...{snippet}...")

        count += 1

        if count >= 10:
            break


print("\n")
print("=" * 90)
print("SEARCH COMPLETE")
print("=" * 90)

print(
    "\nWe specifically need evidence identifying the actual 44 patients."
)
print(
    "Finding only 'N=44' or 'testing set' is not enough to reconstruct "
    "the original split."
)
