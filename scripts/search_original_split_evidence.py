#!/usr/bin/env python3

from __future__ import annotations

import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


ROOT = Path.home() / "ovarian-ai"

SEARCH_DIRS = [
    ROOT / "data" / "external" / "elife_28932",
    ROOT / "data" / "external" / "elife-28932-code1",
]

PATTERNS = [
    r"\btraining set\b",
    r"\btesting set\b",
    r"\btraining\b",
    r"\btesting\b",
    r"\btrain\b",
    r"\btest\b",
    r"\b135\b",
    r"\b44\b",
    r"\brandom\b",
    r"\brandomization\b",
    r"\brandomly\b",
    r"\bsplit\b",
    r"\bvalidation\b",
    r"\bvalidation set\b",
]

COMPILED = [
    re.compile(p, re.IGNORECASE)
    for p in PATTERNS
]


def clean_text(text: str) -> str:
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def matching_patterns(text: str) -> list[str]:
    found = []

    for pattern, compiled in zip(PATTERNS, COMPILED):
        if compiled.search(text):
            found.append(pattern)

    return found


def print_context(
    filename: str,
    text: str,
    label: str,
) -> None:

    cleaned = clean_text(text)

    matches = []

    for compiled in COMPILED:
        for match in compiled.finditer(cleaned):

            start = max(
                0,
                match.start() - 180,
            )

            end = min(
                len(cleaned),
                match.end() + 260,
            )

            context = cleaned[start:end]

            matches.append(context)

            if len(matches) >= 8:
                break

        if len(matches) >= 8:
            break

    if not matches:
        return

    print()
    print("=" * 100)
    print(f"{filename} [{label}]")
    print("=" * 100)

    for i, context in enumerate(matches, start=1):
        print()
        print(f"CONTEXT {i}")
        print(context)


def process_docx(path: Path) -> None:

    try:
        with zipfile.ZipFile(path) as z:
            xml = z.read("word/document.xml")
    except Exception as exc:
        print(f"[ERROR] {path.name}: {exc}")
        return

    try:
        root = ET.fromstring(xml)
    except Exception as exc:
        print(f"[ERROR] XML parse failed for {path.name}: {exc}")
        return

    text_parts = []

    for elem in root.iter():
        if elem.tag.endswith("}t"):
            text_parts.append(elem.text or "")

    text = " ".join(text_parts)

    print_context(
        path.name,
        text,
        "DOCX",
    )


def process_xml(path: Path) -> None:

    try:
        text = path.read_text(
            encoding="utf-8",
            errors="ignore",
        )
    except Exception as exc:
        print(f"[ERROR] {path.name}: {exc}")
        return

    print_context(
        path.name,
        text,
        "XML",
    )


def main() -> None:

    print("=" * 100)
    print("eLIFE ORIGINAL TRAIN/TEST SPLIT EVIDENCE SEARCH")
    print("=" * 100)

    files = []

    for directory in SEARCH_DIRS:

        if not directory.exists():
            continue

        for path in sorted(directory.rglob("*")):

            if not path.is_file():
                continue

            if path.suffix.lower() in {
                ".docx",
                ".xml",
            }:
                files.append(path)

    print()
    print(f"Files searched: {len(files)}")

    for path in files:

        if path.suffix.lower() == ".docx":
            process_docx(path)

        elif path.suffix.lower() == ".xml":
            process_xml(path)

    print()
    print("=" * 100)
    print("SEARCH COMPLETE")
    print("=" * 100)
    print()
    print("No model was trained.")
    print()
    print("Look especially for evidence containing:")
    print("  - 135 training / 44 testing")
    print("  - explicit sample IDs")
    print("  - random assignment")
    print("  - seed values")
    print("  - a list/table identifying training and testing patients")
    print()


if __name__ == "__main__":
    main()
