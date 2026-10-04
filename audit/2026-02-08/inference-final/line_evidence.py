"""Reusable pure static evidence helpers; no project imports or deserialization."""
import hashlib
import json
from pathlib import Path


def line_evidence(path, start, end):
    """Return a bounded, numbered UTF-8 source span (inclusive, 1-based)."""
    if start < 1 or end < start:
        raise ValueError('invalid span')
    lines = Path(path).read_text(encoding='utf-8').splitlines()
    return '\n'.join(f'{i}: {lines[i-1]}' for i in range(start, min(end, len(lines)) + 1))


def static_identity(path):
    data = Path(path).read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def extracted_spans(evidence_path, spans):
    """Number already-extracted evidence without importing audited code."""
    entries = json.loads(Path(evidence_path).read_text(encoding='utf-8'))
    index = {e['path'].replace('\\', '/'): e for e in entries}
    return {p: '\n'.join(f'{i}: {index[p]["lines"][i-1]}' for a,b in ranges for i in range(a,b+1)) for p,ranges in spans.items()}
