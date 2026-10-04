"""Pure bounded source evidence; no project imports or execution."""
from pathlib import Path
from hashlib import sha256


def evidence(path, first=1, last=40):
    if first < 1 or last < first or last - first >= 200:
        raise ValueError('Request 1..200 lines')
    raw = Path(path).read_bytes()
    if len(raw) > 1_000_000:
        raise ValueError('Source exceeds 1 MB bound')
    lines = raw.decode('utf-8-sig').splitlines()
    return {'path': str(path), 'sha256': sha256(raw).hexdigest(),
            'bytes': len(raw), 'total_lines': len(lines),
            'lines': [(i, lines[i-1]) for i in range(first, min(last, len(lines))+1)]}
