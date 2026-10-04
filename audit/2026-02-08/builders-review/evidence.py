"""Pure static evidence helpers: no project imports, execution or writes."""
from hashlib import sha256
from pathlib import Path


def evidence(path, start=1, stop=None):
    """Return full-byte hash and requested 1-based source lines."""
    raw = Path(path).read_bytes()
    lines = raw.decode('utf-8').splitlines()
    if start < 1 or (stop is not None and stop < start):
        raise ValueError('Invalid line interval')
    stop = len(lines) if stop is None else min(stop, len(lines))
    return {'path': str(path), 'bytes': len(raw), 'sha256': sha256(raw).hexdigest(),
            'total_lines': len(lines),
            'lines': [(i, lines[i - 1]) for i in range(start, stop + 1)]}
