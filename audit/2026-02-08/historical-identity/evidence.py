"""Pure reusable evidence function: no I/O, imports of project code, or execution."""
import hashlib


def source_evidence(source: bytes, first: int = 1, last: int | None = None) -> dict:
    lines = source.decode('utf-8').splitlines()
    stop = len(lines) if last is None else last
    if not 1 <= first <= stop <= len(lines):
        raise ValueError('invalid inclusive line range')
    return {'sha256': hashlib.sha256(source).hexdigest(), 'bytes': len(source),
            'lines': len(lines), 'excerpt': '\n'.join(
                f'{i}: {lines[i-1]}' for i in range(first, stop + 1))}
