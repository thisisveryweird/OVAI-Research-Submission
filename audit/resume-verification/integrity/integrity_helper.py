"""Offline read-only SHA256 verifier. Outputs JSON to stdout; no project imports."""
import argparse
import hashlib
import json
from pathlib import Path


def digest(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()


def verify(root, entries):
    root = Path(root).resolve()
    rows = []
    for name, expected in entries.items():
        p = (root / name).resolve()
        if not p.is_relative_to(root):
            rows.append(dict(path=name, expected=expected, status='unsafe_path'))
            continue
        try:
            actual = digest(p)
            status = 'match' if actual == expected else 'mismatch'
        except FileNotFoundError:
            actual, status = None, 'missing'
        except OSError as e:
            actual, status = None, 'read_error:' + str(e)
        rows.append(dict(path=name, expected=expected, actual=actual, status=status))
    return rows


if __name__ == '__main__':
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument('root')
    a.add_argument('manifest')
    args = a.parse_args()
    rows = verify(args.root, json.loads(Path(args.manifest).read_text(encoding='utf-8')))
    print(json.dumps(rows, indent=2))
    raise SystemExit(0 if all(r['status'] == 'match' for r in rows) else 1)
