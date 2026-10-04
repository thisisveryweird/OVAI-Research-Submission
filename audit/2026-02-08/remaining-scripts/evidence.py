"""Read-only source fingerprint utility; never imports inspected source."""
from pathlib import Path
import hashlib


def fingerprint(path):
    data = Path(path).read_bytes()
    return {"path": Path(path).as_posix(), "bytes": len(data),
            "lines": len(data.decode("utf-8").splitlines()),
            "sha256": hashlib.sha256(data).hexdigest()}


def numbered_source(path):
    return "\n".join(f"{i}: {line}" for i, line in enumerate(
        Path(path).read_text(encoding="utf-8").splitlines(), 1))
