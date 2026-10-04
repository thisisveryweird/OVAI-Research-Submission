"""Pure byte/text evidence helpers; never import inspected code or load models."""
from hashlib import sha256

def source_evidence(content: bytes, start=1, end=None):
    lines = content.decode('utf-8').splitlines()
    end = len(lines) if end is None else end
    if not 1 <= start <= end <= len(lines):
        raise ValueError('invalid inclusive line range')
    return {'sha256': sha256(content).hexdigest(), 'bytes': len(content),
            'lines': len(lines), 'evidence': [{'line': i, 'text': lines[i-1]}
                                             for i in range(start, end+1)]}

def count_labels(rows, field, allowed):
    """Count exact supplied labels; no inferred disease categories."""
    counts = {label: 0 for label in allowed}
    excluded = 0
    for row in rows:
        label = row.get(field)
        if label in counts:
            counts[label] += 1
        else:
            excluded += 1
    return {'counts': counts, 'excluded': excluded}
