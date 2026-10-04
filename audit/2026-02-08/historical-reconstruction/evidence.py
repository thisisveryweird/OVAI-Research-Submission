"""Pure static evidence helper: no source import or execution, no writes."""
from hashlib import sha256

def describe_bytes(path, data, ranges):
    """Caller supplies raw bytes and actually reviewed inclusive line ranges."""
    lines = data.decode('utf-8').splitlines()
    for start, end in ranges:
        if not 1 <= start <= end <= len(lines):
            raise ValueError('Invalid reviewed range')
    return {'path': str(path), 'bytes': len(data), 'sha256': sha256(data).hexdigest(),
            'lines': len(lines), 'reviewed_ranges': ranges,
            'complete': set(range(1, len(lines)+1)).issubset(
                {i for a,b in ranges for i in range(a,b+1)})}

def excerpt(text, start, end):
    lines = text.splitlines()
    if not 1 <= start <= end <= len(lines):
        raise ValueError('Invalid excerpt range')
    return '\n'.join(f'{i}: {lines[i-1]}' for i in range(start,end+1))
