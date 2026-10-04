"""Pure offline source evidence: no imports or execution of inspected sources."""
from hashlib import sha256

def source_evidence(content: bytes, needles=()):
    text = content.decode('utf-8')
    lines = text.splitlines()
    return {'sha256': sha256(content).hexdigest(), 'bytes': len(content),
            'lines': len(lines), 'matches': [
                {'line': i, 'text': line, 'needle': needle}
                for i, line in enumerate(lines, 1) for needle in needles if needle in line]}
