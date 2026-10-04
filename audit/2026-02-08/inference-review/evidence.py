"""Pure offline audit evidence: hashes/line-numbered UTF-8 text or binary string hits.
Never deserializes artifacts or imports project modules.
"""
import argparse
import hashlib
import json
import pathlib
import re


def evidence(path):
    p = pathlib.Path(path)
    b = p.read_bytes()
    result = {'path': str(p), 'bytes': len(b), 'sha256': hashlib.sha256(b).hexdigest()}
    if p.suffix.lower() == '.joblib':
        result['mode'] = 'static_binary_unexecuted'
        result['hits'] = [{'offset': m.start(), 'text': m.group().decode('ascii')} for m in re.finditer(rb'[ -~]{5,}', b) if any(k in m.group() for k in (b'feature', b'sklearn', b'Classifier', b'Scaler', b'version'))]
    else:
        text = b.decode('utf-8-sig')
        result['mode'] = 'full_text_read_unexecuted'
        result['lines'] = text.splitlines()
    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('paths', nargs='+')
    parser.add_argument('--output')
    args = parser.parse_args()
    rows = [evidence(p) for p in args.paths]
    serialized = json.dumps(rows, indent=2)
    if args.output:
        pathlib.Path(args.output).write_text(serialized, encoding='utf-8')
    else:
        print(serialized)
