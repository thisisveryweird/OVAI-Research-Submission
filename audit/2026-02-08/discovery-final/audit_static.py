"""Pure bounded static evidence helper. Does not import/execute inspected sources.
Usage (separately authorized): python audit_static.py REPOSITORY_ROOT
Prints hashes, function spans and concise scientific-contract anchors as JSON.
No writes, recursion, dependency installation, model loading or network.
Hashing is NOT substantive read coverage; that is a human audit assertion.
"""
from pathlib import Path
import ast
import hashlib
import json
import sys

FILES = (
    'run_phase6_feature_selection.py', 'run_phase7_clinical_experiments.py',
    'run_phase8_model_experiments.py', 'run_phase9_data_augmentation.py',
    'test_high_capacity_models.py', 'run_pairwise_ratio_experiment.py',
    'run_rank_transform_experiment.py', 'analyze_gse94533_errors.py',
)
ANCHORS = ('fit_transform', 'fit_resample', 'selected =', 'SIGNATURE =',
           'scoring=', 'cv=', 'probability >=', 'probabilities >=',
           'metrics(labels, probability)', 'meta_train_columns',
           'stage_I_vs_non_cancer', 'best_augmented', 'mirna_200b_ca125',
           'mirna_296_age')

def inspect_bytes(raw):
    text = raw.decode('utf-8-sig')
    lines = text.splitlines()
    tree = ast.parse(text)
    return {
        'sha256': hashlib.sha256(raw).hexdigest(), 'bytes': len(raw),
        'lines': len(lines), 'coverage': 'static extraction only; not human read certification',
        'functions': [{'name': n.name, 'start': n.lineno, 'end': n.end_lineno}
                      for n in ast.walk(tree)
                      if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))],
        'anchors': [{'line': i, 'excerpt': line.strip()[:240]}
                    for i, line in enumerate(lines, 1)
                    if any(a in line for a in ANCHORS)],
    }

def audit(root):
    root = Path(root).resolve()
    out = {}
    for name in FILES:
        path = root / 'scripts' / name
        try:
            out['scripts/' + name] = inspect_bytes(path.read_bytes())
        except (OSError, UnicodeError, SyntaxError) as exc:
            out['scripts/' + name] = {'error': str(exc)}
    return out

if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Expected repository root')
    print(json.dumps(audit(sys.argv[1]), indent=2))
