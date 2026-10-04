"""Pure standard-library read-only structural audit; no project imports or execution."""
import json
from pathlib import Path
import sys


def check(root):
    root = Path(root).resolve()
    inventory = json.loads((root / 'audit/2026-02-08/governance/inventory.json').read_text(encoding='utf-8-sig'))
    current = [p.relative_to(root).as_posix() for p in root.rglob('*')
               if p.is_file() and not any(x in p.relative_to(root).parts for x in ('.git', 'node_modules', '.venv'))]
    product = [p for p in current if not p.startswith('audit/')]
    keys = ('license', 'copying', 'data_sources', 'results_index', 'lock', 'pytest', 'tox', 'test_')
    print(json.dumps({'inventory_entries': len(inventory), 'current_product_files': len(product),
                     'governance_test_lock_candidates': [p for p in product if any(k in Path(p).name.lower() for k in keys)],
                     'ci_paths': [p for p in product if p.startswith(('.github/', '.gitlab', '.circleci/'))],
                     'tests_directory': (root / 'tests').exists()}, indent=2))
    for p in product:
        if p.endswith(('.py', '.sh', '.bat', '.command', '.md', '.txt', '.yaml', '.toml')):
            for n, line in enumerate((root / p).read_text(encoding='utf-8-sig', errors='replace').splitlines(), 1):
                if any(s in line for s in ('Path.home()', '~/ovarian-ai', '/ovarian-ai', 'pip install', 'pytest')):
                    print(f'{p}:{n}:{line.strip()}')


if __name__ == '__main__':
    check(sys.argv[1] if len(sys.argv) > 1 else '.')
