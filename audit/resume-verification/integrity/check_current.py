"""New verification driver; never runs earlier audit generators."""
import csv
import json
import os
import re
import runpy
import subprocess
from pathlib import Path

root = Path.cwd()
out = root / 'audit/resume-verification/integrity'
h = runpy.run_path(str(out / 'integrity_helper.py'))
verify, digest = h['verify'], h['digest']
def save(name, data):
    with (out / name).open('x', encoding='utf-8') as f:
        json.dump(data, f, indent=2)

manifest = json.loads((root / 'audit/2026-02-08/review/HASH_MANIFEST.json').read_text())
review_files = sorted(p.relative_to(root).as_posix() for p in (root / 'audit/2026-02-08/review').rglob('*') if p.is_file())
listed = re.findall(r'^- `([^`]+)`', (root / 'audit/2026-02-08/review/FILES_CHANGED.md').read_text(), re.M)
audit_files = sorted(p.relative_to(root).as_posix() for p in (root / 'audit/2026-02-08').rglob('*') if p.is_file())
historical = json.loads((root / 'audit/2026-02-08/integration/PRESERVATION_CHECKS.json').read_text())
commands = [['git', '-c', 'core.fsmonitor=false', 'rev-parse', 'HEAD'], ['git', '-c', 'core.fsmonitor=false', 'status', '--short', '--untracked-files=all'], ['git', 'diff', '--no-ext-diff', '--stat'], ['git', 'diff', '--no-ext-diff', '--name-only'], ['git', 'diff', '--no-ext-diff', '--cached', '--name-only']]
git = []
for command in commands:
    p = subprocess.run(command, capture_output=True, text=True, env=dict(os.environ, GIT_OPTIONAL_LOCKS='0'), timeout=30)
    git.append(dict(command=command, code=p.returncode, stdout=p.stdout, stderr=p.stderr))
save('evidence.json', dict(manifest=verify(root, manifest), review_files=review_files, manifest_uncovered=sorted(set(review_files)-set(manifest)), file_list_count=len(listed), file_list_missing=sorted(set(listed)-set(audit_files)), file_list_unlisted=sorted(set(audit_files)-set(listed)), session_files=[p for p in audit_files if '/session/' in p], handoff_reference_search={name:[p for p in audit_files if p.endswith('/'+name)] for name in ['FILES_CREATED.csv','handoff_inventory.py']}, historical_source='audit/2026-02-08/integration/PRESERVATION_CHECKS.json', historical=verify(root,{r['path']:r['initial_sha256'] for r in historical}), git=git))
synthetic = out / 'synthetic'
synthetic.mkdir()
with (synthetic / 'sample.txt').open('x') as f:
    f.write('offline synthetic bytes')
entries = {'sample.txt':digest(synthetic/'sample.txt'), 'missing.txt':'0'*64, '../escape':'0'*64}
rows = verify(synthetic,entries)
bad = verify(synthetic,{'sample.txt':'0'*64})
assert [r['status'] for r in rows] == ['match','missing','unsafe_path']
assert bad[0]['status']=='mismatch'
save('synthetic_results.json',dict(cases=rows+bad,passed=4))
save('preservation_recheck.json',verify(root,json.loads((out/'initial_snapshot.json').read_text())))
print('Evidence and four synthetic cases saved; initial snapshot rechecked.')
