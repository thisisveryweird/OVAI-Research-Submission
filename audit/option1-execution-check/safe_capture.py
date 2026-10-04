import hashlib
import json
import pathlib
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[2]
INPUTS = ['scripts/train_mirna_baseline.py', 'scripts/demo_ovarian_cases.py', 'scripts/predict_ovarian_risk.py', 'scripts/ovarian_diagnostic_portal.py', 'Launch_Diagnostic_Portal.bat', 'Launch_Diagnostic_Portal.command', 'tests/research_inference/test_containment.py', 'tests/research_inference/ast_helpers.py']

def capture(args):
    try:
        p = subprocess.run(args, cwd=ROOT, capture_output=True, text=True, timeout=90)
        return dict(command=args, stdout=p.stdout, stderr=p.stderr, exit_code=p.returncode)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return dict(command=args, stdout='', stderr=str(exc), exit_code=None)

def save_new(folder, name, data):
    with (folder / name).open('x', encoding='utf-8') as f:
        f.write(data)

def main():
    folder = pathlib.Path(__file__).parent
    names = ['test_output.txt', 'git_evidence.txt', 'inputs_sha256.json']
    if any((folder / n).exists() for n in names):
        raise FileExistsError('Output already exists; refusing overwrite')
    save_new(folder, names[0], json.dumps(capture(['py', '-3.11', '-B', '-m', 'unittest', 'discover', '-s', 'tests/research_inference', '-v']), indent=2))
    commands = [['git', 'rev-parse', 'HEAD'], ['git', 'status', '--short', '--untracked-files=all'], ['git', 'diff', '--name-only'], ['git', 'diff', '--check']]
    save_new(folder, names[1], json.dumps([capture(c) for c in commands], indent=2))
    save_new(folder, names[2], json.dumps({p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in INPUTS}, indent=2))

if __name__ == '__main__':
    main()
