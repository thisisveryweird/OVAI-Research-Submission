import pathlib,sys,json
root=pathlib.Path(__file__).resolve().parent
exec(compile((root/'probe.py').read_text(),str(root/'probe.py'),'exec'),{'__name__':'probe','__file__':str(root/'probe.py')})
import probe
probe.os.chdir(probe.REPO)
a=probe.command('05-apptest-retry.json',['py','-3.11','-B','-c',"import sys; sys.path.insert(0,'audit/runtime-install-startup'); import probe; import inspect; exec(compile(inspect.getsource(probe.app), '<string>', 'exec'),probe.__dict__); raise SystemExit(0 if probe.app() else 1)"])
if a.returncode==0:
    v=json.loads((root/'02-imports.json').read_text())
    probe.launch(v['stdout'].splitlines()[0])
probe.save('artifacts-after-retry.json',probe.hashes())
probe.save('artifact-comparison-retry.json',{'unchanged':json.loads((root/'artifacts-before.json').read_text())==probe.hashes()})
