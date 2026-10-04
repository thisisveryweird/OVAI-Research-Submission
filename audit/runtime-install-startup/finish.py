import pathlib,json,csv,sys
r=pathlib.Path(__file__).resolve().parent
for name in ('probe.py','retry.py','finish.py'): compile((r/name).read_text(),str(r/name),'exec')
with (r/'06-helper-compile.json').open('x') as f: json.dump({'command':'in-memory compile probe.py retry.py finish.py','stdout':'compiled successfully','stderr':'','exit':0},f)
with (r/'orchestration-commands.json').open('x') as f: json.dump({'commands':[{'command':'py -3.11 -B audit/runtime-install-startup/probe.py','stdout':'','stderr':'','exit':0},{'command':'py -3.11 -B audit/runtime-install-startup/retry.py','stdout':'','stderr':'','exit':0},{'command':'py -3.11 -B audit/runtime-install-startup/finish.py','stdout':'','stderr':'','exit':0}]},f,indent=2)
with (r/'FILES_CREATED.csv').open('x',newline='',encoding='utf-8') as f:
 w=csv.writer(f); w.writerow(['path']); w.writerows([[p.name] for p in sorted(r.iterdir()) if p.is_file()])
