import os, sys, subprocess, json, pathlib, hashlib, time, socket, urllib.request
ROOT=pathlib.Path(__file__).resolve().parent
REPO=ROOT.parent.parent
os.environ.update(PYTHONDONTWRITEBYTECODE='1',PIP_DISABLE_PIP_VERSION_CHECK='1',STREAMLIT_BROWSER_GATHER_USAGE_STATS='false')
def save(name,data):
    with (ROOT/name).open('x',encoding='utf-8') as f: json.dump(data,f,indent=2,default=str)
def command(name,args):
    p=subprocess.run(args,cwd=REPO,capture_output=True,text=True,errors='replace')
    save(name,dict(command=args,stdout=p.stdout,stderr=p.stderr,exit=p.returncode))
    return p
def hashes():
    return {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((REPO/'production_artifacts').iterdir()) if p.is_file()}
def app():
    from streamlit.testing.v1 import AppTest
    at=AppTest.from_file('scripts/ovarian_diagnostic_portal.py').run(timeout=30)
    errors=[x.value for x in at.error]; exceptions=[x.value for x in at.exception]
    markdown=[x.value for x in at.markdown]
    buttons=[{'label':x.label,'disabled':x.disabled,'value':x.value} for x in at.button]
    successes=[x.value for x in at.success]
    result=dict(exceptions=exceptions,errors=errors,buttons=buttons,successes=successes,research_title=any('Ovarian Research Workbench' in x and 'Research-only' in x for x in markdown),schema_loaded=any('Frozen artifact set loaded' in x for x in successes))
    result['passed']=not errors and not exceptions and result['research_title'] and result['schema_loaded'] and any(x['label']=='Execute Research Model' and x['disabled'] and not x['value'] for x in buttons)
    save('apptest-result.json',result)
    print(json.dumps(result)); return result['passed']
def launch(exe):
    port=None
    for candidate in (8501,8502):
        with socket.socket() as s:
            try: s.bind(('127.0.0.1',candidate)); port=candidate; break
            except OSError: pass
    if port is None: save('launch.json',{'blocked':'Both ports occupied; no existing process changed'}); return
    args=[exe,'-B','-m','streamlit','run','scripts/ovarian_diagnostic_portal.py','--server.address=127.0.0.1',f'--server.port={port}','--server.headless=true','--server.fileWatcherType=none','--browser.gatherUsageStats=false']
    out=(ROOT/'server.stdout.log').open('x'); err=(ROOT/'server.stderr.log').open('x')
    flags=subprocess.DETACHED_PROCESS|subprocess.CREATE_NEW_PROCESS_GROUP
    p=subprocess.Popen(args,cwd=REPO,stdin=subprocess.DEVNULL,stdout=out,stderr=err,creationflags=flags,close_fds=True)
    out.close(); err.close()
    url=f'http://127.0.0.1:{port}'; checks=[]; healthy=False
    for i in range(20):
        try:
            current=[]
            for path in ('/_stcore/health','/'):
                with urllib.request.urlopen(url+path,timeout=2) as r: current.append(dict(path=path,status=r.status,body=r.read().decode('utf-8')[:1000]))
            checks.append(current); healthy=all(x['status']==200 for x in current)
            if healthy: break
        except Exception as e: checks.append(str(e)); time.sleep(.5)
    save('launch.json',dict(command=args,interpreter=exe,pid=p.pid,url=url,healthy=healthy,checks=checks,process_exit=p.poll()))
if __name__=='__main__':
    os.chdir(REPO)
    if len(sys.argv)>1 and sys.argv[1]=='app': sys.exit(0 if app() else 1)
    save('artifacts-before.json',hashes())
    p=command('01-install.json',['py','-3.11','-m','pip','install','streamlit','openpyxl'])
    v=command('02-imports.json',['py','-3.11','-B','-c','import sys,streamlit,openpyxl,numpy,pandas,sklearn,joblib; import importlib.metadata as m; print(sys.executable); print(sys.version); print({n:m.version(n) for n in ("streamlit","openpyxl","numpy","pandas","scikit-learn","joblib")})'])
    command('03-tests.json',['py','-3.11','-B','-m','unittest','discover','-s','tests/research_inference','-v'])
    a=command('04-apptest.json',['py','-3.11','-B',str(ROOT/'probe.py'),'app'])
    if a.returncode==0 and v.returncode==0: launch(v.stdout.splitlines()[0])
    save('artifacts-after.json',hashes())
    save('artifact-comparison.json',{'unchanged':json.loads((ROOT/'artifacts-before.json').read_text())==hashes()})
