"""Offline audit metadata checker. Never imports product code or loads models."""
import csv, hashlib, json, re, subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
BASE=ROOT/'audit/2026-02-08'
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def writecsv(name,data):
    with (OUT/name).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
def check():
    c=rows(OUT/'PDF_CHECKLIST.csv'); expected=re.findall(r'- (PDF-P\d+-R\d+):', (BASE/'pdf/recommendation-inventory.md').read_text(encoding='utf-8'))
    checks=[]
    for r in c:
        p=ROOT/r['source_path']; lines=p.read_text(encoding='utf-8').splitlines();a,b=map(int,r['inclusive_lines'].split('-'))
        checks.append({'id':r['recommendation_id'],'range_valid':1<=a<=b<=len(lines),'specific_evidence_nonempty':bool(r['evidence'] and r['limits'] and r['next_action']),'source_sha256':digest(p),'excerpt':'\n'.join(lines[a-1:b])})
    proofs={}
    for p in BASE.glob('*/COVERAGE.csv'):
        if p.parent.name in ('integration','review'):continue
        for r in rows(p):
            path=r.get('path','').replace('\\','/'); q=ROOT/path; h=r.get('sha256',r.get('current_sha256',''));status='; '.join(r.get(k,'') for k in ['status','read_status','coverage','read_coverage','scope','actual_review','final_substantive_read'])
            acceptable=any(s in status.lower() for s in ['complete','full source','entire source read','all extracted source lines read'])
            if q.is_file() and h and digest(q)==h and acceptable:
                proofs.setdefault(path,[]).append({'coverage':p.relative_to(ROOT).as_posix(),'status':status,'ranges':r.get('reviewed_ranges') or r.get('read_ranges') or r.get('read_coverage') or ('1-'+str(r.get('lines',r.get('extracted_source_line_count','unknown')))), 'proof_boundary':'Inherited complete-read assertion with matching current bytes; this checker does not establish substantive reading or clinical validity','sha256':h})
    scripts=sorted(p.relative_to(ROOT).as_posix() for p in (ROOT/'scripts').iterdir() if p.suffix in ('.py','.sh'))
    research=sorted(p.relative_to(ROOT).as_posix() for p in (ROOT/'research').rglob('*.py'))
    preservation=[]
    for r in json.loads((BASE/'integration/PRESERVATION_CHECKS.json').read_text(encoding='utf-8')):
        p=ROOT/r['path']; preservation.append({'path':r['path'],'initial_sha256':r['initial_sha256'],'current_sha256':digest(p) if p.is_file() else None,'unchanged':p.is_file() and digest(p)==r['initial_sha256']})
    registry=rows(OUT/'ENDPOINT_REGISTRY.csv')
    git={arg:subprocess.run(['git','--no-optional-locks']+arg.split(),cwd=ROOT,capture_output=True,text=True).stdout for arg in ['status --short','diff --stat','diff --name-only']}
    result={'recommendation_count':len(c),'ids_exact_once':sorted(expected)==sorted(r['recommendation_id'] for r in c) and len(set(r['recommendation_id'] for r in c))==87,'citations':checks,'semantic_review':'Explicit ID-based reasoned dispositions; range checks do not certify scientific truth. Owner-only recommendations cite bounded negative source evidence, not implemented clinical workflows.','registry_count':len(registry),'registry_ids_unique':len(set(r['family_id'] for r in registry))==len(registry),'script_count':len(scripts),'script_proofs':{p:proofs.get(p,[]) for p in scripts},'script_proof_gaps':[p for p in scripts if not proofs.get(p)],'research_count':len(research),'research_proofs':{p:proofs.get(p,[]) for p in research},'research_proof_gaps':[p for p in research if not proofs.get(p)],'preservation_count':len(preservation),'preservation':preservation,'all_original_hashes_unchanged':all(r['unchanged'] for r in preservation),'git_read_only':git,'own_tests':'Only this stdlib metadata checker; no clinical tests, product execution/import, fitting or model loading.','text_table_evidence':{'reports-review':'17 MD + 2 JSON;15 CSV 146684 rows; TABLE_CHECKS and DISPLAY_CHECKS','research-reports':'18 MD +2 JSON;26 CSV1889 rows; nested workbooks/cells and PMML/SOFT evidence; TABLE_CHECKS','limitation':'Inherited reports document complete text/table reads, not new rendered figures, raw clinical adjudication or model reproduction.'},'preliminary_superseded':['integration/PDF_TO_REPOSITORY_CHECKLIST.csv','integration/ENDPOINT_REGISTRY.csv'],'approval':False}
    (OUT/'QUALITY_REVIEW.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    artifacts=[]
    for p in sorted(BASE.rglob('*')):
        if p.is_file() and p.name not in ('FILES_CHANGED.md','HASH_MANIFEST.json'):
            purpose='Generated helper cache; not product source' if '__pycache__' in p.parts or p.suffix=='.pyc' else ('Authoritative review document/helper' if p.parent==OUT else 'Prior audit evidence; preliminary unless finalized subsystem report; retained unchanged by this review')
            artifacts.append(f'- `{p.relative_to(ROOT).as_posix()}` — {purpose}.')
    artifacts.extend(['- `audit/2026-02-08/review/FILES_CHANGED.md` — Exact audit artifact path inventory.','- `audit/2026-02-08/review/HASH_MANIFEST.json` — Final review hashes, excluding itself.'])
    (OUT/'FILES_CHANGED.md').write_text('# Audit artifacts, not product changes\n\nThis review writes only review/. Prior audit artifacts are listed for full accounting, not claimed as this reviewer edits. Source PDF `Claude Methodology and Errors for Ovarian Cancer Project.pdf` is OWNER-added, not assistant-modified. No original source edits.\n\n'+'\n'.join(artifacts)+'\n',encoding='utf-8')
    (OUT/'HASH_MANIFEST.json').write_text(json.dumps({p.relative_to(ROOT).as_posix():digest(p) for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='HASH_MANIFEST.json'},indent=2),encoding='utf-8')
    assert result['ids_exact_once'] and all(x['range_valid'] and x['specific_evidence_nonempty'] for x in checks)
    assert result['registry_ids_unique'] and result['all_original_hashes_unchanged']
    assert len(scripts)==101 and len(research)==8 and not result['script_proof_gaps'] and not result['research_proof_gaps']
    print(json.dumps({k:result[k] for k in ['recommendation_count','ids_exact_once','registry_count','script_count','script_proof_gaps','research_count','research_proof_gaps','preservation_count','all_original_hashes_unchanged']}))
if __name__=='__main__':check()
