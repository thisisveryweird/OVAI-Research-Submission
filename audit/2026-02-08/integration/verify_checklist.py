"""Reusable stdlib-only verifier: inventory IDs exactly once and source ranges exist."""
import csv,re,json
from pathlib import Path

def verify(inventory,checklist,root):
    expected=re.findall(r'^- (PDF-P\d+-R\d+):',Path(inventory).read_text(encoding='utf-8-sig'),re.M)
    with Path(checklist).open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
    seen=[r['recommendation_id'] for r in rows]; errors=[]
    if sorted(expected)!=sorted(seen) or len(seen)!=len(set(seen)): errors.append('IDs missing, duplicated or unexpected')
    allowed={'implemented_verified','partial','unresolved','obsolete','unsupported','needs_external_data_owner_decision'}
    for r in rows:
        if r['disposition'] not in allowed: errors.append(r['recommendation_id']+': bad disposition')
        p=Path(root)/r['source_path']
        n=len(p.read_text(encoding='utf-8-sig').splitlines()) if p.is_file() else 0
        for span in r['source_lines'].split(';'):
            bounds=[int(x) for x in span.split('-')]
            if not 1<=bounds[0]<=bounds[-1]<=n: errors.append(r['recommendation_id']+': missing/out-of-range '+str(p)+':'+span)
        if not r['next_verifiable_step']: errors.append(r['recommendation_id']+': empty next step')
    return {'expected':len(expected),'mapped':len(rows),'errors':errors,'passed':not errors,'clinical_test':False}
if __name__=='__main__':
    root=Path.cwd();base=root/'audit/2026-02-08'
    result=verify(base/'pdf/recommendation-inventory.md',base/'integration/PDF_TO_REPOSITORY_CHECKLIST.csv',root)
    print(json.dumps(result,indent=2));raise SystemExit(0 if result['passed'] else 1)
