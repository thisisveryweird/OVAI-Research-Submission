from pathlib import Path
import csv,re
B=Path('audit/2026-02-08')
rows=list(csv.DictReader((B/'review/PDF_CHECKLIST.csv').open(encoding='utf-8-sig')))
ids={'PDF-P08-R05','PDF-P08-R06','PDF-P16-R02','PDF-P17-R03','PDF-P18-R01','PDF-P20-R03','PDF-P21-R02','PDF-P06-R03'}
seen=set()
for r in rows:
 if r['status']!='implemented_verified' and r['recommendation_id'] not in ids: continue
 k=(r['source_path'],r['inclusive_lines'])
 if k in seen: continue
 seen.add(k); p=Path(k[0]); a,z=map(int,k[1].split('-')); ls=p.read_text(encoding='utf-8',errors='replace').splitlines()
 print('\nSOURCE',p,a,z)
 hits=range(a,z+1) if z-a<60 else [i for i in range(a,z+1) if re.search('stage|split|train_test|predict|threshold|filter|fit|cv|label|feature|barcode',ls[i-1],re.I)][:30]
 for i in hits: print(f'{i}: {ls[i-1]}')
for n in ['inference-final','core-final','external-final','historical-reconstruction','historical-identity','transfer-stage-review','builders-review']:
 p=B/n/'REPORT.md'; ls=p.read_text(encoding='utf-8').splitlines(); print('\nREPORT',p)
 hits=[i for i,s in enumerate(ls) if re.search('seven|nine|7 |9 |all.positive|all.negative|Stage|split|fold|align|methyl|TCGA',s,re.I)]
 for i in hits[:12]: print(f'{i+1}: {ls[i][:1000]}')
