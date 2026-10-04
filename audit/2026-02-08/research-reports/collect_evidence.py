"""Read-only local research evidence collector. Writes only beside this script."""
import csv,gzip,hashlib,json,re,statistics,struct,zipfile
import xml.etree.ElementTree as E
from pathlib import Path
from collections import Counter,defaultdict
from table_evidence import inspect_table,auc
OUT=Path(__file__).parent
checks={}; coverage=[]; tables={}; claims=[]
def cov(p,kind,status,limit=''):
 b=p.read_bytes();coverage.append([p.as_posix(),kind,len(b),hashlib.sha256(b).hexdigest(),status,limit])
for p in sorted(Path('research').rglob('*')):
 if not p.is_file():continue
 if p.suffix in ('.md','.json'):
  text=p.read_text(encoding='utf-8');cov(p,'text','complete bounded text processing and substantive review')
  if p.suffix=='.json':checks[p.as_posix()]={'parsed':json.loads(text)}
  for n,line in enumerate(text.splitlines(),1):
   if line.strip():claims.append({'path':p.as_posix(),'line':n,'text':line,'classification':'source report claim; not independently clinical-adjudicated'})
 elif p.suffix=='.csv':
  check,rows=inspect_table(p);checks[p.as_posix()]=check;tables[p.name]=rows;cov(p,'CSV','all rows/cells parsed; schema missingness duplicate numeric/category checks')
 elif p.suffix=='.py':cov(p,'source','comparison only: previously finalized research-sources audit','not reaudited here')
 elif p.suffix=='.png':
  b=p.read_bytes();checks[p.as_posix()]={'signature_valid':b[:8]==b'\x89PNG\r\n\x1a\n','width_height':struct.unpack('>II',b[16:24]),'visual_status':'UNVERIFIED; metadata only; not clinical evidence'};cov(p,'PNG','metadata only','No rendering or visual verification')
 else:cov(p,'raw','full container/hash; XLSX cells or SOFT metadata parsed','no formula execution; no visual inspection; clinical adjudication incomplete')
# Existing external evidence is reused as comparison, not treated as fresh primary proof.
for name in ('binary_inventory.json','soft_metadata.json'):
 p=Path('audit/2026-02-08/external')/name
 obj=json.loads(p.read_text(encoding='utf-8'))
 checks['reused/'+name]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'top_level_type':type(obj).__name__,'entries':len(obj),'status':'existing evidence reused; originals inspected below'}
raw={};sheets={}
NS='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
for p in Path('research').rglob('*.xlsx'):
 with zipfile.ZipFile(p) as z:
  ss=[]
  if 'xl/sharedStrings.xml' in z.namelist():
   ss=[''.join(s.itertext()) for s in E.fromstring(z.read('xl/sharedStrings.xml')).findall(NS+'si')]
  wb=E.fromstring(z.read('xl/workbook.xml'));rels=E.fromstring(z.read('xl/_rels/workbook.xml.rels'));targets={r.attrib['Id']:r.attrib['Target'] for r in rels}
  info={'magic':p.read_bytes()[:4].hex(),'zip_bad_member':z.testzip(),'sheets':[]}
  for sh in wb.find(NS+'sheets'):
   target=targets[sh.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']];target=target.lstrip('/') if target.startswith('/') else 'xl/'+target
   root=E.fromstring(z.read(target)); rows=[];types=Counter();formulas=0
   for row in root.findall('.//'+NS+'sheetData/'+NS+'row'):
    vals={}
    for c in row.findall(NS+'c'):
     types[c.attrib.get('t','n')]+=1
     if c.find(NS+'f') is not None:formulas+=1
     v=c.find(NS+'v');value=v.text if v is not None else ''
     if c.attrib.get('t')=='s' and value:value=ss[int(value)]
     elif c.attrib.get('t')=='inlineStr':value=''.join(c.find(NS+'is').itertext())
     vals[re.sub(r'\d','',c.attrib['r'])]=value
    rows.append(vals)
   cols=list(rows[0]);header=[rows[0][c] for c in cols];data=[[r.get(c,'') for c in cols] for r in rows[1:]]
   sheets[(p.name,sh.attrib['name'])]=(header,data)
   counts={h:dict(Counter(r[i] for r in data)) for i,h in enumerate(header) if len(set(r[i] for r in data))<=35}
   info['sheets'].append({'name':sh.attrib['name'],'dimension':root.find(NS+'dimension').attrib if root.find(NS+'dimension') is not None else None,'data_rows':len(data),'columns':len(header),'headers':header,'cell_types':dict(types),'formula_cells':formulas,'blank_counts_by_column':{str(i+1)+':'+h:sum(r[i]=='' for r in data) for i,h in enumerate(header)},'unique_first_column':len(set(r[0] for r in data)),'exact_duplicate_rows':len(data)-len(set(map(tuple,data))),'categorical_counts':counts})
  raw[p.as_posix()]=info
for p in Path('research').rglob('*.zip'):
 with zipfile.ZipFile(p) as z:
  members=[]
  for name in z.namelist():
   b=z.read(name);root=E.fromstring(b)
   elems=list(root.iter());members.append({'member':name,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'root':root.tag,'networks':[x.attrib for x in elems if x.tag.endswith('NeuralNetwork')],'layers':[x.attrib for x in elems if x.tag.endswith('NeuralLayer')],'inputs':[x.attrib for x in elems if x.tag.endswith('NeuralInput')],'outputs':[E.tostring(x,encoding='unicode') for x in elems if x.tag.endswith('NeuralOutput')],'connections':sum(x.tag.endswith('Con') for x in elems),'probability_output_fields':[x.attrib for x in elems if x.tag.endswith('OutputField') and x.attrib.get('feature')=='probability']})
  raw[p.as_posix()]={'zip_bad_member':z.testzip(),'members':members}
for p in Path('research').rglob('*.gz'):
 samples=[];series=[];current=None;line_count=0
 with gzip.open(p,'rt',encoding='utf-8',errors='replace') as f:
  for line in f:
   line_count+=1;line=line.rstrip('\n')
   if line.startswith('^SAMPLE = '):current={'GSM':line.split(' = ',1)[1]};samples.append(current)
   elif line.startswith('!Sample_') and current is not None:
    k,v=line.split(' = ',1);current.setdefault(k,[]).append(v)
   elif line.startswith('!Series_'):series.append(line)
 labels=Counter();stages=Counter();titles=[]
 for s in samples:
  titles+=s.get('!Sample_title',[])
  chars=s.get('!Sample_characteristics_ch1',[])
  d=next((x.split(':',1)[1].strip() for x in chars if x.lower().startswith('disease state:')),'MISSING');labels[d]+=1
  if d=='ovarian cancer':stages[next((x.split(':',1)[1].strip() for x in chars if x.lower().startswith('stage:')),'MISSING')]+=1
 raw[p.as_posix()]={'decompressed_lines_read':line_count,'series_fields':series,'GSM_records':len(samples),'unique_GSM':len(set(s['GSM'] for s in samples)),'unique_titles':len(set(titles)),'disease_labels':dict(labels),'ovarian_stages':dict(stages),'sample_metadata_field_counts':dict(Counter(k for s in samples for k in s))}
# Reconcile workbook IDs/splits, histology conflict, assay ranges and saved outputs without evaluating networks.
recon={};h,d=sheets[('NIHMS2021079-supplement-1.xlsx','miRNA')]; hp,dp=sheets[('NIHMS2021079-supplement-1.xlsx','protein + metadata')];hh,dh=sheets[('NIHMS2021079-supplement-1.xlsx','histologies')]
stage=defaultdict(set)
for r in dh:stage[r[0]].add(r[2])
recon['2024']={'unique_ids':len(set(r[0] for r in d)),'ordered_ID_cohort_diagnosis_match':all(a[:3]==b[:3] for a,b in zip(d,dp)) and len(d)==len(dp),'split_diagnosis':dict(Counter(r[1]+'|'+r[2] for r in d)),'histology_unique':len(stage),'histology_exact_duplicate_rows':len(dh)-len(set(map(tuple,dh))),'conflicts':{k:sorted(v) for k,v in stage.items() if len(v)>1},'unambiguous_stages':dict(Counter(next(iter(v)) for v in stage.values() if len(v)==1)),'cancers_without_histology':sum(r[2]=='1' and r[0] not in stage for r in d)}
panel=tables['2024_adnexal_mass_mirna_panel_inventory.csv'];ranges=[]
for row in panel:
 i=h.index(row['workbook_column']);v=[float(r[i]) for r in d];ranges.append({'feature':row['2024_identifier'],'min_equal':abs(min(v)-float(row['minimum']))<1e-10,'median_equal':abs(statistics.median(v)-float(row['median']))<1e-10,'max_equal':abs(max(v)-float(row['maximum']))<1e-10,'missing':sum(r[i]=='' for r in d)})
recon['panel_ranges']=ranges
pred=tables['2022_gockley_signature_best_10_predictions.csv'];metricchecks=[]
for fn in ('2022_gockley_signature_best_10_metrics.csv','2022_gockley_published_reproduction_metrics.csv'):
 for m in tables[fn]:
  ep=m['endpoint'];split=m.get('subset',m.get('split'));rr=[r for r in pred if (split=='All' or r['set']==split) and (not ('stage_I_II' in ep or ep=='strict_stage_I_II') or r['diagnosis']!='Cancer Early' or r['stage_group']=='Stage I/II') and (not ep.endswith('benign_mass') or r['diagnosis']!='Control')]
  conf=Counter(('tp' if r['predicted_class']=='malignant' else 'fn') if r['diagnosis']=='Cancer Early' else ('fp' if r['predicted_class']=='malignant' else 'tn') for r in rr)
  scores={'margin':auc(rr,'score_margin_malignant'),'malignant':auc(rr,'output_malignant')}
  metricchecks.append({'file':fn,'endpoint':ep,'split':split,'n':len(rr),'confusion':dict(conf),'auc':scores,'reported':m,'margin_auc_matches':abs(scores['margin']-float(m.get('roc_auc_margin',m.get('roc_auc_output_margin'))))<1e-12,'malignant_auc_matches':abs(scores['malignant']-float(m.get('roc_auc_malignant_output',scores['malignant'])))<1e-12})
recon['saved_prediction_metrics']=metricchecks
recon['prediction_margin_class_checks']={'margin_mismatches':sum(abs(float(r['output_malignant'])-float(r['output_benign_control'])-float(r['score_margin_malignant']))>1e-12 for r in pred),'class_mismatches':sum((float(r['score_margin_malignant'])>0)!=(r['predicted_class']=='malignant') for r in pred),'duplicate_sample_ids':len(pred)-len(set(r['sample_id'] for r in pred))}
recon['transportability_arithmetic']=[{'row':r,'balanced_accuracy_matches':abs((float(r['sensitivity'])+float(r['specificity']))/2-float(r['balanced_accuracy']))<1e-12} for r in tables['transportability_summary.csv']]
# Parse every nested workbook inventory cell and every JSON preview; do not mistake preview for full raw data.
recon['nested_inventory']=[{'sheet':r['sheet'],'fields':{k:len(json.loads(v)) if isinstance(json.loads(v),(list,dict)) else None for k,v in r.items() if k.endswith('_json')}} for r in tables['2024_adnexal_mass_workbook_inventory.csv']]
checks['raw_primary_inspection']=raw;checks['reconciliation']=recon
(OUT/'TABLE_CHECKS.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
(OUT/'REPORT_LINE_EVIDENCE.json').write_text(json.dumps(claims,indent=2),encoding='utf-8')
with (OUT/'COVERAGE.csv').open('w',encoding='utf-8',newline='') as f:
 w=csv.writer(f);w.writerow(['path','type','bytes','sha256','read_status','limitations']);w.writerows(coverage)
print('text',sum(r[1]=='text' for r in coverage),'tables',len(tables),'rows',sum(len(t) for t in tables.values()),'raw',len(raw),'metric_checks',len(metricchecks));print(json.dumps(recon['2024']));print('range failures',sum(not all(r[k] for k in ('min_equal','median_equal','max_equal')) for r in ranges));print('metric failures',sum(not m['margin_auc_matches'] or not m['malignant_auc_matches'] for m in metricchecks));print('GSE',[(k,v.get('GSM_records'),v.get('ovarian_stages')) for k,v in raw.items() if k.endswith('gz')])
