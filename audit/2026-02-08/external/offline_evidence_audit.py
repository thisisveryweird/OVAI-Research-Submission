"""Read-only offline evidence inventory and saved binary-prediction checks.
No acquisition imports, network, model evaluation/fitting, or source writes.
Outputs must be NEW files beneath this script's directory. Stdlib only.
Usage: python audit/2026-02-08/external/offline_evidence_audit.py ROOT
"""
import collections,csv,gzip,hashlib,io,json,pathlib,sys,zipfile,xml.etree.ElementTree as ET
OUT=pathlib.Path(__file__).resolve().parent
ROOT=pathlib.Path(sys.argv[1]).resolve() if len(sys.argv)>1 else OUT.parents[2]
def save(name,text):
    with (OUT/name).open('x',encoding='utf-8',newline='') as f:f.write(text)
def xlsx(blob):
    ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            strings=[''.join(n.itertext()) for n in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si',ns)]
        sheets=[]
        for name in z.namelist():
            if name.startswith('xl/worksheets/sheet') and name.endswith('.xml'):
                rows=[]
                for row in ET.fromstring(z.read(name)).findall('.//m:row',ns):
                    cells={}
                    for c in row.findall('m:c',ns):
                        v=c.find('m:v',ns); value=v.text if v is not None else ''
                        if c.get('t')=='s' and value:value=strings[int(value)]
                        if c.get('t')=='inlineStr':value=''.join(c.find('m:is',ns).itertext())
                        cells[c.get('r')]=value
                    rows.append(cells)
                sheets.append({'sheet':name,'rows':rows})
        return {'workbook_xml':z.read('xl/workbook.xml').decode(),'sheets':sheets}
def prediction_check(path,text):
    rows=list(csv.DictReader(io.StringIO(text))); cols=list(rows[0]) if rows else []
    label=next((c for c in ['y_true','label','target','true_label','actual','y'] if c in cols),None)
    pred=next((c for c in ['y_pred','prediction','predicted_label','pred_label','pred'] if c in cols),None)
    result={'path':str(path.relative_to(ROOT)),'rows':len(rows),'columns':cols,'label_column':label,'prediction_column':pred}
    if label and pred:
        counts=collections.Counter(); invalid=0
        for r in rows:
            try:a,b=float(r[label]),float(r[pred])
            except (ValueError,TypeError):invalid+=1;continue
            if a not in (0,1) or b not in (0,1):invalid+=1;continue
            counts[(int(a),int(b))]+=1
        tn,fp,fn,tp=[counts[k] for k in [(0,0),(0,1),(1,0),(1,1)]]
        result.update(tn=tn,fp=fp,fn=fn,tp=tp,invalid=invalid,sensitivity=tp/(tp+fn) if tp+fn else None,specificity=tn/(tn+fp) if tn+fp else None)
    else:result['status']='No unambiguous stored binary truth/prediction columns; no threshold guessed.'
    return result
files=sorted(p for d in ['research','reports'] for p in (ROOT/d).rglob('*') if p.is_file())
files+=sorted(p for p in (ROOT/'scripts').glob('*.py') if any(t in p.name for t in ['external','gse106817','pmml','original','historical','split']))
coverage=[]; texts=[]; binary=[]; predictions=[]; soft=[]
for p in files:
    blob=p.read_bytes(); item={'path':str(p.relative_to(ROOT)),'bytes':len(blob),'sha256':hashlib.sha256(blob).hexdigest()}
    try:
        if p.suffix=='.xlsx':
            binary.append({'path':item['path'],'content':xlsx(blob)});item['status']='all worksheet cells decoded (formula caches only; no formula execution)'
        elif p.suffix=='.zip':
            entries=[]
            with zipfile.ZipFile(io.BytesIO(blob)) as z:
                for n in z.namelist():
                    b=z.read(n); e={'name':n,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
                    if n.lower().endswith('.xlsx'):e['content']=xlsx(b)
                    else:
                        try:e['text']=b.decode('utf-8')
                        except UnicodeDecodeError:e['status']='binary member not decoded'
                    entries.append(e)
            binary.append({'path':item['path'],'entries':entries});item['status']='all ZIP members read; see binary_inventory.json'
        elif p.suffix=='.gz':
            text=gzip.decompress(blob).decode('utf-8'); blocks=[]; current=None
            for line in text.splitlines():
                if line.startswith('^'):
                    current={'entity':line,'metadata':[]};blocks.append(current)
                elif line.startswith('!') and current is not None:current['metadata'].append(line)
            soft.append({'path':item['path'],'lines':len(text.splitlines()),'blocks':blocks});item['status']='entire gzip decompressed/read; all SOFT metadata retained; expression rows not evaluated'
        elif p.suffix.lower() in ['.png','.jpg','.jpeg']:
            binary.append({'path':item['path'],'signature_hex':blob[:32].hex(),'status':'image bytes hashed; visual interpretation not performed by this tool'});item['status']='image read, metadata only'
        else:
            text=blob.decode('utf-8-sig');texts.append('\n### '+item['path']+'\n'+ '\n'.join(f'{i}: {l}' for i,l in enumerate(text.splitlines(),1)))
            item['status']='entire text read';item['lines']=len(text.splitlines())
            if p.suffix=='.csv' and 'prediction' in p.name:predictions.append(prediction_check(p,text))
    except Exception as exc:item['status']='UNREADABLE: '+repr(exc)
    coverage.append(item)
save('coverage.json',json.dumps(coverage,indent=2));save('text_evidence.txt','\n'.join(texts));save('binary_inventory.json',json.dumps(binary,indent=2));save('soft_metadata.json',json.dumps(soft,indent=2));save('saved_prediction_checks.json',json.dumps(predictions,indent=2))
print(json.dumps({'files':len(coverage),'text_files':len(texts),'binary_files':len(binary),'gzip_files':len(soft),'predictions':predictions},indent=2))
