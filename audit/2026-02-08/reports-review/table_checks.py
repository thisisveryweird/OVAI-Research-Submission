"""Offline stdlib archived-table checks. No fitting or project imports.
Usage: python table_checks.py REPORT_DIRECTORY NEW_OUTPUT_DIRECTORY
Refuses overwriting outputs. Patient independence cannot be inferred from IDs.
"""
import csv, io, json, hashlib, math, sys
from pathlib import Path
from collections import Counter, defaultdict

def metrics(rs, decision=None, threshold=None):
    pc='probability' if 'probability' in rs[0] else 'predicted_probability'
    pairs=[(float(r[pc]),int(r['true_label'])) for r in rs]
    pos=sum(y for p,y in pairs); neg=len(rs)-pos
    sortedpairs=sorted(pairs); ranks=0; i=0
    while i<len(pairs):
        j=i+1
        while j<len(pairs) and sortedpairs[j][0]==sortedpairs[i][0]: j+=1
        ranks+=sum(y for p,y in sortedpairs[i:j])*(i+1+j)/2; i=j
    auc=(ranks-pos*(pos+1)/2)/(pos*neg) if pos*neg else None
    tp=fp=0; ap=0; i=0; descending=sorted(pairs,reverse=True)
    while i<len(pairs):
        j=i+1
        while j<len(pairs) and descending[j][0]==descending[i][0]: j+=1
        add=sum(y for p,y in descending[i:j]); tp+=add; fp+=j-i-add
        if pos: ap+=add/pos*tp/(tp+fp)
        i=j
    out={'n':len(rs),'positive':pos,'negative':neg,'roc_auc':auc,'pr_auc':ap if pos else None,'brier_score':sum((p-y)**2 for p,y in pairs)/len(rs)}
    if decision or threshold is not None:
        c=Counter((int(r['true_label']),int(r[decision]) if decision else int(float(r[pc])>=threshold)) for r in rs)
        tn,fp,fn,tp=(c[k] for k in ((0,0),(0,1),(1,0),(1,1)))
        se=tp/pos if pos else None; sp=tn/neg if neg else None
        out.update(tn=tn,fp=fp,fn=fn,tp=tp,sensitivity=se,specificity=sp,accuracy=(tn+tp)/len(rs),balanced_accuracy=(se+sp)/2 if se is not None and sp is not None else None,geometric_mean=math.sqrt(se*sp) if se is not None and sp is not None else None)
    else: out['operating_metrics']='incapable: no saved decision or threshold; no borrowed threshold'
    return out

def check(path):
    raw=path.read_bytes(); reader=csv.DictReader(io.StringIO(raw.decode('utf-8-sig'))); rows=list(reader); fields=reader.fieldnames
    result={'path':path.as_posix(),'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw),'rows':len(rows),'schema':fields,'duplicate_exact_rows':len(rows)-len({tuple(r.get(k) for k in fields) for r in rows}),'columns':{}}
    for k in fields:
        vals=[r.get(k) for r in rows]; freq=Counter(str(v) for v in vals); nums=[]
        for v in vals:
            try: nums.append(float(v))
            except (ValueError,TypeError): pass
        finite=[v for v in nums if math.isfinite(v)]
        col={'missing':sum(v is None or str(v).strip().lower() in ('','na','nan','none','null','-') for v in vals),'unique':len(freq),'numeric_count':len(nums),'nonfinite':len(nums)-len(finite)}
        if len(freq)<=100 or k in ('true_label','binary_label','selected_folds','selection_frequency','fold','folds'): col['frequency']=dict(freq)
        if finite:
            col.update(min=min(finite),max=max(finite))
            if any(t in k for t in ('probability','auc','accuracy','sensitivity','specificity','frequency','brier','ks_statistic','ks_pvalue')): col['outside_0_1']=sum(not 0<=v<=1 for v in finite)
        result['columns'][k]=col
    keys=[k for k in ('experiment','train_dataset','test_dataset','normalization','dataset','task','endpoint','pipeline','selection_method','signature_size','model','augmentation','imputation','threshold') if k in fields]
    groups=defaultdict(list)
    for r in rows: groups[tuple(r[k] for k in keys)].append(r)
    summaries=[]
    for key,rs in groups.items():
        g={'group':dict(zip(keys,key)),'rows':len(rs)}
        idcol=next((k for k in ('sample_id_normalized','sample_index','feature','miRNA','dataset') if k in fields),None)
        if idcol: g.update(id_column=idcol,unique_ids=len({r[idcol] for r in rs}),duplicate_ids=len(rs)-len({r[idcol] for r in rs}))
        if 'fold' in fields: g['fold_counts']=dict(Counter(r['fold'] for r in rs))
        if 'true_label' in fields:
            g['label_counts']=dict(Counter(r['true_label'] for r in rs))
            dec=next((k for k in ('prediction','predicted_class') if k in fields),None)
            g['descriptive_metrics']=metrics(rs,dec)
            if 'fold' in fields:
                # 0.5 is an explicit audit scenario, not an inferred original threshold.
                fs=[metrics([r for r in rs if r['fold']==f],threshold=0.5) for f in sorted({r['fold'] for r in rs})]
                g['fold_metrics_explicit_0_5_scenario']=fs
                g['fold_means_explicit_0_5_scenario']={k:sum(x[k] for x in fs)/len(fs) for k in fs[0] if isinstance(fs[0][k],(int,float)) and all(x[k] is not None for x in fs)}
        if 'selected_folds' in fields:
            folds=5  # protocol denominator; stability 'folds' stores selected fold IDs, not denominator
            if 'folds' in fields:
                g['fold_list_count_mismatches']=sum(len([v for v in r['folds'].split(',') if v])!=int(r['selected_folds']) for r in rs)
            g.update(selected_total=sum(int(r['selected_folds']) for r in rs),nonzero_features=sum(int(r['selected_folds'])>0 for r in rs),frequency_mismatches=sum(abs(float(r['selection_frequency'])-int(r['selected_folds'])/folds)>1e-10 for r in rs),selected_outside_fold_range=sum(not 0<=int(r['selected_folds'])<=folds for r in rs),top_features=[{'feature':r['feature'],'selected_folds':int(r['selected_folds'])} for r in sorted(rs,key=lambda r:(-int(r['selected_folds']),r['feature']))[:15]])
        if all(k in fields for k in ('tn','fp','fn','tp','n')):
            g['confusion_count_mismatches']=sum(sum(float(r[k]) for k in ('tn','fp','fn','tp'))!=float(r['n']) for r in rs if all(r[k] for k in ('tn','fp','fn','tp','n')))
        summaries.append(g)
    result['group_keys']=keys; result['groups']=summaries
    return result

if __name__=='__main__':
    src,out=map(Path,sys.argv[1:]); out.mkdir(parents=True,exist_ok=True)
    checks=[check(p) for p in sorted(src.glob('*.csv'))]
    with (out/'TABLE_CHECKS.json').open('x',encoding='utf-8') as f: json.dump(checks,f,indent=2,ensure_ascii=False)
    with (out/'COVERAGE.csv').open('x',encoding='utf-8',newline='') as f:
        w=csv.writer(f); w.writerow(['path','sha256','bytes','lines_or_rows','full_read','analysis_type','gap'])
        for p in sorted(src.iterdir()):
            if p.suffix not in ('.csv','.md','.json'): continue
            raw=p.read_bytes(); text=raw.decode('utf-8-sig')
            if p.suffix=='.json': json.loads(text)
            checkrow=next((x for x in checks if x['path']==p.as_posix()),None)
            w.writerow([p.as_posix(),hashlib.sha256(raw).hexdigest(),len(raw),checkrow['rows'] if checkrow else len(text.splitlines()),True,'all-row schemas/missingness/duplicates/ranges/group metrics/stability' if checkrow else 'full text/JSON claim review','original inputs/execution provenance not established'])
