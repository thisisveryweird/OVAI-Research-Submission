"""Read-only stdlib saved-decision audit. No project imports or threshold borrowing."""
import csv, hashlib, json, sys
from pathlib import Path
from collections import defaultdict

def audit(path):
    raw=path.read_bytes(); rows=list(csv.DictReader(raw.decode('utf-8-sig').splitlines()))
    fields=list(rows[0]) if rows else []
    keys=[k for k in ('experiment','train_dataset','test_dataset','normalization','threshold','dataset','model','augmentation') if k in fields]
    groups=defaultdict(list)
    for r in rows: groups[tuple(r[k] for k in keys)].append(r)
    output=[]
    for key, rs in groups.items():
        result={'group':dict(zip(keys,key)), 'n':len(rs)}
        decision=next((k for k in ('prediction','predicted_class') if k in fields),None)
        if decision:
            counts={'tn':0,'fp':0,'fn':0,'tp':0}; mismatch=0
            for r in rs:
                y=int(r['true_label']); d=int(r[decision]); assert y in (0,1) and d in (0,1)
                counts[{(0,0):'tn',(0,1):'fp',(1,0):'fn',(1,1):'tp'}[y,d]]+=1
                if 'threshold' in r: mismatch+=d!=int(float(r['probability'])>=float(r['threshold']))
            tn,fp,fn,tp=(counts[k] for k in ('tn','fp','fn','tp'))
            result.update(counts,sensitivity=tp/(tp+fn) if tp+fn else None,specificity=tn/(tn+fp) if tn+fp else None,threshold_decision_mismatches=mismatch,decision_column=decision)
        else: result['status']='No saved decision or row threshold; operating metrics unknown. No OOF threshold borrowed.'
        output.append(result)
    return {'path':str(path),'sha256':hashlib.sha256(raw).hexdigest(),'group_keys':keys,'groups':output}
if __name__=='__main__':
    print(json.dumps({'command':'python '+' '.join(sys.argv),'positive_class':1,'audits':[audit(Path(p)) for p in sys.argv[1:]]},indent=2))
