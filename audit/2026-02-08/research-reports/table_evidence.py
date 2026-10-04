"""Offline stdlib table evidence: no project imports or model/formula execution."""
import csv, hashlib, json, math
from collections import Counter
from pathlib import Path

def inspect_table(path):
    raw=Path(path).read_bytes()
    rows=list(csv.reader(raw.decode('utf-8-sig').splitlines()))
    header=rows.pop(0)
    columns={}
    for i,h in enumerate(header):
        vals=[r[i] if i<len(r) else '' for r in rows]
        missing=Counter(v for v in vals if v.strip().lower() in ('','na','nan','n/a','null','none') or v.lower().startswith('unknown'))
        nums=[]
        for v in vals:
            try:
                n=float(v)
                if math.isfinite(n): nums.append(n)
            except ValueError: pass
        counts=Counter(vals)
        columns[h]={'semantic_missing_tokens':dict(missing),'unique':len(counts),'numeric_n':len(nums),'min':min(nums) if nums else None,'max':max(nums) if nums else None,'values':dict(counts) if len(counts)<=35 else None}
    return {'path':str(path),'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw),'row_count':len(rows),'columns':header,'widths':dict(Counter(map(len,rows))),'duplicate_header_names':[h for h,n in Counter(header).items() if n>1],'exact_duplicate_rows':len(rows)-len(set(map(tuple,rows))),'column_checks':columns},[dict(zip(header,r)) for r in rows]

def auc(rows,score):
    pos=[float(r[score]) for r in rows if r['diagnosis']=='Cancer Early']
    neg=[float(r[score]) for r in rows if r['diagnosis']!='Cancer Early']
    return sum((a>b)+.5*(a==b) for a in pos for b in neg)/(len(pos)*len(neg))

if __name__=='__main__':
    import sys
    print(json.dumps(inspect_table(sys.argv[1])[0],indent=2))
