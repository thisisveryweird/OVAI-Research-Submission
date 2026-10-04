"""Compare every displayed numeric metric row to its archived summary row."""
import csv,json
from pathlib import Path
O=Path(__file__).parent
pairs={'clinical_ablation_results':('clinical_ablation_results',['endpoint','pipeline','model','imputation'],['roc_auc_mean','balanced_accuracy_mean','sensitivity_mean','specificity_mean','brier_score_mean']), 'feature_selection_stability':('feature_selection_results',['task','selection_method','signature_size'],['roc_auc_mean','balanced_accuracy_mean','sensitivity_mean','specificity_mean','accuracy_mean']), 'model_ensemble_comparison':('model_ensemble_results',['model'],['roc_auc_mean','pr_auc_mean','accuracy_mean','balanced_accuracy_mean','sensitivity_mean','specificity_mean','brier_score_mean']), 'data_augmentation_results':('data_augmentation_results',['model','augmentation'],['roc_auc_mean','pr_auc_mean','accuracy_mean','balanced_accuracy_mean','sensitivity_mean','specificity_mean','brier_score_mean','geometric_mean_mean'])}
out=[]
for md,(csvname,keys,cols) in pairs.items():
 with open('reports/'+csvname+'.csv',newline='') as f: rows=list(csv.DictReader(f))
 for line,text in enumerate(Path('reports/'+md+'.md').read_text(encoding='utf-8').splitlines(),1):
  if not text.startswith('|'):continue
  cells=[s.strip() for s in text.strip('|').split('|')]
  matching=[r for r in rows if [r[k] for k in keys]==cells[:len(keys)]]
  if not matching:continue
  r=matching[0]; diffs=[]
  for col,cell in zip(cols,cells[len(keys):]):
   diff=abs(float(r[col])-float(cell)); tolerance=.5*10**(-len(cell.split('.')[-1]))+1e-12
   diffs.append({'metric':col,'display':cell,'stored':r[col],'rounding_match':diff<=tolerance})
  out.append({'report':md,'line':line,'key':{k:r[k] for k in keys},'checks':diffs})
# Stability text expressions are independently checked against all table rows.
stability=[]
for md,tab in [('clinical_ablation_results','clinical_ablation_stability'),('feature_selection_stability','feature_selection_stability')]:
 import re
 with open('reports/'+tab+'.csv') as f: rows=list(csv.DictReader(f))
 task=method=None
 for line,text in enumerate(Path('reports/'+md+'.md').read_text().splitlines(),1):
  if text.startswith('### '):task=text[4:]
  if text.startswith('#### '):method=text[5:]
  if not text.startswith('- ') or '(5/5)' not in text and not re.search(r'\([1-4]/5\)',text):continue
  if md=='clinical_ablation_results':
   pipeline,model,imputation=text[2:].split(':')[0].split('/');rs=[r for r in rows if r['endpoint']==task and r['pipeline']==pipeline and r['model']==model and r['imputation']==imputation]
  else:
   size=re.search(r'Top (\d+)',text)[1];rs=[r for r in rows if r['task']==task and r['selection_method']==method and r['signature_size']==size]
  lookup={r['feature']:int(r['selected_folds']) for r in rs}
  for feature,n in re.findall(r'([^;: ]+) \((\d)/5\)',text):stability.append({'report':md,'line':line,'feature':feature,'claimed':int(n),'stored':lookup.get(feature),'match':lookup.get(feature)==int(n)})
with (O/'DISPLAY_CHECKS.json').open('x') as f:json.dump({'numeric_rows':out,'stability_expressions':stability},f,indent=2)
print(len(out),'numeric rows',sum(not c['rounding_match'] for r in out for c in r['checks']),'mismatches;',len(stability),'stability expressions',sum(not c['match'] for c in stability),'mismatches')
