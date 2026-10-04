"""Descriptive reconciliation only; uses table_checks audit functions, never project code."""
import csv,json,re,hashlib
from pathlib import Path
from table_checks import metrics
from collections import defaultdict
OUT=Path(__file__).parent
checks=json.loads((OUT/'TABLE_CHECKS.json').read_text(encoding='utf-8'))
recon=[]; claims=[]
fields=['source','heading','line','model','cohort','endpoint','features','threshold','validation','counts','source_reported','recomputed_archived','status','paired_fold_comparability','calibration_DCA_limits']
def claim(**kw):
 d=dict.fromkeys(fields,'unknown'); d.update(paired_fold_comparability='unknown across experiment families; same ID/fold labels are not original run provenance',calibration_DCA_limits='Brier is not calibration proof; no portal OOF rows, no validated DCA utility'); d.update(kw); claims.append(d)
for name,pred in [('model_ensemble_results','model_ensemble_predictions'),('data_augmentation_results','data_augmentation_predictions')]:
 table=next(t for t in checks if t['path']==f'reports/{pred}.csv')
 with open(f'reports/{name}.csv',encoding='utf-8',newline='') as f:
  for line,r in enumerate(csv.DictReader(f),2):
   g=next(g for g in table['groups'] if all(r[k]==v for k,v in g['group'].items()))
   computed=g['fold_means_explicit_0_5_scenario']; diffs={k:float(v)-computed[k[:-5]] for k,v in r.items() if k.endswith('_mean') and k[:-5] in computed}
   recon.append({'source':name,'group':g['group'],'differences_reported_minus_fold_mean':diffs,'max_abs_difference':max(map(abs,diffs.values()))})
   claim(source=f'reports/{name}.csv',line=line,model=r['model'],cohort='GSE94533',endpoint='Stage I/II vs benign/control; borderline excluded',features='fold-selected Top-5 + age + CA125',threshold='0.5 explicit scenario from source audit, not stored row threshold',validation='5 outer folds; arithmetic fold means',counts=json.dumps(g['label_counts']),source_reported=json.dumps(r),recomputed_archived=json.dumps(computed),status='matches archived probabilities at 0.5 scenario' if max(map(abs,diffs.values()))<1e-10 else 'difference')
phase=json.loads(Path('reports/gse106817_phase5_results.json').read_text())
pt=next(t for t in checks if t['path']=='reports/gse106817_phase5_predictions.csv')
for r in phase['experiments']:
 g=next(g for g in pt['groups'] if all(str(r[k])==v for k,v in g['group'].items() if k!='threshold') and abs(float(g['group']['threshold'])-r['threshold'])<1e-10)
 m=g['descriptive_metrics']; delta={k:r[k]-m[k] for k in ('n','positive','negative','roc_auc','pr_auc','accuracy','balanced_accuracy','sensitivity','specificity','brier_score')}
 recon.append({'source':'phase5','group':g['group'],'differences':delta,'confusion_matches':all(r['confusion_matrix'][k]==m[k] for k in ('tn','fp','fn','tp'))})
 claim(source='reports/gse106817_phase5_results.json',heading='experiments',model='source logistic / '+r['normalization'],cohort=r['train_dataset']+' -> '+r['test_dataset'],endpoint='invasive cancer vs benign/control; title-derived external negatives exclude borderline',features='14 historical shared miRNAs',threshold=r['threshold'],validation='separated bidirectional source fit; within_dataset modes transductive',counts=f"{r['positive']} positive/{r['negative']} negative",source_reported=json.dumps(r),recomputed_archived=json.dumps(m),status='matches saved grouped rows' if max(map(abs,delta.values()))<1e-10 else 'difference')
# Every substantive report line is inventoried, including metadata prose with no row-level proof.
for p in sorted(Path('reports').glob('*.md')):
 heading=''
 for line,text in enumerate(p.read_text(encoding='utf-8').splitlines(),1):
  if text.startswith('#'): heading=text
  elif text.strip() and not re.fullmatch(r'[| :\-]+',text):
   claim(source=p.as_posix(),heading=heading,line=line,source_reported=text,status='source-reported; see REPORT and RECONCILIATION for bounded verification, not independently reproduced')
# Summary-only CSVs retained as claims, not represented as recomputed OOF evidence.
for name in ('clinical_ablation_results','feature_selection_results','error_analysis_subgroups','2024_adnexal_mass_compatibility','dataset_compatibility','2024_adnexal_mass_feature_map'):
 with open(f'reports/{name}.csv',encoding='utf-8',newline='') as f:
  for line,r in enumerate(csv.DictReader(f),2): claim(source=f'reports/{name}.csv',line=line,model=r.get('model',r.get('selection_method','unknown')),endpoint=r.get('endpoint',r.get('task',r.get('analysis','unknown'))),cohort=r.get('dataset','GSE94533' if 'results' in name or 'subgroups' in name else '2024 literature cohort'),features=r.get('signature_size',r.get('feature_space','unknown')),source_reported=json.dumps(r),status='summary-only; structural/range checks, original OOF recomputation incapable')
with (OUT/'CLAIMS.csv').open('x',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(claims)
# Subgroups binary metrics are 0.5 source behavior, saved classes may differ.
with open('reports/error_analysis_predictions.csv',encoding='utf-8',newline='') as f: rows=list(csv.DictReader(f))
errors=[]
for model in sorted({r['model'] for r in rows}):
 rs=[r for r in rows if r['model']==model]; a=metrics(rs,'predicted_class'); b=metrics(rs,threshold=0.5)
 errors.append({'model':model,'saved_decision_metrics':a,'explicit_0_5_metrics':b,'decision_disagreements':sum(int(r['predicted_class'])!=int(float(r['predicted_probability'])>=.5) for r in rs),'unique_ids':len({r['sample_id_normalized'] for r in rs}),'missing_age':sum(not r['age'] for r in rs),'missing_ca125':sum(not r['ca125'] for r in rs),'false_negative_histology':dict(__import__('collections').Counter(r['histology'] for r in rs if r['true_label']=='1' and r['predicted_class']=='0'))})
with (OUT/'RECONCILIATION.json').open('x',encoding='utf-8') as f: json.dump({'numeric_reconciliation':recon,'error_models_separated':errors},f,indent=2)
print('claims',len(claims),'comparisons',len(recon))
print('max differences',max(x.get('max_abs_difference',max(map(abs,x.get('differences',{'a':0}).values()))) for x in recon))
for e in errors: print(e['model'],e['decision_disagreements'],e['saved_decision_metrics']['fn'],e['saved_decision_metrics']['fp'],e['explicit_0_5_metrics']['balanced_accuracy'])
