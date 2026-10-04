"""Read-only prior evidence inventory and metadata eligibility reconstruction."""
import csv,json,re,hashlib
from pathlib import Path
root=Path('audit/2026-02-08/external-final')
prior=Path('audit/2026-02-08/external')
coverage=json.loads((prior/'coverage.json').read_text())
substantive={'scripts/bulk_test_external_datasets.py','scripts/evaluate_gse106817_external.py','scripts/evaluate_gse106817_phase5.py','reports/bulk_external_test_results.md','research/cohort_expansion/reports/GSE211692_external_validation_audit.md','research/cohort_expansion/reports/2024_adnexal_mass_partial_signature_replication.md','research/cohort_expansion/reports/2024_adnexal_mass_published_model_reproduction.md','research/cohort_expansion/reports/2022_gockley_score_semantics_audit.md'}
with (root/'COVERAGE.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.writer(f);w.writerow(['path','sha256','bytes','prior_ingestion_status','final_substantive_read','gap'])
 for r in coverage:
  p=r['path'].replace('\\','/');s=p in substantive
  w.writerow([p,r['sha256'],r['bytes'],r['status'],'complete text substantive read' if s else 'ingestion only / evidence reuse','' if s else 'substantive full review not complete; binary visuals not rendered' if not p.endswith(('.py','.md','.json','.csv')) else 'substantive full review gap'])
soft=json.loads((prior/'soft_metadata.json').read_text())
out=[]
for dataset in soft:
 counts={'positive':0,'negative':0,'excluded':0};samples=[]
 for b in dataset['blocks']:
  if not b['entity'].startswith('^SAMPLE'):continue
  selected=[x.split(' = ',1)[-1] for x in b['metadata'] if any(x.startswith('!Sample_'+k+' =') for k in ('title','source_name_ch1','characteristics_ch1','description'))]
  value=' '.join(selected).lower()
  y=None
  if 'borderline' not in value:
   if re.search(r'ovarian\s+(cancer|carcinoma)|ovarian\s+tumou?r',value):y=1
   elif re.search(r'\b(benign|healthy|control|non[- ]?cancer|normal)\b',value):y=0
  counts['positive' if y==1 else 'negative' if y==0 else 'excluded']+=1
 out.append({'path':dataset['path'],'bulk_parser_label_counts_before_feature_eligibility':counts,'note':'Feature eligibility not verified; mirrors metadata regex without importing project.'})
(root/'metadata_eligible_counts.json').write_text(json.dumps(out,indent=2))
