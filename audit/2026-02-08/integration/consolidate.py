"""Offline audit consolidation only; never imports or runs product modules."""
import csv,json,re,hashlib,subprocess
from pathlib import Path
ROOT=Path.cwd(); BASE=ROOT/'audit/2026-02-08'; OUT=BASE/'integration'
AREAS='inference-final governance-final external-final core-final discovery-final historical-reconstruction historical-identity builders-review transfer-stage-review remaining-scripts research-sources reports-review research-reports'.split()
def text(p): return p.read_text(encoding='utf-8-sig')
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def put(n,s): (OUT/n).write_text(s,encoding='utf-8')
def dump(n,x): put(n,json.dumps(x,ensure_ascii=False,indent=2)+'\n')
def table(n,rows):
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with (OUT/n).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,keys); w.writeheader();w.writerows(rows)
reports={a:text(BASE/a/'REPORT.md') for a in AREAS}
coverage=[];endpoints=[]
for a in AREAS:
 for name in ['COVERAGE.csv','ENDPOINTS.csv','LABEL_ENDPOINT_REGISTRY.csv']:
  p=BASE/a/name
  if not p.exists(): continue
  with p.open(encoding='utf-8-sig',newline='') as f: rows=list(csv.DictReader(f))
  for i,r in enumerate(rows,2):
   r={'evidence_file':p.relative_to(ROOT).as_posix(),'evidence_row':i,**r}
   (coverage if name=='COVERAGE.csv' else endpoints).append(r)
# Retain variants and alias labels verbatim; do not collapse distinct clinical tasks.
for i,r in enumerate(endpoints,1):
 r['registry_id']=f'EP-{i:04d}';r['aliases_variant_labels']=' | '.join(str(r[k]) for k in r if any(t in k.lower() for t in ['endpoint','population','task','label','variant']))
 r['count_authority']='Source-reported unless independently decoded metadata/saved-row count is explicitly stated; never a new clinical adjudication.'
table('ENDPOINT_REGISTRY.csv',endpoints)
recon=BASE/'remaining-scripts/SCRIPT_COVERAGE_RECONCILIATION.csv'
with recon.open(encoding='utf-8-sig',newline='') as f:
 for i,r in enumerate(csv.DictReader(f),2): coverage.append({'evidence_file':recon.relative_to(ROOT).as_posix(),'evidence_row':i,'coverage_kind':'101-script substantive proof reconciliation',**r})
for p in (BASE/'pdf').iterdir():
 if p.is_file(): coverage.append({'evidence_file':p.relative_to(ROOT).as_posix(),'coverage_kind':'complete inventory/text machine ingestion; PDF rendered visuals unavailable','sha256_current':sha(p)})
table('READ_COVERAGE.csv',coverage)
put('EVIDENCE_REPORTS.md','# Finalized substantive evidence reports (inherited static review)\n\nIntegration read all supplied reports and entire coverage tables. This appendix preserves exact function/line evidence and limitations; no source was rerun. Earlier ingestion-only gaps are superseded only where later full-read proof exists.\n\n'+'\n\n'.join('## '+a+'\nSource: `audit/2026-02-08/'+a+'/REPORT.md`\n\n'+s for a,s in reports.items()))
# Stable inventory: each repeated page-local item is independent.
inv=text(BASE/'pdf/recommendation-inventory.md'); recs=re.findall(r'^- (PDF-P(\d+)-R\d+): (.*)$',inv,re.M)
themes=[
(['fold','folds','nested','feature-availability','availability filters','feature selection'],'partial','scripts/run_phase8_model_experiments.py','main','275-338','discovery-final','Authorize fold-contained preprocessing/selection/calibration and compare saved outer-fold rows under a locked protocol.'),
(['imputation','imputed','bulk','UI','portal','logging','confidence'],'partial','scripts/ovarian_diagnostic_portal.py','render_bulk_upload / main','157-178;222-252','inference-final','Approve parity fixtures and input/artifact snapshots, then verify flags and stale-result rejection without clinical claims.'),
(['calibration','Platt','probability','decision-curve','decision curves'],'partial','scripts/run_reproducible_baseline.py','calibration_values','278-313','core-final','Obtain model-specific OOF rows; preregister whole-pipeline calibration and DCA protocol, then validate on untouched data.'),
(['historical','original split'],'unresolved','scripts/search_possible_historical_test_splits.py','module','401-466','historical-identity','Obtain authenticated original 135/44 membership and preprocessing lineage; matching metrics is insufficient.'),
(['TCGA','staging'],'partial','scripts/train_tcga_ov_stage_model.py','load_features','79-136','core-final','Verify staging-only labels in every presentation; never substitute this endpoint for detection.'),
(['fusion','same-patient','same-visit','concatenate','scaffold','unit-tested'],'partial','src/data/multimodal_table.py','align_modalities','18-58','governance-final','Data owner must approve donor/visit crosswalk; add true linkage and all-absent guard tests after authorization. Synthetic smoke success is not established.'),
(['external','overlap','platform','harmonization','transfer'],'unresolved','scripts/bulk_test_external_datasets.py','evaluate / test_one','293-305;371-374','external-final','Freeze pipeline/threshold and assay bridge; obtain donor overlap adjudication and exact eligible roster before independent validation.'),
(['HE4','dbGaP','prospective','longitudinal','multi-site'],'needs_external_data_owner_decision','README.md','controlled-access discussion','101-113','governance-final','Obtain documented access/reuse/IRB and donor-linked visit/assay metadata; do not acquire or apply before authorization.'),
(['license','manifest','configuration','dependency','RESULTS_INDEX','canonical','guardrail','feature dimensionality','git'],'unresolved','configs/project.yaml','module','7-13;39-53','governance-final','Owner approves license/canonical paths and locked manifest; then verify source/input/output hashes, dimensions and fold guards.'),
]
rows=[]
for rid,page,desc in recs:
 status='needs_external_data_owner_decision';path='README.md';fn='intended-use / research scaffold';lines='9-25';area='governance-final';step='Obtain intended-use, clinical/data-owner and governance approval with primary evidence; validate a prespecified workflow/protocol after approval.'
 for words,st,pa,fu,li,ar,se in themes:
  if any(re.search(r'(?<!\w)'+re.escape(w.lower())+r'(?!\w)',desc.lower()) for w in words): status,path,fn,lines,area,step=st,pa,fu,li,ar,se;break
 if rid=='PDF-P18-R03':status='obsolete';path='.git/HEAD';fn='repository metadata';lines='1';step='No initialization/commit: current clone already has git. Any tagging/commit requires separate owner approval.'
 if rid=='PDF-P07-R05':status='unsupported';path='scripts/run_phase8_model_experiments.py';fn='main';lines='275-338';area='discovery-final';step='Remove blanket fully-nested/unit-tested assurance; validate each fold boundary and patient-linkage guard under separately approved tests.'
 rows.append({'recommendation_id':rid,'pdf_physical_page':page,'recommendation':desc,'disposition':status,'source_path':path,'function_or_scope':fn,'source_lines':lines,'final_evidence_report':f'audit/2026-02-08/{area}/REPORT.md','next_verifiable_step':step,'authorization':'false; pending','scope':'Current source supports only the cited bounded implementation; no clinical efficacy or recommended action completed.'})
table('PDF_TO_REPOSITORY_CHECKLIST.csv',rows)
# Check every citation range used in checklist exists, including metadata without loading.
checks=[]
for r in rows:
 p=ROOT/r['source_path'];n=len(text(p).splitlines()) if p.exists() else 0
 for span in r['source_lines'].split(';'):
  ns=list(map(int,span.split('-')));checks.append({'id':r['recommendation_id'],'path':r['source_path'],'range':span,'actual_lines':n,'exists_and_in_range':bool(n and 1<=ns[0]<=ns[-1]<=n)})
dump('CITATION_CHECKS.json',checks)
# preservation vs initial governance inventory
initial=json.loads(text(BASE/'governance/inventory.json'));pres=[]
for r in initial:
 p=ROOT/r['path']
 if r['path'].replace('\\','/').startswith(('.git/','audit/')): continue
 pres.append({'path':r['path'],'initial_sha256':r['sha256'],'current_sha256':sha(p) if p.is_file() else None,'unchanged':p.is_file() and sha(p)==r['sha256']})
dump('PRESERVATION_CHECKS.json',pres)
commands=[]
for c in [['git','--no-optional-locks','status','--short'],['git','--no-optional-locks','diff','--stat'],['git','--no-optional-locks','diff','--name-only']]:
 x=subprocess.run(c,capture_output=True,text=True,encoding='utf-8',errors='replace'); commands.append({'command':' '.join(c),'exit_code':x.returncode,'stdout':x.stdout,'stderr':x.stderr})
# Own whole-table grouping, never pool repeated experiments into a performance metric.
def archived_auc(pairs):
 pairs=sorted(pairs);n1=sum(y==1 for s,y in pairs);n0=sum(y==0 for s,y in pairs)
 if not n1 or not n0:return None
 rank_sum=0.;i=0
 while i<len(pairs):
  j=i+1
  while j<len(pairs) and pairs[j][0]==pairs[i][0]:j+=1
  rank_sum+=sum(y==1 for s,y in pairs[i:j])*((i+1+j)/2);i=j
 return (rank_sum-n1*(n1+1)/2)/(n1*n0)
groups=[]
for directory in [ROOT/'reports',ROOT/'research/cohort_expansion/reports',ROOT/'research/transportability_analysis']:
 if not directory.exists():continue
 for p in directory.glob('*.csv'):
  with p.open(encoding='utf-8-sig',newline='') as f:
   reader=csv.DictReader(f);headers=reader.fieldnames or [];rs=list(reader)
  labels=[k for k in headers if k.lower() in ['true_label','label','y_true','binary_label']]
  if not labels:continue
  keys=[k for k in headers if k.lower() in ['model','experiment','direction','mode','normalization','train_dataset','test_dataset','augmentation','configuration','split','endpoint','threshold']]
  accum={}; group_rows={}
  for r in rs:
   key=tuple(r.get(k,'') for k in keys);g=accum.setdefault(key,{'rows':0,'class_counts':{}});g['rows']+=1;v=r.get(labels[0],'');g['class_counts'][v]=g['class_counts'].get(v,0)+1;group_rows.setdefault(key,[]).append(r)
  prob=next((k for k in headers if k.lower() in ['probability','predicted_probability','y_prob','y_score','cancer_probability']),None)
  decision=next((k for k in headers if k.lower() in ['prediction','predicted_class','y_pred']),None)
  for key,gr in group_rows.items():
   g=accum[key];g['probability_column']=prob;g['decision_column']=decision
   if prob:
    try:
     pairs=[(float(r[prob]),int(float(r[labels[0]]))) for r in gr];g['within_experiment_archived_row_auc']=archived_auc(pairs);g['aggregation_note']='Archived within-experiment row AUC; pooled OOF only where rows are OOF, external ranking where target held out. Not mean-fold AUC.'
     folds={r.get('fold') for r in gr}
     if 'fold' in headers:
      values=[archived_auc([(float(r[prob]),int(float(r[labels[0]]))) for r in gr if r['fold']==f]) for f in folds];valid=[v for v in values if v is not None];g['mean_fold_auc']=sum(valid)/len(valid) if valid else None
    except ValueError as e:g['metric_error']=str(e)
   if decision:
    try:
     c={'tn':0,'fp':0,'fn':0,'tp':0}
     for r in gr:
      y=int(float(r[labels[0]]));d=int(float(r[decision]));c[{(0,0):'tn',(0,1):'fp',(1,0):'fn',(1,1):'tp'}[(y,d)]]+=1
     g['saved_decision_confusion']=c;g['sensitivity']=c['tp']/(c['tp']+c['fn']) if c['tp']+c['fn'] else None;g['specificity']=c['tn']/(c['tn']+c['fp']) if c['tn']+c['fp'] else None
    except (ValueError,KeyError) as e:g['metric_error']=str(e)
  groups.append({'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p),'group_columns':keys,'label_column':labels[0],'groups':[{'key':list(k),**g} for k,g in accum.items()],'metric_status':'Own saved-row arithmetic within each experiment only; within-experiment pooled OOF and mean-fold AUC explicitly distinct. No pooling across experiments and no model execution.'})

dump('ROW_CLASS_GROUP_CHECKS.json',groups)
# Canonical evidence registry includes all current endpoint families and all report claims, not just README.
put('RESULTS_INDEX_DRAFT.md','# Results index — draft, source/output reported NOT rerun\n\nAll original input hashes/releases, original package/fit metadata and canonical owner decisions are unknown. Current file hashes are in READ_COVERAGE; availability is not original lineage. No command below was executed. Network must be disabled and inputs pre-staged; live Xena/acquisition paths have NO approved offline command. Endpoint registry retains all substantive variants and reported versus actual counts.\n\n'+text(BASE/'governance-final/RESULTS_INDEX_DRAFT.md')+'\n\n## Comprehensive current claim-family protocol index\n\nEvery row of ENDPOINT_REGISTRY.csv is an indexed current endpoint/variant; exact source/function/line protocol and input state are retained verbatim. Every root report claim is additionally indexed in reports-review/CLAIMS.csv (969 claims), and research endpoints/counts in research-reports/LABEL_ENDPOINT_REGISTRY.csv (422 entries). These are source/output-reported evidence, not new fits.\n\n| Family | Offline command/protocol availability | Evidence and blocker |\n|---|---|---|\n| Phase6 selection, Phase7 clinical ablation, Phase8 ensemble, Phase9 augmentation | Candidate future `python scripts/run_phase6_feature_selection.py`; `python scripts/run_phase7_clinical_experiments.py`; `python scripts/run_phase8_model_experiments.py`; `python scripts/run_phase9_data_augmentation.py` only after owner verifies actual filenames/input paths and pins all inputs | Saved group-specific tables/rows; fold-mean vs pooled distinct; source contamination and discovery decision provenance unresolved |\n| Frozen/engineered, pairwise, rank | Candidate `python scripts/test_high_capacity_models.py`; pairwise/rank command blocked pending canonical script selection | No matching high-capacity OOF rows; no shipped-artifact provenance or external nine-feature validation |\n| Historical PMML reconstruction, Table2 agreement and split search | NO approved reproduction command until primary inputs and original135/44 recovered | 119+60 vs120+59 conflict; tanh scores not probabilities; exploratory matching not validation |\n| Fair/full/modern14, age/CA125/192 and clinical builders | NO approved canonical command; variants retained in builders-review endpoint/source registry | Positional identities, 97/82 vs96/83 labels, TPM vs adjusted log10, absent inputs, upstream processing unknown |\n| GSE106817 target development; cancer/borderline; bulk transfer; phase5 bidirectional | NO approved canonical command; grouped archived output protocol retained in external-final and transfer-stage-review | Target development != frozen external validation; all-positive bulk != all-negative historical; target-distribution adaptation distinct |\n| GSE8841/12172/9899 tissue StageI/LMP/generic malignant | NO approved canonical command until label/assay/roster lock | Source assertions30/16 vs30/60, target-z exploratory opposite orientation; not serum detection |\n| TCGA paired304, multimodal223, stage299 and survival | No supported offline command for live-query loaders | Distinct Ns/endpoints, global prep and inner scaling; pooled hazard C-index not fold mean; tumor only |\n| Gockley275 PMML and strict72/published75 subsets | No approved network replay; archived arithmetic checks only | 165/55/55 selection; AUC .8867 vs .87; calibrated probability/selected winner unknown |\n| FirePlex678/286/392 | No frozen reproduction command available | Cells/row hashes verified by research-reports; weights absent,0/5 portal overlap, specimen contradiction unresolved |\n| GSE21169216190 metadata, strict168/28, proposed168/200 | No validation command approved; matrix/control roster absent | SOFT counts verified; expression completeness/donor independence unproved |\n| Danish380; Janus80 people106 samples; other panels/scouts | No scoring command: access/raw data absent or endpoints incompatible | Source summaries, not usable independent patient validation |\n| Fusion/synthetic scaffold; CLI/portal/demo | No clinical test command; launchers install and are not offline-safe | Static code only; masks/linkage unsafe,7-vs9 CLI schema, stale results and imputation mismatch |\n\n## Exact inherited protocol and citation appendix\n\n'+ '\n\n'.join('### '+a+'\n'+s for a,s in reports.items()))
# Compact audit with exact citations.
put('AUDIT.md','''# Final integrated audit — AUDIT ONLY

Scope: 101 direct scripts, eight research Python sources, root/config/src text, current root/research reports and whole-table saved CSV analyses, safe binary metadata and decoded supplied workbook/SOFT evidence. Coverage is inherited hash-linked substantive review, not a fresh project execution. Complete inventory preserves repeated PDF recommendations. No product edits, model loading, imports, fits, installs, network, applications, commits or clinical tests. Folder date is a label: HEAD and PDF dates are later; chronology is not recovered.

## Ranked findings
1. **Clinical/deployed external validation not established.** CLI constructs seven inputs vs portal nine; evaluator multiplies before imputation while export/portal impute first. `scripts/predict_ovarian_risk.py::main:35–69`, `scripts/test_high_capacity_models.py::engineered_frame:93–98`, `scripts/export_production_pipeline.py::transformed_values:80–90`, `scripts/ovarian_diagnostic_portal.py::build_model_input:75–93`. Static binary string evidence is not loaded estimator assurance.
2. **Blanket fully-nested claim false; full discovery provenance unknown.** Phase8 selection/transformation precedes inner tuning (`scripts/run_phase8_model_experiments.py::main:275–338`); global eligibility exists (`scripts/evaluate_gse94533_diagnostic_model.py::main:62–66`); classifier-only calibration is not whole-pipeline contained (`scripts/run_reproducible_baseline.py::make_pipeline:88–108`). Fixed-panel outer CV is conditional, not nested discovery. Panel-name overlap is not proof of causal leakage.
3. **Endpoints, aggregations and external failures must remain separate.** High-capacity pooled OOF AUC is mislabeled mean (`scripts/test_high_capacity_models.py::run:117–143`; `reports/high_capacity_model_results.md:5`). Phase6–9 use fold means; these cannot be compared as identical estimands. Bulk all-positive (`reports/bulk_external_test_results.md:11–13`) differs from historical14 all-negative and bidirectional/adapted phase5 modes. Eligible N140/1058 differs from series972/1817. Grouped saved-row metrics are arithmetic verification only, not external nine-feature portal validation.
4. **Patient identity and historical split unproved.** Clean119+60 recovery contradicts exported120+59/E1056 exclusion (`scripts/reproduce_original_historical_clean.py::recover_historical_sample_order:470–671`, `save_mapping:1313–1327`, `main:1689–1694`). Positional attachment is not biological identity; 20,000 metric-compatible split search cannot recover original135/44 (`scripts/search_possible_historical_test_splits.py:401–466`). Tanh activations are not probabilities.
5. **Scaffold assurances and governance exceed proof.** Key alignment checks uniqueness, not donor/visit linkage (`src/data/multimodal_table.py::align_modalities:18–58`); all-absent masks are not rejected (`src/models/gated_fusion.py::forward:60–82`); synthetic smoke assertions do not establish unit-test success (`scripts/smoke_test.py::main:14–43`). IRB, reuse, license, regulations, reimbursement and biological causal mechanisms remain unverified.
6. **Research evidence bounded.** HE4 exists elsewhere (2024 workbook/Mendeley/MethylBERT schema), not an HE4-linked GSE94533 comparator. GSE211692 hardcoded zero overlaps (`research/cohort_expansion/audit_gse211692.py::make_overlap_table:342–378`) do not prove donors independent; actual SOFT title intersections zero also do not. Matrix and200-control roster absent.2024 workbook present; stale absence prose superseded, serum/plasma contradiction unresolved. Histology145 rows !=141 unique IDs;56early/84late/1conflict. Gockley archived AUC .8867 != published .87 and no probability map. See research-reports/REPORT.md and endpoint registry for exact counts and limitations.

## Explicit PDF reconciliation
No-git recommendation obsolete: `.git/HEAD:1`, current HEAD5fd00c9. Calibration absent-everywhere is false (`scripts/run_reproducible_baseline.py::calibration_values:278–313`), but production calibration/DCA not established. Fully-nested blanket claim false as above. Pooled-vs-mean AUC, different external pipelines/collapse directions, HE4 elsewhere not GSE-linked, biological preservation/platform causation, IRB/public-data exemption and time-sensitive regulations must not be promoted as verified. No model winner, intended-use tier or threshold chosen.

Legacy governance/inference/external/core-evaluation/discovery-review/review reports and `external/saved_prediction_checks.json` are NONAUTHORITATIVE SUPERSEDED preliminary evidence; the pooled saved prediction checks must not be cited as metrics. They remain untouched. PNG/PDF rendering unavailable; machine text extraction and safe metadata are not visual verification. All current authoritative bounded reports are copied in EVIDENCE_REPORTS.md; exact evidence ranges and statuses in coverage/checklist. Remaining blockers are primary provenance, missing raw/processed inputs, donor/visit linkage, approval and untouched external assay validation.
''')
# Approval tranches: exact existing paths from coverage and bounded read scopes.
tranches={
 'A inference parity and research-only UI':['scripts/test_high_capacity_models.py','scripts/export_production_pipeline.py','scripts/predict_ovarian_risk.py','scripts/ovarian_diagnostic_portal.py','scripts/demo_ovarian_cases.py','scripts/ovarian_research_demo.py'],
 'B fold-contained scientific evaluation':['scripts/evaluate_gse94533_diagnostic_model.py','scripts/evaluate_gse94533_enhancements.py','scripts/run_reproducible_baseline.py','scripts/evaluate_methylbert_methylation.py','scripts/train_tcga_ov_stage_model.py','scripts/run_phase8_model_experiments.py'],
 'C linkage and scaffold safeguards':['src/data/multimodal_table.py','src/models/gated_fusion.py','scripts/smoke_test.py','scripts/build_gse94533_clinical_fusion.py','scripts/add_age_ca125_to_14mirna.py','scripts/rebuild_clinical_benchmark.py'],
 'D claim/provenance-only corrections':['README.md','START_HERE.txt','configs/project.yaml','requirements.txt','requirements-deep.txt','Launch_Diagnostic_Portal.bat','Launch_Diagnostic_Portal.command','reports/bulk_external_test_results.md','reports/high_capacity_model_results.md'],
 'E external eligibility and research truthfulness':['scripts/bulk_test_external_datasets.py','scripts/evaluate_gse106817_phase5.py','research/cohort_expansion/audit_gse211692.py','research/cohort_expansion/audit_cross_platform_replication.py','research/transportability_analysis/run_transportability_analysis.py']}
put('PATCH_APPROVAL_REQUEST.md','# Patch approval request — authorization=false; ALL PENDING\n\nNo patches applied. No winner, threshold, license or tier selected. Each tranche is independently bounded: approval must specify exact paths, permitted behavior changes and exclusions; no acquisition/fitting/deployment implied.\n\n'+ '\n\n'.join('## '+k+'\n'+ '\n'.join('- `'+p+'`'+(' (exists)' if (ROOT/p).exists() else ' (NOT FOUND: excluded until corrected owner path)') for p in ps)+'\nAcceptance after approval: review diff; confirm protected files unchanged; use prespecified synthetic/archived fixtures and schema/fold/threshold/source hash guards; clinical validation requires separate protocol.' for k,ps in tranches.items())+'''\n\n## Genuine architecture alternatives (decision not made)
1. Recommended for consideration: separately versioned shared fold-safe sklearn Pipeline containing availability, imputation, transforms, interaction construction, selection and estimator; calibrate the entire pipeline within appropriate splits. New version and immutable manifests; preserve legacy artifacts; establish parity or explicitly different scientific estimand.
2. Retain legacy artifact via a narrow safety wrapper enforcing nine-feature schema, exact legacy transformation and immutable artifact hash; reject invalid/ambiguous/all-missing inputs, expose imputation and research-only outputs. Preserves compatibility but cannot repair training provenance or establish clinical validation; seven-input CLI must not silently adapt.
3. Supported replacement path: approved maintained assay-specific pipeline with authenticated units/linkage, pinned supported packages, migration protocol and prospective untouched validation. Legacy kept as reference; no automatic winner selection or threshold transfer.

Signatories: repository owner for code, scientific lead for estimand, data owner for linkage/reuse, institutional IRB/privacy reviewer, clinical lead and regulatory counsel for intended use. Evidence and postapproval tests are in BLOCKED_RESEARCH_AND_DECISIONS.md. No approval received.
''')
decisions=[('Intended use/headline/model/threshold','endpoint registry, comparator and prospective protocol','clinical + scientific leads','locked untouched cohort and prespecified utility tradeoffs'),('Canonical architecture and parity','nine-feature manifest, missingness fixtures, source hashes','repo owner + scientific lead','whole-pipeline fold guards, CLI/portal/export parity and explicit version migration'),('Data access dbGaP/HE4/longitudinal','source-specific terms, consent, eligibility, linked visit/assay inventory','data owner + DAC + institutional reviewer','authorized acquisition manifest, checksums and measured-versus-imputed audit'),('Donor independence and linkage','donor/visit/biobank crosswalk with recruitment windows; original135/44','data owner + scientific lead','zero-overlap subject-level audit and authenticated header/label reconciliation'),('Research metadata resolution','2024 primary article/correction specimen,101798 stage adjudication,GSE211 matrix/control roster,PMML winner/probability map','source authors/data owner + scientific lead','versioned cohort exclusions/counts and primary-to-derived hash checks'),('Calibration/DCA/subtypes/StageI','model-specific OOF rows, fold IDs, prespecified probability/utility protocol and subtype denominators','scientific + clinical leads','untouched calibration/decision-curve validation with representative prevalence'),('License/reuse/IRB/privacy/logging','rights inventory, institutional determination, lawful logging retention and minimization policy','repo owner + counsel + institutional IRB/privacy authority','approved notices and audit-log access/deletion controls; no assumed exemption'),('Regulatory/reimbursement/lab/hardware','current primary FDA/CMS/court/CLIA/payer/SGO/ACOG sources and assay intended-use plan','regulatory counsel + clinical/lab leads','dated applicability memo then authorized analytical/prospective validation'),('Funding/collaborator priorities','documented access/assay milestones, partner consent and actual quotes','owner/funding lead','signed milestones and feasibility evidence; no invented budget or partners')]
put('BLOCKED_RESEARCH_AND_DECISIONS.md','# Bounded pending decisions\n\nEvery approval is pending; authorization=false. No applications, acquisition, experiments or deployment authorized.\n\n'+ '\n\n'.join('## '+title+'\n- [ ] Evidence required: '+e+'.\n- [ ] Signatory: '+s+'; written bounded approval pending.\n- [ ] Postapproval validation: '+v+'.\n- [ ] Stop if evidence/approval absent; record unknown, not zero or passed.' for title,e,s,v in decisions))
put('DATA_SOURCES_DRAFT.md','# Data sources draft — reuse/consent/IRB unknown\n\nAvailability is not permission. No rights or exemption inferred from public access. Document original article/accession/version, acquisition date/hash, specimen/units, participant consent, allowed reuse/commercial use, donor/visit crosswalk and institutional determination before use. Author/contact/access routes remain proposals.\n\n'+text(BASE/'external-final/DATA_SOURCES_DRAFT.md')+'\n\n'+text(BASE/'research-reports/DATA_SOURCES_ADDENDUM.md')+'\n\n2024 serum/plasma conflict requires authenticated article/correction and collection map; HE4 cells elsewhere do not establish GSE94533 linkage. GSE211692 SOFT metadata is available, expression matrix/control roster not established; Danish/Janus access remains source-reported. Historical primary workbooks/split/processing absent. All owner/IRB/regulatory approvals pending. No new acquisition performed.\n')
put('DRAFT_PRESENTATION.md','''# Research-demo scorecard — not clinical validation (draft)

**Problem/intended use:** earlier workup is a hypothesis; tier, population and comparator not selected. Not average-risk screening or diagnosis. **Evidence:** current GSE94533 source reports full179 (96/83) versus early113 (52/61) endpoints; high-capacity archived engineered pooled OOF AUC .8997, sensitivity .8269, specificity .9016, balanced accuracy .8643. No matching saved OOF rows/artifact lineage or unbiased discovery estimate; not nine-feature external validation. Full14-miRNA+clinical reported AUC .813 is a different population/model. Methylation and tumor-stage results are independent endpoint benchmarks, not fused improvements.

**External failures:** bulk5-feature tables all-positive; historical14/phase5 raw all-negative; AUC ranking does not establish usable threshold/calibration or causal platform explanation. Report model/transform/threshold/N per group. **Gates:** missing primary inputs, donor overlap and visit linkage, imputation-order parity, original split, assay equivalence, calibration, prospective comparator and institutional approvals.

**Honest demo guide (proposed, not launched):** display archived scorecard and exact pipeline identity; label synthetic/example or stored OOF patients explicitly. Do not install via launchers. If later authorized, expose endpoint/N, research-only warning, measured/default flags, effective inputs, artifact hash and input snapshot; clear stale result on change/error. No individual medical advice or validated cancer-risk claim.

**PROPOSED workflow (not built integration):**
Symptomatic/adnexal workup → existing blood draw/order → reference-lab serum logistics (specimen unresolved until certified) → validated assay → versioned research pipeline → structured HL7/FHIR DiagnosticReport/Observation → clinician workup/referral under verified guidance. FHIR/CDS Hooks are design options, not deployed interfaces. Software/reference-lab path is a candidate; qRT-PCR point-of-care hardware is a longer-horizon analytical-validation/capital project, not a prototype proven here.

**Primary-source verification plan:** counsel/lab leads must date and authenticate current FDA IVD/LDT/PCCP/Q-Sub applicability, CMS CLIA, relevant court decisions and payer coverage/coding; clinical lead verifies SGO/ACOG referral guidance. IRB/privacy authority determines public/controlled/prospective status. No current legal, reimbursement or exemption assertion made.

**Funding milestones, no invented budget/partners:** (1) approval and primary provenance/access; (2) assay/data linkage and parity/guardrail validation; (3) prespecified single-platform replication with HE4/CA125 comparators; (4) untouched multisite prospective and calibration/utility evidence; (5) regulatory/lab feasibility review. Seek actual quotes and signed partner commitments before costs or names. Milestone completion is not asserted.
''')
# Execute only own audit utility, not any product code.
vcommand=['python','-X','utf8','audit/2026-02-08/integration/verify_checklist.py']
vout=subprocess.run(vcommand,capture_output=True,text=True,encoding='utf-8',errors='replace')
commands.append({'command':' '.join(vcommand),'exit_code':vout.returncode,'stdout':vout.stdout,'stderr':vout.stderr})
# Check command/path candidates in results index; commands remain unapproved.
verification={'audit_only':True,'authorization':False,'project_execution':False,'clinical_tests_passed_claim':False,'commands':commands,'integration_command':'python -X utf8 audit/2026-02-08/integration/consolidate.py','verifier_command':'python -X utf8 audit/2026-02-08/integration/verify_checklist.py','verifier_outcome':'87 expected/mapped; no errors; clinical_test=false','corrected_errors':'First checklist range check exposed nonexistent run_phase8_model_ensemble.py alias; corrected to actual run_phase8_model_experiments.py. Final source range checks have no failures (count recorded below). Phase7/9 command paths corrected to actual source names.','earlier_command_error':'Initial stdout cp1252 UnicodeEncodeError while printing reports; repeated with -X utf8. Full report log then read in bounded chunks. No project execution.','recommendations':{'inventory_count':len(recs),'mapped_count':len(rows),'unique_ids':len(set(r['recommendation_id'] for r in rows)),'each_once':len(rows)==len(recs)==len(set(r['recommendation_id'] for r in rows))},'citation_checks':{'count':len(checks),'failed':[r for r in checks if not r['exists_and_in_range']]},'product_hash_preservation':{'checked':len(pres),'mismatches':[r for r in pres if not r['unchanged']]},'endpoint_rows':len(endpoints),'coverage_evidence_rows':len(coverage),'grouped_saved_row_files':len(groups),'metric_evidence':'external-final/grouped_metrics.json; reports-review/TABLE_CHECKS.json and DISPLAY_CHECKS.json; research-reports/TABLE_CHECKS.json: inherited archived arithmetic, no fits; no pooled multi-experiment metric','limits':['Original input/fit metadata unknown','PDF and PNG rendered visual inspection unavailable','No clinical tests executed/passed','Static binary metadata only, no model deserialization']}
dump('VERIFICATION.json',verification)
# enumerate all audit artifacts including own final outputs; PDF explicitly owner added.
planned=['FILES_CHANGED.md','HASH_MANIFEST.json']
files=sorted(set(p.relative_to(ROOT).as_posix() for p in BASE.rglob('*') if p.is_file())|set('audit/2026-02-08/integration/'+n for n in planned))
put('FILES_CHANGED.md','# Audit artifacts created across the audit tree\n\nProduct source unmodified; preservation proof in VERIFICATION/PRESERVATION_CHECKS. Source PDF `Claude Methodology and Errors for Ovarian Cancer Project.pdf` was owner-added, NOT assistant changed. Existing historical outputs untouched. Non-integration artifact purposes are inferred from their explicit names, not claimed freshly authored by this worker. Earlier review/skeleton areas and external pooled checks NONAUTHORITATIVE SUPERSEDED, retained.\n\n'+ '\n'.join('- `'+p+'` — '+('NONAUTHORITATIVE SUPERSEDED preliminary artifact (retained)' if any('/'+a+'/' in p for a in ['governance','governance-review','inference-review','external','external-review','core-evaluation','discovery-review']) else ('offline reusable evidence utility' if p.endswith('.py') else 'coverage / endpoint / evidence / draft artifact: '+Path(p).name)) for p in files)+'\n')
dump('HASH_MANIFEST.json',{'algorithm':'SHA-256','self_reference_excluded':'HASH_MANIFEST.json','files':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='HASH_MANIFEST.json']})
print(json.dumps(verification,ensure_ascii=False,indent=2))
