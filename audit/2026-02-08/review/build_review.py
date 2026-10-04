"""Materialize explicitly reviewed ID dispositions; not a keyword classifier."""
from review_audit import *
# Each evidence anchor has a precise and limited meaning. Recommendation IDs are assigned explicitly below.
E={
'use':('scripts/ovarian_diagnostic_portal.py','render_result',201,219,'UI displays Stage I/II benchmarks and probability, not comparator or screening study evidence.','No intended-use agreement, comparator trial or clinical utility demonstrated.'),
'link':('src/data/multimodal_table.py','align_modalities',18,58,'ID-column and duplicate-key checks followed by clinical-universe reindexing; independent modalities are not concatenated into invented rows.','Matching strings do not verify donor, visit, specimen, label concordance or consent.'),
'endpoint':('scripts/test_high_capacity_models.py','load_training_data / run',56,144,'Five frozen names; early endpoint construction; identical outer folds compare frozen and two-product models; metrics pool OOF scores.','Reported N113=52/61 remains source/report evidence; no unbiased discovery provenance or independent clinical validation.'),
'external':('scripts/bulk_test_external_datasets.py','evaluate / test_one',293,374,'Threshold .5; saved bulk report sensitivity1 specificity0 means all-positive; source failure text is misleading.','Five-feature bulk model is not historical14 or nine-input portal; frozen failure mechanism is not proved platform causation.'),
'phase5':('scripts/evaluate_gse106817_phase5.py','run_direction',181,212,'Separate direction-specific 14-feature models and training-derived thresholds; saved raw .69 arm TN2788 FN320 is all-negative.','Within-target transformations in other arms are transductive; no nine-input portal validation.'),
'methyl':('scripts/evaluate_methylbert_methylation.py','load_data / evaluate',44,96,'Independent processed methylation evaluation and pooled OOF metrics.','Not a linked miRNA improvement; input absent and globally learned eligibility remains.'),
'stage':('scripts/train_tcga_ov_stage_model.py','load_features / main',89,196,'Barcode deduplication and tumor Stage III/IV versus I/II CV; eligibility learned globally.','Not cancer detection, blood screening, or Stage-I-only evidence.'),
'fusion':('scripts/smoke_test.py','main',14,43,'Random synthetic tensors exercise fusion forward; fixed missing-row assertion has a vacuous disjunct.','Source exists, not proof executed unit tests; no clinical training, all-absent safety or linked cohort.'),
'border':('scripts/bulk_test_external_datasets.py','label_from_metadata',137,150,'Metadata label rules exclude borderline and choose ovarian-positive/general negative terms.','Not universal separation: GSE94533 binary all-cancer endpoint includes borderline negatives; regex can misclassify no-cancer/other diseases.'),
'cal':('scripts/run_reproducible_baseline.py','make_pipeline',88,108,'Preprocessing and ANOVA precede classifier-only calibration.','Calibration implementation exists, but calibration-fold isolation and deployed probability validity remain unresolved; no DCA.'),
'discovery':('scripts/run_phase8_model_experiments.py','main',275,338,'Selection and transformed arrays precede inner model search; outer test is not used to fit those transforms.','Inner validation is contaminated; frozen high-capacity panel discovery is not independently reconstructed.'),
'filter':('scripts/evaluate_gse94533_diagnostic_model.py','main',62,96,'Full-frame missingness/constant eligibility precedes outer folds; searched pipeline handles later imputer/scaler/ANOVA.','Not fully fold-contained evaluation; no claim every model is nested.'),
'overlap':('src/data/multimodal_table.py','align_modalities',31,48,'Rejects duplicate modality keys and aligns by IDs.','This is not cross-study donor/visit overlap proof; accession/title separation alone is insufficient.'),
'error':('scripts/analyze_gse94533_errors.py','main',162,207,'Errors use saved decisions; subset metric routine uses .5; Stage-I CSV has cancer-only rows.','Seven saved models are not necessarily portal model; subtype FN counts need subtype denominators; regeneration was not run.'),
'manifest':('scripts/run_reproducible_baseline.py','run_model / main persistence',392,460,'Prediction rows persist model/fold/label/probability/threshold/decision and limited seed/version/path metadata.','No complete input hashes, code revision, sklearn versions or canonical run proof; all output claims remain source-attributed.'),
'access':('README.md','controlled-access roadmap',101,113,'README describes controlled access and separate assays.','No submitted access application, approval, HE4-linked data, assay permission or prospective protocol established.'),
'reg':('scripts/ovarian_diagnostic_portal.py','render_result',201,219,'Probability and suggested interpretation displayed in research UI.','No CLIA/FDA/payer determination, EHR interface, CDS guidance validation, workflow study or institutional approval evidenced.'),
'ui':('scripts/ovarian_diagnostic_portal.py','render_bulk_upload / render_result',157,219,'Bulk exports ID/probability/status; manual clinical defaults available; benchmarks hardcoded.','No per-row imputation flags/effective-input provenance or artifact-bound validated benchmark cards.'),
'companion':('audit/2026-02-08/pdf/extracted-pages.txt','physical PDF page15 companion introduction',539,551,'Inventory preserves recommendations and technical companion claims as unverified document assertions.','Text inventory is not PDF visual fidelity, clinical fact checking or authorization.'),
'fixed':('scripts/test_high_capacity_models.py','SIGNATURE / engineered_frame / run',31,144,'Five literal signature features and fixed products are evaluated conditional on chosen panel.','Outer CV does not independently validate the historical feature-choice process; untouched cohort proof missing.'),
'history':('scripts/reproduce_original_historical_clean.py','recover_historical_sample_order',470,671,'Conditional recovery requires119 numeric+60 E samples and179 unique IDs.','Does not recover original135/44 split; actual primary inputs absent; exported120/59 provenance contradicts algorithm.'),
'git':('audit/2026-02-08/governance-final/REPORT.md','Checkout and runtime',7,13,'Current HEAD and read-only git commands documented.','Current clone has git; original delivered snapshot context may differ; no commit/tag authorized.'),
'license':('audit/2026-02-08/governance-final/REPORT.md','Full-tree and PDF reconciliation',18,21,'Full-tree path inventory reports no LICENSE/COPYING, tests/CI, locks, DATA_SOURCES or RESULTS_INDEX.','Structural absence is not legal license interpretation; source reuse terms require per-source review.'),
'guards':('scripts/train_mirna_baseline.py','main',299,381,'Global filters; undefined patient_id_column and calibrate_threshold calls; test labels select threshold.','No successful baseline execution or guardrail suite established; review did not run source.'),
'paths':('scripts/ovarian_research_demo.py','module',4,36,'Research demo reads stored OOF results through Path.home()/ovarian-ai.','Not config-driven reproducibility or new inference; 45 hardcoded-home scripts are a structural count.'),
'lock':('requirements.txt','dependency list',1,18,'Dependencies have lower bounds rather than exact locked hashes.','Current interpreter is not original experiment environment; serialized compatibility not certified.'),
'ca':('scripts/ovarian_diagnostic_portal.py','sidebar_inputs',180,198,'Age/CA125 widgets default to stored medians; pending checkbox maps to NaN.','Median-valued unchecked inputs are not proof measured; source cohort missingness needs actual clinical table.'),
'logging':('scripts/ovarian_diagnostic_portal.py','main',222,252,'Session stores probability without immutable input snapshot, artifact hash or timestamp.','Not governed durable audit logging; stale probability can survive changed inputs or failed assessment.'),
'compare':('research/cohort_expansion/audit_2022_gockley_score_semantics.py','module and score semantics',1,80,'Published-network score-semantic audit distinguishes reconstruction from probability calibration.','Different circulating panels/assays and endpoints require exact feature/arm/score mapping; no head-to-head clinical superiority.'),
'prospective':('README.md','data separation roadmap',101,113,'Roadmap names future controlled/independent data sources.','No harmonized multi-site collection, repeated draws, same-visit comparator measurements or funding authorization.'),
}
# ID|anchor|status|recommendation-specific reviewed disposition / exact next action.
S='''01-01|use|needs_external_data_owner_decision|Agree clinical decision timepoint and comparator evidence before claiming earlier diagnosis.
02-01|use|needs_external_data_owner_decision|Owner must sign one intended-use tier and workflow-friction statement.
02-02|use|needs_external_data_owner_decision|Obtain comparator trial protocol or prospective symptom-cohort protocol for the selected tier.
03-01|use|needs_external_data_owner_decision|Draft Tier1 lead, funded Tier2 proposal, and explicit average-risk-screening exclusion.
05-01|link|implemented_verified|Retain independent cohorts; add donor/visit linkage adjudication before any fusion row construction.
05-02|endpoint|partial|Publish model-specific endpoint, subset, aggregation and source-attributed N beside every headline.
06-01|endpoint|partial|Pair full179 and narrow113 tables, declaring borderline handling and source-only counts.
06-02|endpoint|needs_external_data_owner_decision|Review interaction delta .0003 AUC and sensitivity/specification tradeoff; do not choose model automatically.
06-03|external|partial|Show confusion-derived external collapse direction, threshold and five-feature model identity next to AUC.
07-01|external|partial|Draft a fundable transfer-gap statement without claiming proven platform-shift causation.
07-02|methyl|implemented_verified|Keep methylation endpoint/table separate; verify independent study provenance before scientific claims.
07-03|stage|implemented_verified|Preserve tumor-staging label; prohibit use of TCGA metric as detection evidence.
07-04|fusion|partial|Label fusion synthetic scaffold; replace unit-tested claim with source-only smoke coverage until tests authorized.
07-05|discovery|partial|Qualify fixed seeds and row separation; remove blanket nested-CV rigor claim and describe inner contamination.
08-01|border|partial|Document accession screening and source-specific borderline rules, not a universal three-class guarantee.
08-02|endpoint|needs_external_data_owner_decision|Select a headline only after publishing complete endpoint registry and validation limitations.
08-03|cal|unresolved|Authorize whole-pipeline isolated calibration assessment; existing constructors do not validate portal risk.
08-04|cal|unresolved|Prespecify clinical utility/prevalence and compare decision curves with refer-all/refer-none after valid predictions.
08-05|discovery|unresolved|Design outer-fold discovery including candidate interactions and inner-contained tuning; compare conditional frozen panel.
08-06|filter|unresolved|Move eligibility fit inside each inner/outer training pipeline and test all-empty feature handling.
08-07|overlap|needs_external_data_owner_decision|Request donor/visit rosters and namespaces for each external cohort; title/GSM checks are insufficient.
08-08|error|partial|Bind subgroup outputs to model/run hash and regenerate only in a separately approved evaluation.
09-01|error|partial|Compute subtype denominators and sensitivity intervals from adjudicated histology; FN frequencies alone insufficient.
09-02|manifest|partial|Add command, code/input/artifact hashes, software versions and output lineage to each numerical run.
09-03|access|needs_external_data_owner_decision|Verify phs000897 eligibility and terms; owner authorizes application; do not state submitted.
09-04|access|needs_external_data_owner_decision|Acquire consented same-visit HE4 comparator cohort only after source/ethics review.
09-05|external|needs_external_data_owner_decision|Request original assay/run/batch metadata; separate target-adapted analyses from frozen transfer.
09-06|link|partial|Validate same donor and visit, specimen and label agreement before multimodal alignment; retain independent tables otherwise.
09-07|prospective|needs_external_data_owner_decision|First obtain repeated-draw timestamps and clinical outcomes; do not fabricate longitudinal trends.
09-08|prospective|needs_external_data_owner_decision|Owner/lab/clinical partners prespecify prospective single-platform Tier1 protocol and statistical plan.
10-01|reg|needs_external_data_owner_decision|Obtain current legal/regulatory advice on LDT/IVD and update/PCCP strategy.
10-02|reg|needs_external_data_owner_decision|Align indication and reimbursement with evidence-backed workup rather than screening promises.
10-03|reg|needs_external_data_owner_decision|Seek consultant/Q-Sub advice before locking prospective evidence collection.
11-01|reg|needs_external_data_owner_decision|Validate reflex/parallel blood-order workflow with clinical partners; no new-visit savings claim yet.
11-02|reg|unresolved|Specify governed HL7/FHIR results interface and test fixtures; portal is not EHR integration.
11-03|reg|needs_external_data_owner_decision|Clinicians review current SGO/ACOG referral guidance before CDS/labels are designed.
11-04|reg|needs_external_data_owner_decision|Validate specimen handling, stability and reference-lab logistics analytically.
11-05|reg|needs_external_data_owner_decision|Frame software/reference-lab and qRT-PCR hardware as future alternatives requiring validation/capital.
11-06|endpoint|unresolved|Draft one-page scorecard with N, endpoint, aggregation, external model-specific failure and provenance limits.
12-01|external|unresolved|Add costed transfer-gap closure milestones to scorecard; not a funded executed plan.
12-02|ui|unresolved|Block clinically framed demo until bulk imputation/effective-input flags and endpoint/N limitations are exposed.
12-03|reg|unresolved|Draft clearly proposed regulatory/reimbursement and clinical workflow diagrams for expert review.
12-04|access|needs_external_data_owner_decision|Owner approves funding ask and partners for access, HE4, replication and CLIA pathway; no acquisition here.
12-05|external|unresolved|Rehearse all-positive bulk versus all-negative selected phase5 distinction without asserted biological cause.
13-01|endpoint|needs_external_data_owner_decision|Owner chooses headline conditional estimate and shows every competing endpoint/model table.
13-02|cal|unresolved|Develop fold-safe calibration and prespecified decision-curve plan, not a plot from unvalidated portal scores.
14-01|ui|unresolved|Attach endpoint/N and per-row clinical imputation flags to manual/bulk output.
14-02|fusion|partial|Label every fusion mention untrained synthetic scaffold; remove verified unit-tested wording.
14-03|reg|unresolved|Prepare proposal-only workflow/regulatory/reimbursement diagrams with unresolved approvals visible.
14-04|access|needs_external_data_owner_decision|Verify access application requirements; owner may apply; claim in-progress only with dated evidence.
14-05|companion|implemented_verified|Companion recommendation inventory consulted; keep allegations conditional and reconcile source-specific evidence.
15-01|endpoint|partial|Place full179/narrow113 composition and conditional panel provenance beside every percentage.
16-01|fixed|partial|Describe frozen-panel CV conditional on prior panel choice, not unbiased new-patient discovery performance.
16-02|discovery|unresolved|Put interaction-space selection and tuning inside nested folds; canonical estimate requires approved protocol.
17-01|fixed|needs_external_data_owner_decision|Obtain genuinely untouched donor-disjoint, assay-compatible cohort with locked pipeline before evaluation.
17-02|endpoint|needs_external_data_owner_decision|Owner/clinician reviews missing-data parity and sensitivity-specificity objective before interaction model choice.
17-03|filter|unresolved|Fold-fit missingness/variance eligibility as sklearn transformer; preserve old artifacts unchanged.
17-04|external|partial|Disclose bulk all-positive collapse and selected phase5 all-negative collapse; imbalance accuracy is not utility.
18-01|history|unresolved|Exclude historical independent-test claim until original135/44 roster and frozen score semantics verified.
18-02|stage|implemented_verified|Keep unsuccessful tumor-staging work distinct from cancer-detection evidence.
18-03|git|obsolete|Do not initialize git; owner may later approve version/tag policy for future runs, no commit authorized now.
18-04|license|needs_external_data_owner_decision|Owner/legal chooses explicit code license after dependency/source obligations review.
18-05|guards|unresolved|Implement synthetic guardrails only after approval, starting undefined helpers and data/fold/schema assertions.
19-01|overlap|unresolved|Add donor-disjoint folds and cross-cohort tests plus full pipeline fit containment and artifact dimension parity tests.
19-02|paths|partial|Route every reproducer through resolved configuration paths, removing home-relative checkout assumptions.
19-03|lock|unresolved|Snapshot exact dependencies and persist sklearn/numpy/joblib versions and hashes per versioned artifact.
19-04|manifest|unresolved|Owner selects canonical result scripts; add RESULTS_INDEX lineage and mark variants without moving originals yet.
20-01|access|needs_external_data_owner_decision|Prioritize consented HE4-linked comparator data acquisition subject to protocol and permissions.
20-02|error|partial|Report subtype sensitivity with denominators/uncertainty per model; generalization remains open.
20-03|stage|needs_external_data_owner_decision|Obtain adjudicated Stage-I-only metadata; never infer I from combined I/II or all-malignant transfer.
20-04|ca|partial|Expose measured versus imputed CA125 and verify source missingness counts from original clinical data.
20-05|fusion|partial|Separate scaffold source and synthetic smoke from linked clinical training evidence.
21-01|fusion|unsupported|Replace asserted unit-tested wording with synthetic smoke source present, execution/guardrail coverage unestablished.
21-02|ui|partial|Attach source-reported113=52/61 and external/unbiased-discovery limitation directly to benchmark cards.
21-03|use|needs_external_data_owner_decision|Prespecify clinical operating-point utility and validate tradeoffs; no threshold change approved.
21-04|cal|unresolved|Validate isolated whole-pipeline probability calibration; max(p,1-p) is not clinical confidence.
21-05|ui|unresolved|Return effective Age/CA125, missing/default flags and input snapshot in every bulk row.
21-06|logging|unresolved|Design privacy-governed input/artifact-bound audit log with retention/access policy before implementation.
21-07|access|needs_external_data_owner_decision|Institutional review determines human-subjects/IRB status; public data alone is not exemption proof.
22-01|access|needs_external_data_owner_decision|Write public-data human-subjects statement only after institutional/source review and documented determination.
22-02|license|unresolved|Create source-by-source citation/reuse/consent/access inventory with authoritative terms, no blanket public-license assumption.
22-03|license|needs_external_data_owner_decision|Owner/legal approves code license; do not infer permission from repository visibility.
22-04|manifest|partial|Track each repeated priority separately: endpoint/UI/scaffold/provenance/index/folds/calibration/access/model tradeoff/license; none authorizes execution.
23-01|compare|needs_external_data_owner_decision|Obtain published panel weights, mature arms, preprocessing and comparable endpoint before head-to-head claims.
23-02|cal|unresolved|Calibrate complete fold-isolated pipeline on development only, evaluate untouched probabilities before literal clinical risk wording.
23-03|external|partial|Foreground model-specific observed transfer collapse, separating AUC/decisions from causal interpretation.
23-04|prospective|needs_external_data_owner_decision|Owner/partners fund prespecified harmonized multisite prospective collection rather than implied post-hoc success.'''
inv=(BASE/'pdf/recommendation-inventory.md').read_text(encoding='utf-8')
recommendations=dict(re.findall(r'- (PDF-P\d+-R\d+): ([^\n]+)',inv))
checklist=[]
for entry in S.splitlines():
    short,key,status,action=entry.split('|',3);pg,n=short.split('-');rid=f'PDF-P{pg}-R{n}';p,fn,a,b,e,l=E[key]
    # Exact snippets are read from authoritative bytes, never borrowed unrelated keyword matches.
    ls=(ROOT/p).read_text(encoding='utf-8').splitlines(); b=min(b,len(ls))
    checklist.append({'recommendation_id':rid,'physical_page':int(pg),'recommendation':recommendations[rid],'status':status,'source_path':p,'source_function':fn,'inclusive_lines':f'{a}-{b}','evidence':e,'limits':l,'next_action':action,'review_basis':'Explicit recommendation-specific reasoning; source excerpt metadata verified offline; no runtime validation.'})
assert len(checklist)==87 and set(recommendations)=={r['recommendation_id'] for r in checklist}
writecsv('PDF_CHECKLIST.csv',checklist)
# Preserve ALL heterogeneous inventory aliases, but normalize experiment families from curated subsystem endpoints only.
raw=rows(BASE/'integration/ENDPOINT_REGISTRY.csv')
writecsv('ENDPOINT_ALIASES.csv',[{'alias_id':r['registry_id'],'lineage_file':r['evidence_file'],'lineage_row':r['evidence_row'],'evidence_kind':'heterogeneous source/report/metadata evidence; NOT an experiment count','raw_evidence_json':json.dumps(r,ensure_ascii=False)} for r in raw])
# One endpoint row per pipeline/endpoint; historical positional counters and PDF summary aliases are explicitly not experiments.
selected=[r for r in raw[:105] if int(r['registry_id'].split('-')[1]) not in set(range(35,36))|set(range(39,57))|{76,92,93,94,105}]
registry=[]
for r in selected:
    num=int(r['registry_id'].split('-')[1]); src=r.get('source') or r.get('variant') or r.get('evidence') or r['evidence_file']
    endpoint=r.get('endpoint') or r.get('labels_exclusions') or r.get('labels_exclusions') or r.get('identity_endpoint') or 'See retained source row for exact endpoint coding'
    registry.append({'family_id':f'F-{num:03}','pipeline':src,'endpoint':endpoint,'specimen_platform':r.get('specimen') or 'Unresolved: consult lineage; do not assume serum/platform equivalence','class_coding_borderline':r.get('labels_exclusions') or r.get('exclusions') or r.get('excluded') or 'Cancer usually1; lineage must adjudicate negatives/borderline; not universal','source_frozen_vs_discovery':r.get('features') or r.get('function') or 'Distinct pipeline; fixed panel and target development are not interchangeable','stage_screening':r.get('eligibility') or 'Only explicitly named stage endpoint; I/II is not I; no average-risk screening validation','N_source_reported':r.get('evaluated_n') or r.get('claimed_n') or r.get('asserted_counts') or r.get('reported_or_expected_count') or r.get('expected_reported_counts') or r.get('source_expected_counts') or 'Unresolved; retained lineage contains claims, not adjudicated patients','N_actual_evidence':r.get('verified_local_counts') or r.get('observed_source_data_counts') or r.get('local_established_count') or 'Not independently clinically counted in review; saved-row/metadata counts only if explicitly identified in lineage','measured_clinical_availability':r.get('clinical_measurement_counts') or 'Unresolved measured Age/CA125/HE4; defaults/imputation are not measurement','folds_threshold_model_version':r.get('protocol') or r.get('evaluation') or r.get('aggregation') or 'Source-specific protocol in lineage; no borrowed threshold/model/version','exclusions':r.get('exclusions') or r.get('excluded') or r.get('limitations') or r.get('limitation') or 'Unresolved; see source-specific exclusions','evidence_aliases':r['registry_id'],'lineage':f"{r['evidence_file']} row {r['evidence_row']}",'limits':r.get('limitation') or r.get('limits') or r.get('limitations') or 'Source-reported static evidence; no clinical/model rerun; not proof independent patients'})
# Additional explicit model variants that must never collapse into the portal.
variants=[
('bulk-five','GSE106817/GSE113486/GSE113740','five selected source features; balanced logistic; .5','serum array; platform-specific signals','all ovarian vs regex noncancer; borderline excluded; not stage-filtered','3108/140/1058 saved evaluated; GSE113486972 series total','EP-0001;EP-0002;EP-0003;EP-0009'),
('phase5-A-raw','GSE94533 to GSE106817','historical14 logistic; raw; threshold .69; source CV','TPM source vs array target','all-cancer transfer; no Stage-I validation','3108 saved labels320/2788; TN2788 FP0 FN320 TP0','EP-0001'),
('phase5-A-standard','GSE94533 to GSE106817','historical14 logistic; source standard; .6','TPM source vs array target','distinct preprocessing arm','3108 saved; TN2787 FP1 FN320 TP0','EP-0001'),
('phase5-A-quantile','GSE94533 to GSE106817','historical14 logistic; source quantile; .48','TPM source vs array target','distinct preprocessing arm','3108 saved; TN2788 FP0 FN320 TP0','EP-0001'),
('phase5-B-quantile','GSE106817 to GSE94533','historical14 logistic; source quantile; .48','array source vs TPM target','reverse-direction distinct training cohort','179 saved; TN0 FP83 FN0 TP96','EP-0001'),
('external-selected','GSE106817','separate SelectKBest/C/class-weight search; source OOF threshold','serum assay transfer','all-cancer; source discovery, not frozen14 or portal','3108 report; raw source absent','EP-0001'),
('portal-nine','GSE94533 source only','five frozen miRNAs+Age+CA125+two products; .5; manifest version unresolved','source TPM-like; clinical units require metadata','I/II52 vs benign/control61; borderline excluded','113 source-reported; no external clinical validation','EP-0018;EP-0025'),
('frozen-seven','GSE94533','five frozen miRNAs+Age+CA125; logistic .5; conditional outer5','source miRNA plus clinical','I/II52 vs benign/control61','113 source-reported','EP-0025'),
('engineered-nine-evaluator','GSE94533','products BEFORE imputation; pooled5 OOF; .5','source miRNA plus clinical','I/II vs benign/control','113 source-reported; different missingness transform from portal','EP-0025'),
('cli-seven','selected artifact directory','seven columns no manifest; .5; incompatible default nine scaler','research inference interface, not experiment','no validated clinical endpoint','N not applicable','EP-0018'),
('fusion-synthetic','random synthetic tensors','gated network synthetic smoke only','synthetic multimodal','no clinical classes or screening','no patient N','EP-0069'),
('GSE211-strict196','GSE211692','candidate five-name assay compatibility; not portal validation','serum GPL21263 signed normalized array','I/II168 vs ovarian benign28; controls separate','196 metadata/report; donor independence unresolved','EP-0005;EP-0101'),
('GSE211-proposed368','GSE211692','proposed selected-control endpoint','serum GPL21263','I/II168 vs selected200; public roster unknown','368 proposed, not reconstructed','EP-0006'),
('GSE211-regex1026','GSE211692','bulk regex; not strict endpoint','serum GPL21263','400 cancer/626 regex negatives; misses no cancer, pools benign diseases','1026 metadata replay; series16190','EP-0004'),
('FirePlex-full','2024 adnexal','published10 panel only9 exact features available; weights absent','serum FirePlex reported; plasma claim contradiction unresolved','all malignant vs pooled controls, not OVAI validation','678 workbook286/392; train468 internal100 external110','EP-0007;EP-0098'),
('FirePlex-external','2024 adnexal','published split, not recreated model','serum FirePlex; specimen provenance unresolved','external malignant51 vs control59;20 early cancer source reported','110 workbook, no model predictions reproduced','EP-0008'),
('FirePlex-overlap8','2024 adnexal','eight-feature descriptive historical overlap;6/8 directions agree','FirePlex versus historical assay','not frozen14 nor opposite-arm substitution','678 workbook source claims','EP-0098'),
('PMML-published','Gockley','14 to7 tanh to2 tanh; argmax, NOT .5 probability','original frozen affine normalized assay','historical binary cancer vs benign/borderline/control','179 processed; original135/44 unproved','EP-0034;EP-0036;EP-0037;EP-0038'),
]
for i,(name,cohort,model,platform,end,n,aliases) in enumerate(variants,1):
    if name.startswith('phase5-'): continue # Expanded below from all twelve distinct saved grouped arms.
    registry.append({'family_id':f'V-{i:03}','pipeline':name+' / '+cohort,'endpoint':end,'specimen_platform':platform,'class_coding_borderline':end,'source_frozen_vs_discovery':model,'stage_screening':'Only explicitly named stages; no screening utility; I/II not I','N_source_reported':n,'N_actual_evidence':'Saved rows/metadata/workbook only, not adjudicated clinical patients','measured_clinical_availability':'Absent or unresolved; no HE4-linked comparator; imputation not measurement','folds_threshold_model_version':model,'exclusions':'Endpoint-specific; same donor/visit independence and exclusions require source adjudication','evidence_aliases':aliases,'lineage':'ENDPOINT_ALIASES.csv plus finalized subsystem REPORT/ENDPOINTS and external grouped_metrics','limits':'Separate normalized variant; lineage alias can inform multiple pipelines. Not a count of independently executed experiments.'})
# Saved grouped-arm metadata expansion: no probability calculations or clinical/model reruns.
for ai,audit in enumerate(json.loads((BASE/'external-final/grouped_metrics.json').read_text(encoding='utf-8'))['audits'],1):
    for gi,group in enumerate(audit['groups'],1):
        info=group['group']; phase=ai==1
        registry.append({'family_id':f'G-{ai:02}-{gi:02}','pipeline':('historical14 logistic phase5' if phase else 'saved OOF model arm')+' / '+json.dumps(info,sort_keys=True),'endpoint':('direction-specific all cancer versus source-specific noncancer' if phase else 'GSE94533 early I/II vs benign/control61' if ai in (3,4) else 'GSE94533 all cancer96 vs noncancer83 including borderline'),'specimen_platform':'Source/target assay-specific TPM versus array' if phase else 'GSE94533 serum miRNA source-reported; measured clinical lineage unresolved','class_coding_borderline':'Binary1 cancer; phase5 source-specific negatives; error analysis includes borderline; augmentation/ensemble exclude borderline','source_frozen_vs_discovery':'Historical fixed14 direction-specific source-trained; target-adapted modes separately labeled' if phase else 'Saved OOF model arm, not portal artifact or independent discovery replication','stage_screening':'No Stage-I-only or average-risk-screening claim','N_source_reported':str(group['n']),'N_actual_evidence':f"{group['n']} saved prediction rows in this arm; TN={group.get('tn','unknown')} FP={group.get('fp','unknown')} FN={group.get('fn','unknown')} TP={group.get('tp','unknown')}; inherited grouped arithmetic, not new clinical count",'measured_clinical_availability':'Saved covariates/defaults are not raw-source measured availability; HE4 linkage not established','folds_threshold_model_version':json.dumps(info,sort_keys=True)+'; source script version/hash retained in coverage; no recovered complete runtime lock','exclusions':'Source-specific labels; repeated models/modes reuse patients, not new cohorts','evidence_aliases':'EP-0001' if phase else 'EP-0024' if ai==3 else 'EP-0023' if ai==4 else 'EP-0030','lineage':f"audit/2026-02-08/external-final/grouped_metrics.json audits[{ai-1}].groups[{gi-1}]; {audit['path']}",'limits':'Inherited saved metadata; probability-only arms without decision columns do not establish original operating points.'})
# Enrich core curated families with explicit source fields, preserving unknowns without inference.
for item in registry:
    if item['family_id'].startswith('F-'):
        rawrow=next(r for r in raw if r['registry_id'] in item['evidence_aliases'].split(';'))
        item['positive_code_definition']=rawrow.get('positive_code') or rawrow.get('positive') or 'Unresolved; exact retained source row'
        item['negative_code_definition']=rawrow.get('negative_code') or rawrow.get('negative') or 'Unresolved; exact retained source row'
        item['source_evidence_exact']=rawrow.get('evidence') or rawrow.get('saved_evidence') or rawrow.get('source') or 'See alias source'
        item['source_fields_json']=json.dumps(rawrow,ensure_ascii=False)
        n=int(item['family_id'].split('-')[1])
        if 57<=n<=73:
            item['class_coding_borderline']=f"positive={rawrow.get('positive')}; negative={rawrow.get('negative')}; exclusions={rawrow.get('exclusions')}"
            item['stage_screening']=rawrow.get('eligibility')+'; no screening validation'
            item['specimen_platform']='GSE106817 serum microarray' if n in (70,73) else 'GSE94533 serum miRNA TPM source-reported; primary raw assay metadata absent'
        elif 19<=n<=33:
            item['specimen_platform']='GSE94533 serum miRNA TPM plus clinical where specified; external rows assay-specific'
        elif 77<=n<=83:
            item['specimen_platform']='Tumor tissue microarray; not serum detection'
        elif 87<=n<=91:
            item['specimen_platform']='TCGA tumor molecular assays; paired where source says paired; not screening'
    else:
        item['positive_code_definition']='Cancer/malignant positive where binary; synthetic/CLI not applicable; PMML cancer activation oriented, argmax rule'
        item['negative_code_definition']=item['class_coding_borderline']
        item['source_evidence_exact']=item['lineage']
        item['source_fields_json']='Variant-specific report/source references via retained aliases; no clinical availability inferred'
for item in registry:
    item['record_type']='pipeline_endpoint_family' if item['family_id'].startswith('F-') else 'distinct_model_or_saved_arm_variant'
    item['parent_boundary']='Family and child variants are hierarchical, not additive experiments; alias evidence may support multiple related records.'
writecsv('ENDPOINT_REGISTRY.csv',registry)
