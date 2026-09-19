# OVAI Current Pipeline Audit

Audit date: 2026-09-19

Scope: the 83 scripts under `scripts/`, committed data/model/report artifacts, top-level documentation, and the current metric-producing pipelines. No model code or data was changed for this audit.

## Executive Findings

1. The canonical GSE94533 diagnostic result is a patient-level, nested 5-fold OOF evaluation on 179 patients: 96 cancer and 83 non-cancer. The 14-miRNA plus age plus CA-125 logistic model produces 73.74% accuracy and ROC-AUC 0.8131. This is the most coherent current serum diagnostic benchmark, but it is internal resampling within one cohort, not an independent external test.
2. The reported 76.54% accuracy and 77.06% balanced accuracy come from the calibrated-SVM branch of `evaluate_gse94533_enhancements.py`. Its threshold is selected separately inside each outer training fold, so the reported threshold-tuned OOF estimate is not an ordinary fixed-0.5 estimate. It is a valid nested research estimate, but it is still single-cohort OOF performance and should not be called external validation.
3. GSE106817 is an independent serum cohort for ranking assessment, but its current frozen evaluation predicts every sample as non-cancer at the training-derived threshold: accuracy 89.70%, balanced accuracy 0.50, sensitivity 0%, specificity 100%. Its ROC-AUC 0.8232 is the informative result. The accuracy is a class-imbalance artifact and must not be reported as diagnostic success.
4. The methylation result of 88.14% accuracy, 86.90% balanced accuracy, and 0.9414 ROC-AUC is nested OOF evaluation within the MethylBERT supplementary cohort, not an independent test of a model trained elsewhere. The early-stage comparison is 84.13% accuracy, 77.26% balanced accuracy, and 0.8609 ROC-AUC. These data are not patient-linked to GSE94533.
5. The repository does not contain evidence that a claimed historical `>80%` independent-test result was reproduced. The exact historical PMML model is evaluated on all 179 historical patients, which is apparent historical performance. A separate script searches 20,000 compatible 44-patient splits near the published AUC and explicitly says a candidate is not proof of the original split. No committed output identifies the original 44 patients or a verified independent test metric. The current `>80%` numbers that can be traced exactly are the methylation within-cohort OOF accuracy and the invalid GSE106817 all-negative accuracy.
6. There is no valid patient-level multimodal cancer-detection fusion model in the current committed results. GSE94533 miRNA plus clinical values are patient-linked. TCGA RNA/miRNA/clinical features are patient-linked but use a tumor-stage endpoint. Methylation and serum miRNA cohorts remain separate, as required by the repository's scientific rule.

## Canonical Current Pipelines

### GSE94533 serum miRNA diagnosis

Primary scripts:

- Dataset construction: `build_gse94533_clean_dataset.py`, `build_gse94533_clinical_fusion.py`
- Evaluation: `evaluate_gse94533_diagnostic_model.py`
- Artifact: `reports/tables/gse94533_diagnostic_results.json`

Data and endpoint:

- Current GEO annotations define Cancer versus Benign, Borderline, or Control.
- 179 unique patient/sample IDs are aligned between metadata and the cleaned miRNA matrix; 96 cancer and 83 non-cancer.
- One sample, E-1056, is excluded from the current GEO-defined cohort. This differs from the historical eLife reconstruction, which has its own sample-order and publication-label logic.
- The feature matrix contains 2,582 miRNAs before model-layer selection. The clinical fusion table contains the fixed 14 historical miRNAs, age, and CA-125 where available.

Preprocessing and modeling:

- `log1p` is applied to nonnegative miRNA TPM values before cross-validation. This transform has no label dependence, but is performed once across the full matrix rather than represented as a fold-fitted pipeline step.
- The evaluator removes miRNA columns using global non-missing and non-constant checks before splitting. This is unsupervised test-distribution information. It is a small but real transductive dependency and should be moved inside the training fold for a strict external-generalization protocol.
- Imputation, standardization, univariate `SelectKBest`, and logistic-regression fitting are inside a scikit-learn `Pipeline` and are refit in each inner training fold.
- Hyperparameters (`k`, `C`, and class weighting) are selected by inner 3-fold CV. Outer evaluation is shuffled 5-fold stratified OOF with seeds 42 and 43.
- The 0.5 prediction threshold is fixed for the baseline models.

Reported results:

| Layer | Accuracy | Balanced accuracy | ROC-AUC | Interpretation |
|---|---:|---:|---:|---|
| 14 miRNA | 70.39% | 70.68% | 0.7790 | Nested single-cohort OOF |
| All usable miRNA | 69.27% | 69.40% | 0.7687 | Nested single-cohort OOF; high-dimensional layer does not improve |
| 14 miRNA + age + CA-125 | 73.74% | 73.97% | 0.8131 | Nested single-cohort OOF; current serum reference |

The source values are in `reports/tables/gse94533_diagnostic_results.json`. The published current README wording correctly describes these as research metrics, not clinical validation.

### GSE94533 enhancement benchmark

Script: `evaluate_gse94533_enhancements.py`  
Artifact: `reports/tables/gse94533_enhancement_results.json`

The same 179-patient clinical fusion cohort is evaluated with median or KNN imputation, logistic regression, calibrated linear SVM, random forest, XGBoost, ratio features, and a stack. Outer and inner CV use fixed seeds 42 and 43. All candidate preprocessing is in pipelines and model selection is nested within outer folds.

Important interpretation details:

- The calibrated SVM's probability calibration uses an internal 3-fold CV inside the outer training data.
- Its classification threshold is selected from inner-training OOF probabilities separately for each outer fold. The reported `threshold_tuned_metrics` aggregate fold-specific predictions and are therefore not equivalent to choosing 0.58 once on all OOF predictions.
- The threshold is selected to maximize balanced accuracy, and a separate specificity-target threshold is also selected. These are research operating-point analyses, not locked clinical thresholds.
- The lasso stability analysis is fit only on each outer training fold, but it is descriptive and does not define the reported model.
- The stack has multiple flexible learners and was evaluated under nesting, but the observed ranking result does not beat the simpler fusion baseline.

The headline calibrated-SVM threshold-tuned result is 76.54% accuracy, 77.06% balanced accuracy, and 0.8084 ROC-AUC. It remains an internal OOF estimate.

### MethylBERT supplementary methylation benchmark

Scripts: `evaluate_methylbert_methylation.py`, `freeze_methylbert_pipeline.py`  
Artifacts: `reports/tables/methylbert_methylation_results.json`, `models/methylbert_ridge_cancer_vs_healthy.joblib`

The processed HRA011451/PMC11384945 supplementary matrix is evaluated separately from GSE94533. Cancer versus healthy includes 1,872 rows (754 cancer, 1,118 healthy); early-stage versus healthy includes 1,323 rows (205 early stage, 1,118 healthy).

- Imputation, scaling, feature selection, Ridge classification, and hyperparameter selection are inside nested 5-fold/3-fold pipelines.
- Seeds are fixed at 42 and 43.
- Feature availability is filtered globally using non-missing and non-constant checks before the outer split. This uses unlabeled information from held-out rows and should be moved into the fold for a strict protocol.
- The saved joblib is refit on the full cancer-versus-healthy training cohort after model selection. It is a frozen training artifact, not an external validation result.
- The proposed external targets are invalid or unavailable: GSE144365 is not an ovarian methylation target, TCGA has incompatible `cg...` tissue probes, and no feature-aligned HRA011451 external matrix is available.

The 88.14% cancer-versus-healthy accuracy is therefore not evidence that a frozen methylation model has achieved 88% on independent patients.

### TCGA-OV stage research endpoint

Script: `train_tcga_ov_stage_model.py`  
Artifact: `reports/tables/tcga_ov_stage_results.json`

This is not a cancer-detection endpoint. RNA, miRNA, and clinical data are joined by TCGA patient barcode, repeated aliquots are reduced to one sample per patient, and the target is stage III/IV versus stage I/II. Nested 5-fold/3-fold CV contains imputation, scaling, feature selection, and tuning; seeds are fixed at 42 and 43.

The stage cohort is highly imbalanced (22 stage I/II versus 277 stage III/IV). The miRNA-only accuracy of 78.93% is below the 92.64% majority-class accuracy and has balanced accuracy 53.06% and ROC-AUC 0.4889. It is not a successful cancer-detection result.

## Independent-Cohort Checks

### GSE106817

Scripts: `audit_gse106817.py`, `extract_gse106817_validation_cohort.py`, `evaluate_gse106817_external.py`  
Artifact: `reports/tables/gse106817_external_results.json`

- The model is trained on GSE94533 and evaluated on 3,108 GSE106817 samples (320 cancer, 2,788 non-cancer) using mapped platform features.
- GSE106817 is not used for feature selection, normalization fitting, hyperparameter selection, or ComBat. Standard scaling and quantile transformation are fit on training data only.
- The training-side model search selects feature count, regularization, and class weighting. That is allowed for training, but the external script then selects the frozen threshold from cross-validated predictions generated using the already-selected estimator. The threshold is training-only but reuses the search-selected estimator's data for operating-point optimization; a fully nested threshold-selection procedure would make this cleaner.
- The external labels are parsed from sample titles. The script does not perform a cross-accession/patient identity overlap audit against GSE94533.
- The standard-scaling result has ROC-AUC 0.8232 but predicts no positives at the frozen threshold. Accuracy 89.70% is entirely driven by the 89.70% non-cancer prevalence.
- The all-miRNA GSE106817 and cancer-versus-borderline scripts are discovery/within-cohort experiments, not frozen external validation. They use nested CV for their own cohort but must not be presented as independent tests of GSE94533.

### Historical/frozen gene signatures on GSE12172

Scripts: `evaluate_frozen_signature_threshold.py`, `test_frozen_signature_stage1_gse12172*.py`, `test_gse8841_signature_on_gse12172*.py`

These scripts train or freeze a GSE8841-derived signature, fit discovery scaling on GSE8841, and apply the resulting model to GSE12172. The threshold in `evaluate_frozen_signature_threshold.py` is selected from GSE8841 OOF predictions only. The intended direction is scientifically sound, but the committed workspace does not include the referenced GSE8841/GSE12172 data or result artifacts, so these claims cannot be independently verified here. Several near-duplicate final/clean/standalone scripts exist.

## Historical eLife/PMML Provenance

Relevant scripts:

- Exact model/data reconstruction: `reproduce_original_historical_clean.py`, `reproduce_original_historical_standalone.py`, `reproduce_original_pmml.py`, `reproduce_original_pmml_historical_data.py`
- PMML inspection and comparison: `audit_original_pmml_model.py`, `extract_original_14mirna_model.py`, `compare_pmml_to_published_table2.py`, `validate_pmml_against_table2_by_id.py`
- Split investigation: `search_original_split_evidence.py`, `search_split_evidence_clean.py`, `search_possible_historical_test_splits.py`, `search_original_44_patient_split.py`

What is established:

- The historical architecture is a fixed 14 -> 7 tanh -> 2 tanh PMML network. No training occurs during reproduction.
- The clean reconstruction uses the historical supplementary matrix and publication ordering to reconstruct 179 patients and preserves the publication's excluded numeric position 56 logic.
- The exact PMML output is evaluated on the historical 179-patient cohort. That is apparent performance on the data represented by the historical model and cannot be called an independent test.
- The publication describes a 135/44 training/test design with pathology strata. The repository searches source materials for the actual identities but does not recover them.
- `search_possible_historical_test_splits.py` samples 20,000 stratified 44-patient candidates using a fixed RNG seed and ranks candidates near a published AUC of 0.90. Its own output text says a candidate is not proof of the original split. Selecting a candidate because it matches the published metric is a multiple-comparisons exercise and cannot establish independent validation.
- No committed historical result JSON/CSV is present in this checkout that identifies the original 44-patient test IDs and reports a verified independent metric.

Conclusion: the historical `>80%` independent-test claim is unresolved in this repository. It must not be repeated as a validated result. The traceable `>80%` values currently present are the methylation within-cohort OOF accuracy and the GSE106817 all-negative accuracy; neither supports the claimed conclusion.

## Leakage, Overfitting, and Validity Register

| Severity | Finding | Affected surface | Consequence |
|---|---|---|---|
| Critical | GSE106817 accuracy is reported despite 0% sensitivity and balanced accuracy 0.5 | `gse106817_external_results.json` | Accuracy is an invalid headline under severe prevalence shift |
| Critical | Historical PMML full-cohort performance and compatible-split searches can be mistaken for independent testing | PMML and split-search scripts | Apparent or selected-to-match performance is not external validation |
| High | No committed proof of the original historical 44-patient IDs or test result | Historical pipeline | The claimed `>80%` independent result cannot be reproduced or audited |
| High | External GSE106817 label/accession overlap with GSE94533 is not checked | `evaluate_gse106817_external.py` | Patient independence is assumed, not demonstrated |
| Moderate | Global feature availability filters use held-out rows in GSE94533 and methylation evaluation | Current diagnostic and methylation evaluators | Small transductive leakage; must be fold-contained for strict generalization |
| Moderate | GSE106817 frozen threshold is selected after training model search using reused training-side CV predictions | External evaluation | Operating point may be optimistic on the training cohort, even though the target cohort remains untouched |
| Moderate | Some data-build scripts median-impute the complete dataset before any downstream model can use it | `build_tcga_clean_multimodal_dataset.py` and legacy builders | Safe only if never used as model input; otherwise leakage across folds |
| Moderate | Several 192-miRNA and historical baselines use fixed feature lists and a single 5-fold OOF pass without nested hyperparameter selection | `build_age_192mirna_model.py`, `build_age_ca125_192mirna_model.py`, related legacy scripts | Results are exploratory and less protected against model-choice overfitting |
| Moderate | Many legacy scripts use `Path.home()/ovarian-ai` or relative `data` paths | Historical, GSE106817, TCGA-build, and frozen-signature scripts | The current `/workspaces/OVAI` checkout is not reproducible by those commands without path changes |
| Low | Seeds are fixed in most canonical evaluators, but not uniformly across all 83 scripts and dependency-level determinism is not recorded | Repository-wide | Repeated runs may differ in legacy/discovery pipelines |
| Low | Class imbalance is not always paired with balanced accuracy or a prevalence-aware baseline in older scripts | Historical and GSE106817 baselines | Accuracy can be materially misleading |

## Script Inventory and Pipeline Generations

The 83 scripts fall into these functional groups:

- Current GSE94533 ingestion/audit/evaluation: `investigate_gse94533_labels.py`, `build_gse94533_clean_dataset.py`, `build_gse94533_clinical_fusion.py`, `evaluate_gse94533_diagnostic_model.py`, `evaluate_gse94533_enhancements.py`, and the public-data runner.
- Historical PMML/eLife reconstruction and split investigation: `audit_original_pmml_model.py`, `extract_original_14mirna_model.py`, `reproduce_original_*`, `compare_pmml_to_published_table2.py`, `validate_pmml_against_table2_by_id.py`, and `search_*split*` scripts.
- Historical/current fixed-signature baselines: `build_14mirna_benchmark.py`, `build_modern_14mirna_baseline.py`, `build_fair_historical_14mirna_baseline.py`, `build_full_historical_mirna_baseline.py`, and age/CA-125 variants.
- GSE106817 extraction, audit, external evaluation, and within-cohort discovery: `audit_gse106817.py`, `extract_gse106817_validation_cohort.py`, `evaluate_gse106817_external.py`, `build_gse106817_170mirna_baseline.py`, `build_gse106817_allmirna_borderline_discovery.py`, and cancer-versus-borderline variants.
- Methylation: `evaluate_methylbert_methylation.py`, `freeze_methylbert_pipeline.py`, and `audit_methylbert_external_targets.py`.
- TCGA patient-linked stage/multimodal/survival experiments: `download_gdc_tcga_ov.py`, `build_tcga_*`, `train_tcga_ov_stage_model.py`, and audits.
- GSE8841/GSE12172 frozen gene-signature experiments: `build_gse8841_*`, `evaluate_frozen_signature_threshold.py`, `test_*signature*`, and mapping/diagnostic helpers.
- Mendeley/GSE9899 exploratory ingestion and diagnosis scripts: `download_*`, `extract_mendeley_*`, `audit_mendeley_*`, `mendeley_clinical_baselines.py`, `gse9899_gene_diagnostics.py`, and related helpers.
- Demos/smoke tests: `demo_ovarian_cases.py`, `ovarian_research_demo.py`, and `smoke_test.py`.

The main duplication is not duplicated data rows but duplicated generations of historical and frozen-signature scripts: `clean`, `standalone`, `final`, and `*_by_id` variants often implement the same scientific question with different path assumptions and output names. They should be treated as legacy until one canonical command and artifact are designated.

## Required Follow-up Before Model Optimization

1. Add a manifest that names one canonical script, input hashes, package versions, seed, fold assignments, endpoint, threshold policy, and output artifact for every reported result.
2. Move feature-availability filtering and any other data-dependent preprocessing inside outer training folds.
3. Add an explicit sample/patient overlap check between every training and external cohort, including GSE94533 versus GSE106817 and all frozen-signature transfers.
4. Reproduce the historical 44-patient test IDs from primary evidence or permanently downgrade the historical number to apparent/compatible-split evidence.
5. Report external GSE106817 ROC-AUC, balanced accuracy, sensitivity, specificity, confusion matrix, and prevalence together; suppress raw accuracy as a success criterion for this cohort.
6. Do not construct serum-miRNA plus methylation or TCGA multimodal rows without verified patient-level linkage and compatible endpoints.
7. Only after these controls are in place, evaluate early-stage performance and calibration on genuinely held-out, feature-aligned cohorts.
