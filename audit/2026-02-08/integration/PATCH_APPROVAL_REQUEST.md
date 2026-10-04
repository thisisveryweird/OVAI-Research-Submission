# Patch approval request — authorization=false; ALL PENDING

No patches applied. No winner, threshold, license or tier selected. Each tranche is independently bounded: approval must specify exact paths, permitted behavior changes and exclusions; no acquisition/fitting/deployment implied.

## A inference parity and research-only UI
- `scripts/test_high_capacity_models.py` (exists)
- `scripts/export_production_pipeline.py` (exists)
- `scripts/predict_ovarian_risk.py` (exists)
- `scripts/ovarian_diagnostic_portal.py` (exists)
- `scripts/demo_ovarian_cases.py` (exists)
- `scripts/ovarian_research_demo.py` (exists)
Acceptance after approval: review diff; confirm protected files unchanged; use prespecified synthetic/archived fixtures and schema/fold/threshold/source hash guards; clinical validation requires separate protocol.

## B fold-contained scientific evaluation
- `scripts/evaluate_gse94533_diagnostic_model.py` (exists)
- `scripts/evaluate_gse94533_enhancements.py` (exists)
- `scripts/run_reproducible_baseline.py` (exists)
- `scripts/evaluate_methylbert_methylation.py` (exists)
- `scripts/train_tcga_ov_stage_model.py` (exists)
- `scripts/run_phase8_model_experiments.py` (exists)
Acceptance after approval: review diff; confirm protected files unchanged; use prespecified synthetic/archived fixtures and schema/fold/threshold/source hash guards; clinical validation requires separate protocol.

## C linkage and scaffold safeguards
- `src/data/multimodal_table.py` (exists)
- `src/models/gated_fusion.py` (exists)
- `scripts/smoke_test.py` (exists)
- `scripts/build_gse94533_clinical_fusion.py` (exists)
- `scripts/add_age_ca125_to_14mirna.py` (exists)
- `scripts/rebuild_clinical_benchmark.py` (exists)
Acceptance after approval: review diff; confirm protected files unchanged; use prespecified synthetic/archived fixtures and schema/fold/threshold/source hash guards; clinical validation requires separate protocol.

## D claim/provenance-only corrections
- `README.md` (exists)
- `START_HERE.txt` (exists)
- `configs/project.yaml` (exists)
- `requirements.txt` (exists)
- `requirements-deep.txt` (exists)
- `Launch_Diagnostic_Portal.bat` (exists)
- `Launch_Diagnostic_Portal.command` (exists)
- `reports/bulk_external_test_results.md` (exists)
- `reports/high_capacity_model_results.md` (exists)
Acceptance after approval: review diff; confirm protected files unchanged; use prespecified synthetic/archived fixtures and schema/fold/threshold/source hash guards; clinical validation requires separate protocol.

## E external eligibility and research truthfulness
- `scripts/bulk_test_external_datasets.py` (exists)
- `scripts/evaluate_gse106817_phase5.py` (exists)
- `research/cohort_expansion/audit_gse211692.py` (exists)
- `research/cohort_expansion/audit_cross_platform_replication.py` (exists)
- `research/transportability_analysis/run_transportability_analysis.py` (exists)
Acceptance after approval: review diff; confirm protected files unchanged; use prespecified synthetic/archived fixtures and schema/fold/threshold/source hash guards; clinical validation requires separate protocol.

## Genuine architecture alternatives (decision not made)
1. Recommended for consideration: separately versioned shared fold-safe sklearn Pipeline containing availability, imputation, transforms, interaction construction, selection and estimator; calibrate the entire pipeline within appropriate splits. New version and immutable manifests; preserve legacy artifacts; establish parity or explicitly different scientific estimand.
2. Retain legacy artifact via a narrow safety wrapper enforcing nine-feature schema, exact legacy transformation and immutable artifact hash; reject invalid/ambiguous/all-missing inputs, expose imputation and research-only outputs. Preserves compatibility but cannot repair training provenance or establish clinical validation; seven-input CLI must not silently adapt.
3. Supported replacement path: approved maintained assay-specific pipeline with authenticated units/linkage, pinned supported packages, migration protocol and prospective untouched validation. Legacy kept as reference; no automatic winner selection or threshold transfer.

Signatories: repository owner for code, scientific lead for estimand, data owner for linkage/reuse, institutional IRB/privacy reviewer, clinical lead and regulatory counsel for intended use. Evidence and postapproval tests are in BLOCKED_RESEARCH_AND_DECISIONS.md. No approval received.
