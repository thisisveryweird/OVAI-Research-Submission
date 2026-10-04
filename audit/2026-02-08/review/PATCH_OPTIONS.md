# Approval choices — proposals only

**Recommend option B.** No option authorizes runtime installs, retraining, threshold/objective selection, external acquisition, clinical use, publication, file moves, deletes or commits. Existing `production_artifacts/` and saved outputs remain protected; new artifacts need distinct versioned paths and an explicit later training approval.

## A — containment-only research demonstration
Smallest scope: fix fail-fast errors and disclosures, disable incompatible inference and mark old benchmarks conditional. Does not repair scientific evaluation or establish risk calibration.

- `scripts/train_mirna_baseline.py` — define/test intended patient-ID and threshold helpers or fail explicitly; never optimize threshold on final test labels.
- `scripts/demo_ovarian_cases.py` — replace undefined `x` lookup with validated row ID and reject duplicate keys.
- `scripts/predict_ovarian_risk.py` — reject schema mismatch before scaler use instead of silently constructing seven columns for nine-input artifacts.
- `scripts/ovarian_diagnostic_portal.py` — bind result to input/artifact snapshot, disclose source-only endpoint/N and per-row defaults, remove clinical-risk implication.
- `Launch_Diagnostic_Portal.bat` — replace silent per-launch pip install with explicit environment readiness check.
- `Launch_Diagnostic_Portal.command` — replace silent per-launch pip install with explicit environment readiness check.
- NEW `tests/research_inference/` — synthetic parser/schema/snapshot/default tests only, no clinical rows or model loads.

## B — versioned fold-safe pipeline plus safe research inference
Implement A's containment plus shared fit/transform semantics and locked-schema versioned pipeline. Keep fixed-panel conditional evaluation separate from complete discovery evaluation. Imputation order must be a reviewed design decision and identical between evaluation/export/inference; no artifact is replaced by this proposal.

- NEW `src/research_pipeline/eligibility.py` — sklearn fit-only training eligibility with stable output names and explicit all-empty policy.
- NEW `src/research_pipeline/features.py` — shared row-local logs and interactions, with documented order and negative/nonfinite guards.
- NEW `src/research_pipeline/schema.py` — ordered mature-arm/unit schema, duplicate normalization rejection and artifact dimension/version checks.
- NEW `src/research_pipeline/evaluation.py` — nest preprocessing, eligibility, selection, model tuning and whole-pipeline calibration inside appropriate training folds; preserve fold-specific decisions.
- NEW `src/research_pipeline/inference.py` — research-only adapter requiring complete schema/version, effective-input flags and immutable result snapshot.
- NEW `tests/research_pipeline/` — synthetic fit-containment, dimensionality, all-missing, parity, stage-parser, overlap and probability-semantics tests.
- NEW `production_artifacts/research_v2/` — reserved distinct future artifact namespace; create artifacts only after separate training authorization, never overwrite old joblib files.
- `scripts/evaluate_gse94533_diagnostic_model.py` — replace global eligibility with shared fold-fit transformer.
- `scripts/evaluate_gse94533_enhancements.py` — calibrate whole preprocessing pipeline and persist tuned confusion from actual fold decisions with aligned feature names.
- `scripts/run_reproducible_baseline.py` — isolate whole-pipeline calibration/threshold tuning and make summaries agree with saved per-fold decisions.
- `scripts/run_phase8_model_experiments.py` — move supervised selection and fitted preprocessing inside inner searches and meta-OOF training.
- `scripts/test_high_capacity_models.py` — share feature/imputation order, preserve fixed-panel label and separately evaluate fold-local discovery only if authorized.
- `scripts/export_production_pipeline.py` — export shared complete versioned pipeline/schema/hashes without writing legacy artifact names.
- `scripts/predict_ovarian_risk.py` — call shared schema-checked inference instead of unrelated seven-column feature construction.
- `scripts/ovarian_diagnostic_portal.py` — use explicit selected research artifact version and shared inference, effective inputs/default flags, stale-state prevention and bounded labels.
- `scripts/evaluate_methylbert_methylation.py` — fold-fit eligibility and keep independent processed methylation endpoint.
- `scripts/freeze_methylbert_pipeline.py` — share methylation preprocessing and safe post-transform selector bounds without automatic fitting.
- `scripts/train_tcga_ov_stage_model.py` — fold-fit eligibility and assert patient groups while retaining tumor-staging-only interpretation.
- `scripts/bulk_test_external_datasets.py` — exclude training accession, adjudicate stage/diagnosis/borderline labels and label collapse from confusion counts.
- `scripts/test_frozen_signature_stage1_gse12172.py` — repair undefined labels and impossible stage mask; distinguish explicit I from I/II.
- `scripts/test_frozen_signature_stage1_gse12172_final.py` — version/test exact IA/IB/IC parser and endpoint roster before comparison.
- `scripts/test_frozen_signature_gse9899.py` — separate target-normalized score adaptation from source-frozen validation, with polarity documented.
- `research/cohort_expansion/audit_gse211692.py` — replace hardcoded overlap zeros with roster-based bounded checks and unknown donor independence.
- `research/transportability_analysis/run_transportability_analysis.py` — report target-only feature eligibility and target-trained development separately; never select best AUC as frozen transfer success.
- `src/data/multimodal_table.py` — reject invalid IDs/labels and require externally verified donor/visit linkage metadata, not strings alone.
- `src/models/gated_fusion.py` — reject all-absent modality rows and assert finite values/shapes/masks before gating.
- `scripts/smoke_test.py` — replace vacuous assertion with meaningful synthetic gate/mask invariants; no clinical-training claim.
- `scripts/analyze_gse94533_errors.py` — consume actual saved thresholds/decisions and retain subtype denominators, fold/model/run IDs and distinct Stage-I case-only outputs.
- NEW `RESULTS_INDEX.md` — only after owner approval, map each approved number to source, command, output, N, endpoint and immutable hashes; no canonical metric chosen here.
- NEW `tests/metadata/` — source-only manifests/licenses/path inventory checks, not clinical tests.

Historical scripts require a separate primary-header and original-split decision: `scripts/reproduce_original_historical_clean.py` — make exported cohort provenance match119/60 algorithm; `scripts/reproduce_original_historical_standalone.py` — remove incompatible120/59 assumption only after primary workbook review; `scripts/reproduce_original_pmml_historical_data.py` — require explicit header translation and preserve tanh/argmax semantics. No matching-AUC search is accepted as recovered split evidence.

## C — evidence-first freeze
No product edits. Obtain immutable clinical/assay/roster/provenance and published model terms, choose intended-use and validation protocol with owner/clinical/regulatory partners, then revisit B. Safest if demonstration can wait; slower engineering feedback. Use existing `configs/project.yaml` — plan resolved configuration paths in a future approved portability patch; NEW `docs/research_protocol/` — record proposed donor-disjoint endpoints, selection/threshold/calibration plans and legal/IRB determinations. These new documents must not declare approvals or source permissions that have not been obtained.
