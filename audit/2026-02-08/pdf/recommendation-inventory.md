# PDF recommendation inventory

Source: `Claude Methodology and Errors for Ovarian Cancer Project.pdf`. References below are physical PDF pages (1-based), not section numbers. IDs are stable within this inventory. Every item is a **document recommendation/claim, not an independently verified finding or authorization**. No recommended action was executed. Repeated recommendations retain separate page-local IDs to preserve coverage. Descriptive-only pages are explicitly identified.

## PDF page 1
- PDF-P01-R01: Position value around moving diagnosis earlier; avoid duplicating existing CA-125/ultrasound workup unless demonstrating comparator superiority or a genuinely earlier decision point (continues p2).
- Claims: epidemiology, survival, screening-trial conclusions, FDA indications/warnings require independent clinical/regulatory verification.

## PDF page 2
- PDF-P02-R01: Explicitly choose and communicate the intended-use tier; describe workflow friction rather than simply a missing algorithm.
- PDF-P02-R02: For adnexal-mass triage, obtain head-to-head/non-inferiority comparator evidence; for symptomatic pre-imaging triage, prospective symptom-cohort evidence (table continues p3).

## PDF page 3
- PDF-P03-R01: Lead with Tier 1, present Tier 2 as a funded next step, and exclude average-risk screening as a near-term competition claim.
- Claims: GSE94533 clinical suitability and table cohort counts are unverified against current assets.

## PDF page 4
- No new explicit recommendation; continuation of data-source/modality/patient-linkage table. All counts, independent-cohort assertions and modality descriptions remain source claims.

## PDF page 5
- PDF-P05-R01: Never concatenate unrelated-study rows into fabricated multimodal patients; demonstrate this discipline to judges.
- PDF-P05-R02: State endpoint and patient subset beside each headline metric; prepare to explain differences across pipelines (table continues p6).

## PDF page 6
- PDF-P06-R01: Present full-cohort and narrower early-stage results together with clear labels, exclusions and N.
- PDF-P06-R02: Consider near-zero interaction-term AUC benefit and sensitivity/specificity tradeoff before choosing the demo model.
- PDF-P06-R03: Disclose frozen-threshold external failure alongside ranking AUC (continues p7).
- Claims: numerical metrics, nested-CV legitimacy, and platform-scale causation are not verified here; preserved ranking does not independently establish biological validity.

## PDF page 7
- PDF-P07-R01: Present external-transfer failure as a named, understood, fundable open scientific problem.
- PDF-P07-R02: Keep independent methylation results separate from any claim of miRNA improvement.
- PDF-P07-R03: Never cite TCGA tumor staging as cancer-detection evidence.
- PDF-P07-R04: Describe gated fusion as architecture built/unit-tested awaiting patient-linked data, not a trained model.
- PDF-P07-R05: Highlight claimed methodological rigor: nested CV, fixed seeds, anti-fabrication rule and disciplined rejection of unlinked CPTAC data (continues p8).

## PDF page 8
- PDF-P08-R01: Highlight documented external-cohort screening and separate borderline category.
- PDF-P08-R02: Phase 0: choose a canonical headline and show the full endpoint table.
- PDF-P08-R03: Add production-model calibration plot.
- PDF-P08-R04: Add decision-curve analysis against refer-all/refer-none baselines.
- PDF-P08-R05: Repeat evaluation with feature selection inside each outer fold; report beside frozen-signature result.
- PDF-P08-R06: Move global feature-availability filters into training folds.
- PDF-P08-R07: Audit patient/sample-ID overlap with each external cohort.
- PDF-P08-R08: Regenerate subgroup error analysis with every model change (continues p9).

## PDF page 9
- PDF-P09-R01: Make histology-stratified sensitivity a first-class metric.
- PDF-P09-R02: Build a canonical run manifest: script, command, seed and output for every number.
- PDF-P09-R03: Apply for dbGaP phs000897 access early; claimed availability and approval timelines need verification.
- PDF-P09-R04: Acquire/collect an HE4-linked cohort for comparator evaluation.
- PDF-P09-R05: Investigate harmonization only with real batch metadata; consider requesting original run metadata from authors.
- PDF-P09-R06: Build multimodal rows only from same-patient, same-visit measurements; otherwise retain independent benchmarks.
- PDF-P09-R07: Consider longitudinal CA-125/miRNA trends only when repeated-draw data exist.
- PDF-P09-R08: Design a prospective single-platform Tier-1 adnexal-mass study.

## PDF page 10
- PDF-P10-R01: Consider a single CLIA reference-lab LDT track and/or FDA IVD track; plan model-update strategy/PCCP early.
- PDF-P10-R02: Align intended-use and reimbursement strategy to supported workup indications rather than unsupported screening.
- PDF-P10-R03: Seek regulatory consultant or FDA Q-Sub input before designing prospective evidence collection.
- Claims: LDT legal status, dates, regulatory options, PCCP guidance and payer coverage are time-sensitive, unverified claims, not legal advice.

## PDF page 11
- PDF-P11-R01: Fit testing into existing blood draw/order workflow as reflex/parallel add-on, not a new visit.
- PDF-P11-R02: Return structured results through HL7/FHIR into the existing EHR rather than a separate portal.
- PDF-P11-R03: Map result/referral decision support to relevant SGO/ACOG guidance, potentially CDS Hooks, rather than bare binary labels.
- PDF-P11-R04: Use existing serum collection/reference-lab logistics, subject to analytical validation.
- PDF-P11-R05: Lead with software/reference-lab path; describe point-of-care qRT-PCR hardware as a longer-horizon option requiring capital/validation.
- PDF-P11-R06: Prepare one-page results scorecard with N, endpoint and external failure (continues p12).

## PDF page 12
- PDF-P12-R01: Include funded external-gap closure plan on scorecard.
- PDF-P12-R02: Demo single/bulk portal only after exposing endpoint/N and imputation flags.
- PDF-P12-R03: Prepare regulatory/reimbursement and workflow diagrams.
- PDF-P12-R04: Tie funding ask to dbGaP/HE4 acquisition, prospective replication and CLIA lab partnership.
- PDF-P12-R05: Rehearse transparent external-performance explanation, without treating proposed platform-shift causation as independently proven.
- Descriptive glossary begins here, continues p13.

## PDF page 13
- PDF-P13-R01: Choose headline but always show full result table.
- PDF-P13-R02: Add calibration and decision curves.
- Descriptive glossary reviewed: sensitivity, specificity, balanced accuracy, OOF/CV, external validation, calibration, net benefit. Simplifications are source explanations, not a full methodological standard.

## PDF page 14
- PDF-P14-R01: Expose endpoint/N and bulk imputation flags in UI.
- PDF-P14-R02: Label every fusion-network mention as untrained architecture.
- PDF-P14-R03: Prepare regulatory/reimbursement and workflow diagrams.
- PDF-P14-R04: Start dbGaP access application; only say submitted/in progress when actually evidenced.
- PDF-P14-R05: Consult underlying technical-findings companion (embedded beginning p15).

## PDF page 15
- PDF-P15-R01: State N and exact endpoint beside every percentage and show full/narrow result tables together.
- Claims: reported script provenance, model metrics and likely same-cohort feature selection are allegations requiring repository/provenance review; opening cites an audit dated **2026-09-19**, later than this audit folder's 2026-02-08 label.

## PDF page 16
- PDF-P16-R01: Do not present frozen-signature AUC as unbiased new-patient performance if cohort informed discovery.
- PDF-P16-R02: Repeat feature selection within every outer fold using interaction feature space; use nested estimate as canonical.
- Claims: performance trajectory is consistent with leakage but does not by itself prove selection provenance or causal leakage.

## PDF page 17
- PDF-P17-R01: Alternatively validate fixed signature on a genuinely untouched cohort.
- PDF-P17-R02: Revisit deployed interaction/no-interaction model choice; explain clinical sensitivity/specificity tradeoff, not just AUC.
- PDF-P17-R03: Move full-frame availability/variance filters inside each fold.
- PDF-P17-R04: Do not present frozen-threshold external classification as usable despite high imbalance-driven accuracy.
- Source's “confirmed” wording is not independent confirmation by this audit.

## PDF page 18
- PDF-P18-R01: Exclude historical >80% independent-test claim unless original split is recovered; matching AUC alone is insufficient.
- PDF-P18-R02: Do not cite unsuccessful TCGA staging as detection evidence.
- PDF-P18-R03: Source recommends initializing git, committing and tagging metric-producing versions. **Historical/context-mismatched**: active clone has `.git/`; no commits or git changes authorized/executed.
- PDF-P18-R04: Choose an explicit license consistent with business/open-source posture.
- PDF-P18-R05: Add automated scientific guardrail tests (continues p19).
- Claims: absence of tests/license/history in delivered snapshot is unverified for current clone.

## PDF page 19
- PDF-P19-R01: Test zero patient overlap, fold-contained preprocessing and manifest/model feature dimensionality agreement.
- PDF-P19-R02: Route reproducing scripts through existing configuration paths rather than hardcoded home directories.
- PDF-P19-R03: Preserve dependency lock snapshot with artifacts and record sklearn/numpy/joblib versions in manifest.
- PDF-P19-R04: Choose one canonical script per result, archive variants and create RESULTS_INDEX mapping numbers to script/command/output/N/endpoint. No moves executed.
- Claims: 45/100 hardcoded scripts, artifact metadata omissions and duplication require current-code verification.

## PDF page 20
- PDF-P20-R01: Prioritize HE4-linked data acquisition/collection for comparator testing.
- PDF-P20-R02: Report histologic-subtype performance routinely and name generalization across subtypes as an open question.
- PDF-P20-R03: Seek Stage-I-specific metadata/data rather than implying separate Stage-I evaluation.
- PDF-P20-R04: Disclose CA-125 missingness prominently.
- PDF-P20-R05: Distinguish fusion scaffold from real-data training (continues p21).
- Claims: missingness counts, subtype associations and assertion “not a data artifact” are unverified; biological causation is not established here.

## PDF page 21
- PDF-P21-R01: State scaffold is unit-tested, not yet trained.
- PDF-P21-R02: Put N=113, endpoint composition and external-validation limitation next to UI benchmarks; wording must reflect current evidence rather than imply unbiased validation.
- PDF-P21-R03: Consider high-sensitivity/high-specificity operating points with displayed tradeoffs, rather than a single fixed cutoff; requires validation before clinical use.
- PDF-P21-R04: Add calibration/confidence information (implicit recommendation from identified absence).
- PDF-P21-R05: Flag default/imputed Age/CA-125 in bulk output.
- PDF-P21-R06: Add prediction audit logging with appropriate governance.
- PDF-P21-R07: Document human-subjects/IRB status (continues p22); source's “very likely fine” is not an ethics/legal determination.

## PDF page 22
- PDF-P22-R01: Add explicit public-data human-subjects/IRB statement grounded in appropriate review.
- PDF-P22-R02: Document citations/reuse terms per source in DATA_SOURCES.md; verify instead of assuming.
- PDF-P22-R03: Choose code license.
- PDF-P22-R04: Priority list repeats: full endpoint table; UI endpoint/N and imputation; scaffold labeling; version provenance; RESULTS_INDEX; fold-contained selection; calibration/DCA; dbGaP application; interaction-model tradeoff; license. Each remains a proposal, not execution approval.

## PDF page 23
- PDF-P23-R01: Compare with other published circulating-miRNA panels, not only historical internal signature.
- PDF-P23-R02: Add validated probability calibration (Platt/isotonic suggested) before describing output as literal clinical risk.
- PDF-P23-R03: Foreground external platform-transfer failure as a scientific result, preserving distinction between observed metrics and causal/generalizability interpretation.
- PDF-P23-R04: Prioritize prospective harmonized multi-site collection over further post-hoc correction; make it central to funding/collaborator ask.
- Claims: all attempted corrections failed, publishability and prospective superiority are source assessments, not verified conclusions.
