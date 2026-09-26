# 380-Patient Pelvic-Mass Cohort Data-Access Audit

## Study and Clinical Fit

Study: "Diagnostic plasma miRNA-profiles for ovarian cancer in patients with pelvic mass" (Oliveira et al., 2019; DOI `10.1371/journal.pone.0225249`; PMCID `PMC6860451`). The study enrolled 380 presurgical plasma samples: 190 malignant ovarian/pelvic masses and 190 benign pelvic masses. It used a 46-candidate miRNA TaqMan RT-qPCR panel; the supplement enumerates 48 assays including two C. elegans spike-in controls. Discovery and validation cohorts were each reported as 95 malignant/95 benign.

The intended endpoint (malignant versus benign pelvic mass) is directly clinically compatible with OVAI’s triage objective. Specimen is plasma, not serum; platform is targeted TaqMan RT-qPCR with study-specific normalization, so no pooling with GSE94533 serum sequencing is appropriate.

## Patient-Level Availability

The article’s Data Availability statement says relevant data are in the paper and Supporting Information. Public PLOS Figshare record `10322942` attaches only S1-S3 DOCX tables and figures; S1 lists assay IDs/sequences, S2 gives group-level median/minimum/maximum expression summaries, and S3 gives model sensitivity/specificity. No per-patient matrix, sample IDs, or row-level clinical covariates are attached. Exact GEO and SRA article-title searches returned no records; the public Figshare search returned the article and aggregate supplements, not individual measurements. Zenodo search returned no record; Dryad endpoint lookup was unavailable from this environment. No public patient-level deposit was identified.

**Status: access required for individual-level analysis.** The publication does not provide an individual-level data-access route in its data statement; request the raw de-identified matrix and matching clinical labels from the corresponding authors. Do not treat the aggregate S2 table as an ML-ready cohort.

## Published Cohort Counts and Limitations

- Enrolled/reported: 380 total; 190 malignant and 190 benign; discovery 95/95; validation 95/95.
- Table 1 reports 13 Stage I/II cancers in discovery (Stage I=5, Stage II=8) and none in validation; validation cancer cases are Stage III/IV (72 Stage III, 23 Stage IV). Thus only 13 early-stage cases are reported, below the project’s preferred 50.
- The methods also state that 14 cancer and 13 benign samples were excluded for insufficient plasma or poor RNA quality. Per-patient retention and stage among excluded/retained subjects are not available. A simple upper bound after these exclusions is 177 benign plus at most 13 early-stage cancers, at most 190 candidate rows; actual complete-case numbers are likely lower and cannot be verified from summaries.
- Supplement S2 has feature-specific denominators and missing summaries; for example, hsa-miR-200b-3p has no discovery-cohort summary values in S2 but is summarized in validation. No patient-level missingness or exact usable count is recoverable.
- Cohort-level age and CA-125 summaries and histology/stage tables are published, but no individual-level age, CA-125, stage, or histology fields are released. Patient/sample IDs and duplicate/overlap audits are not possible.
- Samples were presurgical. The paper does not provide patient-level treatment status or timing sufficient to certify a treatment-naive analysis.

## Exact Feature Compatibility

The supplement lists 10 exact mature-feature mappings across the requested union of signatures, plus 1 documented name aliases.
- Current production OVAI: 1/5 exact, `hsa-miR-200b-3p` (TaqMan assay 002251; mature sequence `UAAUACUGCCUGGUAAUGA`). Four required production features are not in the panel, so the frozen 5-feature-plus-clinical OVAI model cannot be evaluated.
- Historical 14-feature signature: 6/14 exact panel IDs: miR-23b-3p, miR-29a-3p, miR-92a-3p, miR-200a-3p, miR-200c-3p, and miR-335-5p.
- 2024 paper’s 10-feature list: 4 exact panel IDs (miR-195-5p, miR-200c-3p, miR-122-5p, and miR-424-5p) plus a documented mature-name alias for miR-103-3p to the panel’s miRBase ID miR-103a-3p. The remaining five entries (miR-320e, miR-320c, cel-miR-39-3p, miR-770-5p, and miR-19a-3p) are absent. The panel has cel-miR-54-3p and cel-miR-238-3p, not cel-miR-39-3p.

## Decision

Classify as **YELLOW: clinical-endpoint biological replication / external cohort candidate, access required**. Endpoint and presurgical sampling are strong; individual measurements are not public, early-stage count is only 13, and only one of five frozen OVAI features overlaps. It is not currently a valid external validation of the canonical OVAI model and must not be merged with GSE94533.
