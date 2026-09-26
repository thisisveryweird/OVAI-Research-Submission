# Candidate cohort recommendation for OVAI early-stage expansion

## Executive decision

No green cohort for frozen current-OVAI external validation is currently established.

The most clinically relevant next target is the Danish 380-patient pelvic-mass study, but it remains a data-access project rather than a validated cohort. The public record is aggregate-only and the current production OVAI feature overlap is incomplete. The 2024 FirePlex workbook is valuable for cross-platform biological replication but should remain a separate research dataset, not a training or validation expansion for the canonical OVAI model. JanusRNA is a strong prediagnostic biology candidate, but it is not an adnexal-mass diagnostic comparison cohort.

## 1) Best candidate for expanding early-stage OVAI training

Status: not yet validated; highest-priority access target.

Candidate: Danish 380-patient pelvic-mass study (if patient-level data are obtained from the authors)

Why this is the best target:
- direct malignant-versus-benign pelvic/adnexal mass endpoint
- presurgical plasma samples
- reported stage distribution with 13 Stage I/II cancers in the discovery set
- strong clinical match to the OVAI triage problem
- broader demographic and symptomatic-control context than healthy-only cohorts

Current public limitation:
- no patient-level matrix identified in the public supplements
- only 13 Stage I/II cancers are reported in the published discovery cohort
- only 1/5 exact overlap with the current frozen OVAI feature set (`hsa-miR-200b-3p`)
- no public row-level CA-125, age, histology, or exact sample map

Publicly stated upper bound after exclusions:
- Stage I/II cancers: up to 13
- benign masses: up to 177
- total potential early-stage-plus-benign set: up to 190 rows before feature-complete filtering

Critical note: these are ceiling values from the summary paper; the actual usable row count and exact patient-level retention are unknown until individual-level access is granted.

## 2) Best candidate for benign-mass external validation

Status: the same candidate as above, but only if individual-level data are obtained.

Recommended target: Danish 380-patient pelvic-mass study

Reason:
- it is the only strong benign-mass clinical comparator identified in the research workflow
- it includes malignant and benign pelvic masses rather than healthy controls
- it matches the symptomatic cohort problem OVAI is intended to solve

Public shortfall:
- summary tables only
- validation cohort contains no Stage I/II cancers
- no patient-level data to support direct external validation

## 3) Best candidate for platform-shift testing

Candidate: GSE106817

Why:
- public serum cohort with a large sample size and a long-standing cross-platform comparison history
- relevant for stress testing assay and distribution shifts
- strong candidate for platform-change analysis, not training

Why not as a direct OVAI cohort:
- prior transportability analyses showed substantial platform and distribution shifts
- not clinically matched to the current benign-mass triage problem
- not a valid frozen-OVAI external validation target

## 4) Best prediagnostic candidate

Candidate: JanusRNA

Why:
- controlled-access serum cohort in a prediagnostic setting
- directly relevant to whether a miRNA signal exists before ovarian cancer diagnosis
- useful for timing and biological plausibility work

Why it is not the main training target:
- it is a prediagnostic timing study, not a symptomatic benign-mass comparison
- the public analysis is case-only and lacks benign controls
- data access is constrained and not immediately available

## 5) Which candidates require controlled access

- Danish 380-patient pelvic-mass study: controlled access or author-contact required for patient-level de-identified data
- JanusRNA: controlled access via the appropriate Data Access Committee
- Any dataset with only aggregate supplements and no public matrix should be treated as controlled-access or author-contact pending review

## 6) Which candidates are aggregate-only

- Danish 380-patient pelvic-mass study: aggregate-only in the public record
- 2024 FirePlex workbook: patient-level workbook exists, but benign-mass subtype information is not recoverable at patient level for the pooled controls
- GSE106817: aggregate/public summary context only for the current decision; not a valid patient-level external validation source

## 7) Which candidates are genuinely patient-level

- GSE94533 local training cohort: patient-level, internal benchmark
- 2024 FirePlex workbook: patient-level workbook exists and is usable for a separate FirePlex-specific study
- JanusRNA: patient-level in principle under controlled access
- Danish 380-patient study: not publicly patient-level; requires access negotiation

## 8) Exact number of additional Stage I/II patients potentially obtainable

From the published 380-patient study summary:
- Stage I/II cancers reported in the discovery cohort: 13 total
- Stage I: 5
- Stage II: 8

This is the maximum public early-stage count reported before patient-level feature-complete filtering. Because the paper states exclusions and omits individual retention details, the exact count of usable Stage I/II patients is not recoverable from public supplements.

Potential additional benign-mass controls from the public summary:
- up to 177 benign masses after stated exclusions

Important caveat:
- these are upper-bound counts, not confirmed usable rows
- no exact patient-level matrix is available

## 9) Exact overlap with the current OVAI feature space

For the candidate cohorts that can be assessed from public evidence:
- Danish 380-patient study: 1/5 exact current OVAI feature overlap (`hsa-miR-200b-3p`)
- 2024 FirePlex workbook: 0/5 exact current OVAI feature overlap
- JanusRNA: public summary does not establish exact current-OVAI feature overlap
- GSE106817: not validated against the current OVAI feature set in this workflow
- GSE94533: reference local benchmark; exact overlap is fully contextualized in the local training cohort only

## 10) Final recommendation

The winner for the immediate next step is not a large, public cohort; it is the clinically aligned, access-required Danish pelvic-mass study, provided patient-level de-identified data can be obtained and exact feature compatibility is confirmed.

The recommended sequence is:
1. Request patient-level data from the Danish 380-patient study.
2. Confirm exact assay compatibility with the current OVAI feature set and stage distribution.
3. Preserve the 2024 FirePlex workbook as biological replication only.
4. Keep JanusRNA in the prediagnostic bucket and do not mix it with the diagnostic benign-mass expansion target.
5. Maintain GSE106817 as a platform-shift stress-test resource, not a training source.

This is the dataset strategy that best addresses OVAI’s actual bottleneck: too few clinically representative Stage I/II cancers and benign-mass controls with compatible circulating-miRNA measurements.
