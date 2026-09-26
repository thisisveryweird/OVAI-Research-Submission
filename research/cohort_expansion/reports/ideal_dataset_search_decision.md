# Additional Cohort Search and Decision

## Search Scope and Findings

The search prioritized individual-level circulating-miRNA data from patients presenting with an adnexal/pelvic mass, benign-mass controls, early stage, and quantitative measurements. The existing GEO metadata scout screened 58 candidate series (48 uncached candidates; none passed all prespecified metadata laws), and the remaining GEO queue report found zero candidates suitable for automated cache triage. A direct GEO and SRA search for the 2019 pelvic-mass article returned no accession. PLOS Figshare resolves the article’s supplementary materials but only to assay sequences and aggregate cohort tables. Zenodo returned no article-level data record; Dryad lookup was unavailable from this environment.

Previously audited GSE106817 has large title-derived group counts (3,108 samples; 320 cancer, 29 benign, 2,759 control, 66 borderline, 872 other) but has extensive assay-distribution shift and a failed source-derived threshold in the prior external evaluation. It is not promoted as the ideal clinically matched cohort, and no new transportability experiment was run here.

## Best Target to Pursue

The **Danish 380-patient pelvic-mass cohort** is the best next clinical-data access target: it directly contrasts malignant with benign masses, is presurgical plasma, reports age/CA-125/stage/histology at study level, and measures 46 candidate miRNAs. It is not currently a usable individual-level dataset: public files are aggregate only, Stage I/II counts are just 13, and only one current OVAI feature (miR-200b-3p) is present. Requesting the de-identified patient-level matrix, exact sample inclusion/assay availability, and labels is the appropriate next move.

The 2024 serum FirePlex workbook is the best newly accessible patient-level cohort for **cross-platform biological replication**, with 678 records, 56 unambiguous Stage I/II cancer IDs, and eight historical-miRNA overlaps. Its controls pool benign and healthy participants and none of the current five production miRNAs are measured. It is not an expanded OVAI training cohort.

No GREEN cohort for frozen current-OVAI external validation is presently established. The 380-patient study is YELLOW pending data access and remains feature-incomplete for the current 5-miRNA model; the 2024 workbook is YELLOW for biological replication only; JanusRNA is access-controlled YELLOW for prediagnostic biology.

## Early-Stage and Benign-Mass Capacity

For the 380-patient study, the published table reports 13 Stage I/II cancers (5 Stage I, 8 Stage II) and 190 enrolled benign pelvic masses. After the paper’s stated exclusions (14 cancers and 13 benign samples), no more than 177 benign samples and at most 13 early-stage cancers could remain: an upper bound of 190 early-stage-plus-benign rows. The actual retained Stage I/II count and feature-complete denominator are unknown without individual-level data; the validation cohort itself has zero Stage I/II cases. For current frozen OVAI, the usable count remains zero because four of five production miRNAs are not in its assay panel.

## Answers to the Phase Questions

1. The 2024 cohort supports feature-wise biological replication and external cohort characterization, not direct evaluation or expansion of the frozen OVAI model. A distinct FirePlex-specific model would be a separate future project, not performed here.
2. Exact 2024 overlaps among the historical 14 are: `hsa-miR-23b-3p`, `hsa-miR-29a-3p`, `hsa-miR-32-5p`, `hsa-miR-92a-3p`, `hsa-miR-150-5p`, `hsa-miR-200c-3p`, `hsa-miR-320c`, `hsa-miR-320d`.
Published 2024 feature-list exact matches in the 2024 workbook: `hsa-miR-200c-3p`, `hsa-miR-320c`, `hsa-miR-320e`, `hsa-miR-195-5p`, `cel-miR-39-3p`, `hsa-miR-770-5p`, `hsa-miR-122-5p`, `hsa-miR-19a-3p`, `hsa-miR-424-5p`; absent exact identifier: `hsa-miR-103-3p`. The measured hsa-miR-103a-3p column is not silently recoded to the paper’s hsa-miR-103-3p label in this exact-match count.
3. No valid alias for any of the five current production features was established against the measured 2024 panel. The measured miR-32-5p is a different mature arm from miR-32-3p; miR-200a/200c are not miR-200b.
4. The 380-patient article has no public individual-level measurements identified. Public Figshare files are aggregate tables; individual data require an author data request.
5. The 380 panel contains 1/5 exact current OVAI features: hsa-miR-200b-3p.
6. It is not a viable frozen-OVAI external validation cohort as currently documented (only 1/5 features and no row-level data); it is a high-priority clinically aligned access request and biological replication target.
7. JanusRNA can address prediagnostic timing after controlled-access approval, but not the intended benign-mass diagnostic endpoint from the published case-only subset.
8. Best clinically matched target to pursue: the Danish 380-patient pelvic-mass study. Best additional patient-level dataset currently in hand: the 2024 678-patient workbook for biological replication only.
9. The 380 cohort’s endpoint and presurgical specimen are clinically compatible; its plasma RT-qPCR platform, public-data limitation, and 1/5 feature overlap prevent direct frozen-model validation.
10. Published ceiling: 13 early-stage cancers plus at most 177 retained benign masses (<=190 total after the stated technical exclusions). Actual complete rows and retained early-stage cases are unverified; current full OVAI scoreable count is zero.

No model was trained, optimized, or modified in this phase. No new cohort was merged with GSE94533.
