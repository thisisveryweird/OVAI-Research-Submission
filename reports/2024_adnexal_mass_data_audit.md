# 2024 serum adnexal-mass cohort audit

This audit preserves the canonical OVAI source model and does not use the 2024 adnexal-mass paper to overwrite any existing artifact.

## Summary

The paper title is: "Serum miRNA improves the accuracy of a multivariate index assay for triage of an adnexal mass" (Gynecologic Oncology, 2024; DOI: 10.1016/j.ygyno.2024.08.008).

The public paper metadata confirms the following study-level numbers:

- total subjects: 678
- cancer: 286
- controls: 392
- training subjects: 468
- independent internal validation subjects: 100
- independent external validation subjects: 110
- serum miRNA measurements: 180
- serum protein biomarkers: 7
- age: present
- menopausal status: present

## Current status of the raw data

The individual-level supplementary Data File S1 is not present in this workspace, and the public publisher-facing metadata available from the DOI/article page does not provide the raw patient-level table. The file therefore has not been audited at the row/subject level. No patient-level training or validation model has been fit using this cohort.

A full audit requires the actual Data File S1 or equivalent patient-level supplementary table. Until that is retrieved and reviewed, this dataset remains a literature-level candidate rather than an usable external validation dataset for OVAI.

## Questions that remain unresolved without the raw Data File S1

1. How many rows are in Data File S1?
2. Whether the file contains all 678 subjects or only a subset.
3. Whether individual subject IDs exist.
4. Whether the published training/internal/external split is reconstructible.
5. Which samples are benign adnexal masses vs healthy controls.
6. Which samples are invasive ovarian cancer.
7. Stage for each cancer.
8. Histology.
9. Age.
10. Menopausal status.
11. CA-125.
12. HE4.
13. APO-A1.
14. B2M.
15. transthyretin.
16. transferrin.
17. FSH.
18. Which 180 miRNAs are present.
19. Raw numeric scale and normalization.
20. Missing values.
21. Duplicate subjects or samples.
22. Same-patient repeats or serial collections.
23. Plate, batch, or collection-date information.
24. Patient overlap with GSE94533.

## Local evidence

Candidate files detected in the local workspace:

- none matching the published 2024 adnexal-mass supplementary Data File S1

This means the dataset cannot yet contribute to the patient-level OVAI development/validation pipeline in a leakage-safe manner.

## Decision

This cohort is scientifically relevant and should remain in the development queue, but it is not yet eligible for one of the following:

- direct pooling with GSE94533
- feature selection on external labels
- threshold tuning on external labels
- use as a canonical validation cohort

It should be treated as a high-priority acquisition target for a future, properly audited replication study.

## Reproducibility note

This audit is intentionally conservative. It records the paper-level summary and documents the access gap explicitly instead of inventing patient-level measurements.

