# 2024 Cohort Compatibility and Decision

## Data Integrity and Usable Cohort

The manually supplied workbook is a valid, intact XLSX (1,744,787 bytes; SHA256 `e98f5a60fe3117fc69cd9243bc60f146979b831b20dbb6dee4d8e7ce6affc748`). It contains 678 unique patient/sample IDs and 180 complete numeric assay columns (179 human-target columns plus one cel-miR-39 control). The `miRNA` and `protein + metadata` sheets have exact one-to-one ordered ID matches and matching diagnosis/split values.

The published cohort is reconstructed exactly: 286 cancer and 392 pooled controls; development 191/277, internal validation 44/56, and BWH external validation 51/59. The workbook does not distinguish benign adnexal masses from healthy controls, so subtype-specific control counts and cancer-versus-benign-only performance are unavailable.

## Stage and Control Limits

The histology worksheet has 145 rows but 141 unique IDs, all among cancer patients. Three rows duplicate an otherwise identical record; one ID has both Stage I/II and Stage III/IV entries and is classified as conflicting. After exact duplicate removal, 56 patients have unambiguous Stage I/II, 84 have unambiguous Stage III/IV, one is unresolved due to conflicting entries, and 145 cancers have no row in the histology sheet. Raw worksheet counts (59/86) reproduce the paper, but they are not distinct patient counts. Stages I vs II and III vs IV cannot be separately determined.

Thus the available early-stage comparison is 56 unambiguous Stage I/II cancers versus 392 pooled controls; this is not a benign-mass-only estimate. In the BWH external group, the workbook records 20 Stage I/II and 31 Stage III/IV cancer IDs against 59 pooled controls.

## OVAI Compatibility

The five production features are not present as exact measured mature-miRNA identifiers (0/5). The historical 14-feature signature has 8/14 exact matches; `hsa-miR-335-5p_0_-2` remains an unverified assay alias and is not counted. FirePlex MFI values are not pooled with GSE94533 sequencing values.

**Decision:** the 2024 workbook is not compatible with direct frozen evaluation of the current production OVAI model and is not a transport validation for that signature. No OVAI retraining, feature substitution, threshold tuning, or modeling was performed. See `2024_adnexal_mass_frozen_validation.md` for the scoring stop gate.

## Later Development Eligibility

The workbook is suitable for research-only development of a distinct FirePlex-specific model because it supplies the published development/internal/external labels and complete panel data. Preserve the 468/100/110 split; use only the 468 ASPIRA training subjects for any fitting/selection, and keep internal and BWH external validation untouched until a prespecified model is frozen. The control phenotype is pooled and cannot support benign-mass-only performance claims from this workbook.

The dataset may support a separately identified FirePlex-specific research model after prespecified methods are established. It must remain separate from canonical OVAI artifacts and from GSE94533 values. The published combined-model performance cannot be reproduced faithfully from the workbook alone; details are in `2024_adnexal_mass_published_model_reproduction.md`.
