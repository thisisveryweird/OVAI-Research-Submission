# Frozen GSE94533 OVAI External Validation: 2024 FirePlex

## Preflight Result

**No patient-level predictions were generated.** The frozen production model requires all five production miRNAs, age, and CA-125. The workbook has age and CA-125 but has zero exact measured matches to the five required production miRNAs. Substituting historical markers, using fuzzy aliases, refitting a model, or estimating missing features would not be a frozen validation.

- 2024 patient records: 678.
- Exact production-miRNA coverage: 0/5.
- Rows eligible for frozen production-model scoring: 0.
- Frozen model, scaler, defaults, thresholds, and canonical artifacts were not modified.
- No external labels were used for feature selection, threshold tuning, or model selection.
- This is an input-compatibility stop, not an estimate of OVAI clinical performance on this cohort.

## Required Frozen Inputs

| Required production feature | Exact FirePlex match | Measured |
|---|---:|---:|
| `hsa-miR-200b-3p` | `no` | `no` |
| `hsa-miR-296-3p` | `no` | `no` |
| `hsa-miR-6826-5p` | `no` | `no` |
| `hsa-miR-32-3p` | `no` | `no` |
| `hsa-miR-1255a` | `no` | `no` |

The workbook reports FirePlex MFI values. The production OVAI pipeline was fitted using GSE94533 miRNA values with `log1p(max(value, 0))`, clinical imputation, scaling, and fixed interactions. There is no documented validated transform bridging these assay scales, and the required miRNA features are absent in any case.

## Disposition

Frozen external validation is **not computable from this workbook**. Do not report AUC, sensitivity, specificity, calibration, or threshold performance for the production OVAI model on these data. A future transportability study would first require a cohort that measures the exact five analytes or an independently validated assay bridge; neither can be manufactured from the current 2024 matrix.
