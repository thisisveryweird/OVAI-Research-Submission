# Transportability analysis summary

This research-only analysis uses only the saved OVAI artifacts and does not refit the canonical model or adjust any threshold on external labels. The canonical source artifact remains the frozen production pipeline recorded in [production_artifacts/manifest.json](../../production_artifacts/manifest.json).

## Summary

- Source model: GSE94533 Stage I/II cancer vs benign/control logistic model, encoded in [production_artifacts/manifest.json](../../production_artifacts/manifest.json).
- External stress-test cohort: GSE106817, evaluated using the saved phase-5 results in [reports/gse106817_phase5_results.json](../../reports/gse106817_phase5_results.json).
- Key conclusion: the external failure is not primarily a lack of ranking signal; it is a combination of distribution shift and a probability/threshold mismatch. The ranking signal remains measurable, but the frozen source threshold collapses under cohort shift.

## Evidence from the saved outputs

The best retained GSE94533 -> GSE106817 ranking signal in the saved analysis is:

- normalization: training_derived_standard
- ROC-AUC: 0.7817
- balanced accuracy: 0.4998
- sensitivity: 0.0000
- specificity: 0.9996
- threshold: 0.60

The strongest evidence of threshold failure is:

- normalization: raw_values
- threshold: 0.69
- sensitivity: 0.0000
- specificity: 1.0000
- predicted positive rate: 0.0000

This is exactly the pattern expected when the score distribution shifts but the ranking relationship remains partially preserved.

## Distribution shift diagnostics

The saved shift audit in [reports/gse106817_platform_shift_audit.json](../../reports/gse106817_platform_shift_audit.json) and [reports/gse106817_phase5_distribution.csv](../../reports/gse106817_phase5_distribution.csv) shows severe feature-level shift. The top five KS shifts are:

| feature | KS statistic | median shift (log1p) |
|---|---:|---:|
| hsa-miR-29a-3p | 1.000 | -5.458 |
| hsa-miR-92a-3p | 1.000 | -7.421 |
| hsa-miR-320c | 1.000 | -3.130 |
| hsa-miR-1307-5p | 1.000 | -2.226 |
| hsa-miR-320d | 0.987 | -2.366 |


## Interpretation

- The model still has measurable discrimination on GSE106817, especially after source-fitted scaling; the problem is not pure ranking failure.
- The canonical source threshold cannot be transported unchanged. It causes near-all-negative predictions in the external cohort.
- This is consistent with a platform/assay shift: the miRNA distributions differ substantially, as shown by KS statistics near 1.0 for several features in the audit.
- Therefore, the current computational evidence indicates that no normalization alone reliably fixes transportability without a standardized assay or a truly comparable validation cohort.

## Decision

Did we actually fix transportability, or did we only characterize it?

We only characterized it. The evidence supports a domain-shift bottleneck, not a solved transportability problem. The existing external cohorts are stress tests, not evidence that the canonical model is transportable across assay/platform boundaries.

## Reproducibility metadata

The manifest recorded by this analysis is in [research/transportability_analysis/transportability_manifest.json](transportability_manifest.json).

