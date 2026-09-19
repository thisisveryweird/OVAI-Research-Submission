# GSE106817 Platform-Shift Audit

This audit uses GSE94533 as the training cohort and keeps GSE106817 untouched. It does not fit a transform on external labels and does not add GSE106817 to training.

GSE106817 samples evaluated: 3108. Full title-derived groups: {'cancer': 320, 'benign': 29, 'borderline': 66, 'control': 2759, 'other': 872}.

## What Normalization Can and Cannot Do

The existing external evaluation already fits standard scaling and a source-fitted quantile transform on GSE94533 only. This audit quantifies feature distribution differences and records those frozen results. These transforms can alter numeric scale; they cannot prove that GPL21263 measurements are assay-equivalent to the GSE94533 processed serum miRNA measurements, recover missing batch covariates, or make an external cohort part of the training pool.

## Largest Distribution Differences

- `hsa-miR-29a-3p`: KS=1.000, median log1p shift=-5.458, source-fitted external median z=-6.934.
- `hsa-miR-92a-3p`: KS=1.000, median log1p shift=-7.421, source-fitted external median z=-10.024.
- `hsa-miR-1307-5p`: KS=1.000, median log1p shift=-2.226, source-fitted external median z=-4.621.
- `hsa-miR-320c`: KS=1.000, median log1p shift=-3.130, source-fitted external median z=-3.947.
- `hsa-miR-320d`: KS=0.987, median log1p shift=-2.366, source-fitted external median z=-2.939.

## Frozen External Result

The existing standard-scaling result has ROC-AUC 0.8232 but predicts every evaluated external sample as non-cancer at the source-derived threshold: sensitivity 0%, specificity 100%, balanced accuracy 0.5. The source-fitted quantile result has ROC-AUC 0.6893 and the same threshold failure.

## Decision

GSE106817 remains an untouched external stress test and platform-shift study. It is not added to training, and normalization alone is not accepted as a fix. A valid next experiment requires prespecified assay harmonization, a positive-control operating point defined without target labels, overlap auditing, and replication on another independent serum cohort.

Machine-readable output: `reports/gse106817_platform_shift_audit.json`.
