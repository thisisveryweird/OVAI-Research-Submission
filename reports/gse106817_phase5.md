# Phase 5: GSE106817 Bidirectional Validation

This evaluation uses only shared historical miRNAs and never uses external test labels for preprocessing, model selection, normalization, or threshold selection.

## Experiment A: GSE94533 -> GSE106817

- `raw_values`: ROC-AUC 0.6735, balanced accuracy 0.5000, sensitivity 0.0000, specificity 1.0000, threshold 0.69.
- `log1p_values`: ROC-AUC 0.7439, balanced accuracy 0.6694, sensitivity 0.7438, specificity 0.5951, threshold 0.64.
- `within_dataset_standard`: ROC-AUC 0.7571, balanced accuracy 0.6673, sensitivity 0.7625, specificity 0.5721, threshold 0.56.
- `training_derived_standard`: ROC-AUC 0.7817, balanced accuracy 0.4998, sensitivity 0.0000, specificity 0.9996, threshold 0.60.
- `within_dataset_quantile`: ROC-AUC 0.7160, balanced accuracy 0.6602, sensitivity 0.6250, specificity 0.6955, threshold 0.53.
- `training_derived_quantile`: ROC-AUC 0.7209, balanced accuracy 0.5000, sensitivity 0.0000, specificity 1.0000, threshold 0.48.

## Experiment B: GSE106817 -> GSE94533

- `raw_values`: ROC-AUC 0.5089, balanced accuracy 0.5052, sensitivity 0.0104, specificity 1.0000, threshold 0.45.
- `log1p_values`: ROC-AUC 0.6749, balanced accuracy 0.5971, sensitivity 0.3750, specificity 0.8193, threshold 0.48.
- `within_dataset_standard`: ROC-AUC 0.6783, balanced accuracy 0.5968, sensitivity 0.3021, specificity 0.8916, threshold 0.58.
- `training_derived_standard`: ROC-AUC 0.6788, balanced accuracy 0.5934, sensitivity 0.5000, specificity 0.6867, threshold 0.56.
- `within_dataset_quantile`: ROC-AUC 0.6798, balanced accuracy 0.6515, sensitivity 0.6042, specificity 0.6988, threshold 0.49.
- `training_derived_quantile`: ROC-AUC 0.5063, balanced accuracy 0.5000, sensitivity 1.0000, specificity 0.0000, threshold 0.48.

## Experiment C: Pooled Leave-One-Dataset-Out

Only two compatible serum-miRNA datasets are available. A pooled train/held-out rotation with these two datasets is exactly Experiments A and B, so no additional pooled result is reported as if it were an independent third rotation.

## Experiment D: Platform Shift

The distribution table compares log1p shared-miRNA values using medians and two-sample KS statistics. The raw, within-dataset, training-derived, and quantile modes are evaluated above. Within-dataset normalization uses unlabeled target feature distributions and is therefore explicitly labeled transductive; it does not use target labels.

The external cohort remains untouched. A high ranking score with a failed source-derived threshold is reported as threshold/platform failure, not hidden or converted into diagnostic success.

Machine-readable outputs: `reports/gse106817_phase5_results.json`, `reports/gse106817_phase5_predictions.csv`, and `reports/gse106817_phase5_distribution.csv`.
