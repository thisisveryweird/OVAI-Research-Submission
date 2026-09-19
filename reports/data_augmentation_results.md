# Phase 9 Data Augmentation Results

This study uses only the verified local GSE94533 early-stage endpoint: 52 Stage I/II cancers versus 61 benign/control patients. Borderline tumors, HE4, and external datasets are excluded.

## Leakage Controls

Stratified 5-fold cross-validation was used. Top-5 miRNA selection, Age/CA-125 median imputation, and scaling were fit on each training fold only. SMOTE and BorderlineSMOTE were applied only after those training-fold transformations and never touched validation rows or validation labels.

Synthetic samples were generated only for the positive early-stage training class until the training classes were balanced. Validation rows remain the original biological patients.

## Results

| Model | Augmentation | ROC-AUC | PR-AUC | Accuracy | Balanced accuracy | Sensitivity | Specificity | Brier score | Geometric mean |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| random_forest | borderline_smote | 0.809 | 0.822 | 0.804 | 0.804 | 0.791 | 0.818 | 0.177 | 0.792 |
| random_forest | none | 0.806 | 0.828 | 0.787 | 0.788 | 0.791 | 0.786 | 0.175 | 0.774 |
| random_forest | smote | 0.809 | 0.825 | 0.805 | 0.804 | 0.791 | 0.818 | 0.177 | 0.791 |
| regularized_logistic | borderline_smote | 0.825 | 0.836 | 0.813 | 0.809 | 0.751 | 0.868 | 0.162 | 0.797 |
| regularized_logistic | none | 0.827 | 0.838 | 0.795 | 0.793 | 0.751 | 0.835 | 0.161 | 0.776 |
| regularized_logistic | smote | 0.821 | 0.831 | 0.786 | 0.784 | 0.751 | 0.818 | 0.163 | 0.770 |
| soft_voting | borderline_smote | 0.819 | 0.833 | 0.795 | 0.794 | 0.771 | 0.818 | 0.164 | 0.782 |
| soft_voting | none | 0.825 | 0.840 | 0.804 | 0.803 | 0.771 | 0.835 | 0.163 | 0.791 |
| soft_voting | smote | 0.817 | 0.828 | 0.804 | 0.803 | 0.771 | 0.835 | 0.164 | 0.792 |

## Conclusion

Synthetic augmentation did not demonstrate a clear sensitivity improvement without a material specificity tradeoff.
The largest sensitivity change was `random_forest` with `borderline_smote`: sensitivity delta +0.000, specificity delta +0.032 versus its no-augmentation baseline.

These are internal OOF research estimates. Synthetic augmentation must be re-tested on untouched external serum data before any generalization claim.

Machine-readable aggregate results: `reports/data_augmentation_results.csv`.
Machine-readable fold predictions: `reports/data_augmentation_predictions.csv`.
