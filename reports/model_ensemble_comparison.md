# Phase 8 Model and Ensemble Comparison

This study uses only verified local GSE94533 data. The endpoint is Stage I/II cancer versus benign/control; borderline tumors are excluded. HE4 and all external datasets are excluded.

## Leakage Controls

Stratified 5-fold outer evaluation was used. In every outer training fold, the Top-5 miRNA selector, median imputer, scaler, and inner-CV hyperparameter search were fit on training data only. The stacked logistic meta-model was trained only on inner out-of-fold base-model probabilities; no in-sample base predictions were used.

## Results

| Model | ROC-AUC | PR-AUC | Accuracy | Balanced accuracy | Sensitivity | Specificity | Brier score |
|---|---:|---:|---:|---:|---:|---:|---:|
| stacked_logistic | 0.825 | 0.841 | 0.760 | 0.766 | 0.829 | 0.704 | 0.184 |
| soft_voting | 0.823 | 0.843 | 0.804 | 0.804 | 0.791 | 0.818 | 0.159 |
| regularized_logistic | 0.823 | 0.834 | 0.813 | 0.811 | 0.771 | 0.851 | 0.162 |
| xgboost | 0.823 | 0.844 | 0.796 | 0.798 | 0.811 | 0.786 | 0.166 |
| calibrated_svm | 0.823 | 0.832 | 0.804 | 0.801 | 0.751 | 0.851 | 0.162 |
| random_forest | 0.812 | 0.829 | 0.804 | 0.806 | 0.811 | 0.801 | 0.178 |

## Validation Printout

`stacked_logistic`: meta-training used inner OOF base predictions only; in-sample base predictions were not used.

These are internal out-of-fold research estimates. No external validation claim is made.
