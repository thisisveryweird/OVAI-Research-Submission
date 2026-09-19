# Phase 7 Clinical Ablation Results

This study uses only the verified local GSE94533 cohort. HE4 is not present and is excluded from every pipeline. No external dataset is loaded.

## Protocol

Each configuration uses shuffled stratified 5-fold cross-validation. Age is median-imputed inside each training fold. CA-125 is evaluated with either fold-local median imputation or fold-local KNN imputation. The Top-5 miRNA signature is selected by fold-local L1 logistic regression from the full local miRNA feature space; no global five-feature list is used.

Borderline tumors are included in the all-cancer versus non-cancer endpoint and excluded from the Stage I/II versus benign/control endpoint.

## Results

| Endpoint | Pipeline | Model | Imputation | ROC-AUC | Balanced accuracy | Sensitivity | Specificity | Brier score |
|---|---|---|---|---:|---:|---:|---:|---:|
| all_cancer_vs_non_cancer | clinical_only | calibrated_svm | knn_ca125 | 0.642 | 0.595 | 0.686 | 0.504 | 0.235 |
| all_cancer_vs_non_cancer | clinical_only | calibrated_svm | median | 0.688 | 0.665 | 0.656 | 0.674 | 0.221 |
| all_cancer_vs_non_cancer | clinical_only | regularized_logistic | knn_ca125 | 0.648 | 0.596 | 0.552 | 0.639 | 0.235 |
| all_cancer_vs_non_cancer | clinical_only | regularized_logistic | median | 0.684 | 0.635 | 0.512 | 0.759 | 0.225 |
| all_cancer_vs_non_cancer | combined | calibrated_svm | knn_ca125 | 0.748 | 0.693 | 0.738 | 0.647 | 0.214 |
| all_cancer_vs_non_cancer | combined | calibrated_svm | median | 0.750 | 0.705 | 0.728 | 0.682 | 0.211 |
| all_cancer_vs_non_cancer | combined | regularized_logistic | knn_ca125 | 0.753 | 0.701 | 0.708 | 0.694 | 0.211 |
| all_cancer_vs_non_cancer | combined | regularized_logistic | median | 0.762 | 0.703 | 0.687 | 0.718 | 0.207 |
| all_cancer_vs_non_cancer | mirna_only | calibrated_svm | knn_ca125 | 0.736 | 0.678 | 0.718 | 0.637 | 0.215 |
| all_cancer_vs_non_cancer | mirna_only | calibrated_svm | median | 0.736 | 0.678 | 0.718 | 0.637 | 0.215 |
| all_cancer_vs_non_cancer | mirna_only | regularized_logistic | knn_ca125 | 0.743 | 0.680 | 0.676 | 0.685 | 0.211 |
| all_cancer_vs_non_cancer | mirna_only | regularized_logistic | median | 0.743 | 0.680 | 0.676 | 0.685 | 0.211 |
| early_stage_vs_benign_control | clinical_only | calibrated_svm | knn_ca125 | 0.565 | 0.530 | 0.211 | 0.850 | 0.246 |
| early_stage_vs_benign_control | clinical_only | calibrated_svm | median | 0.665 | 0.601 | 0.349 | 0.853 | 0.235 |
| early_stage_vs_benign_control | clinical_only | regularized_logistic | knn_ca125 | 0.586 | 0.567 | 0.444 | 0.690 | 0.248 |
| early_stage_vs_benign_control | clinical_only | regularized_logistic | median | 0.642 | 0.601 | 0.464 | 0.738 | 0.241 |
| early_stage_vs_benign_control | combined | calibrated_svm | knn_ca125 | 0.813 | 0.774 | 0.713 | 0.836 | 0.162 |
| early_stage_vs_benign_control | combined | calibrated_svm | median | 0.821 | 0.801 | 0.751 | 0.851 | 0.162 |
| early_stage_vs_benign_control | combined | regularized_logistic | knn_ca125 | 0.826 | 0.801 | 0.751 | 0.851 | 0.161 |
| early_stage_vs_benign_control | combined | regularized_logistic | median | 0.827 | 0.793 | 0.751 | 0.835 | 0.161 |
| early_stage_vs_benign_control | mirna_only | calibrated_svm | knn_ca125 | 0.808 | 0.773 | 0.711 | 0.835 | 0.164 |
| early_stage_vs_benign_control | mirna_only | calibrated_svm | median | 0.808 | 0.773 | 0.711 | 0.835 | 0.164 |
| early_stage_vs_benign_control | mirna_only | regularized_logistic | knn_ca125 | 0.815 | 0.766 | 0.731 | 0.801 | 0.165 |
| early_stage_vs_benign_control | mirna_only | regularized_logistic | median | 0.815 | 0.766 | 0.731 | 0.801 | 0.165 |

## Top-5 Selection Stability

Selection frequency is measured across the five folds for each configuration. It is a within-cohort stability measure, not evidence of biological causality or external generalization.

### all_cancer_vs_non_cancer

- mirna_only/regularized_logistic/median: hsa-miR-200c-3p (5/5); hsa-miR-671-5p (4/5); hsa-miR-30d-5p (3/5); hsa-miR-32-3p (3/5); hsa-miR-200b-3p (2/5); hsa-miR-10a-3p (1/5); hsa-miR-1255a (1/5); hsa-miR-1288-3p (1/5); hsa-miR-296-3p (1/5); hsa-miR-3150a-5p (1/5)
- mirna_only/regularized_logistic/knn_ca125: hsa-miR-200c-3p (5/5); hsa-miR-671-5p (4/5); hsa-miR-30d-5p (3/5); hsa-miR-32-3p (3/5); hsa-miR-200b-3p (2/5); hsa-miR-10a-3p (1/5); hsa-miR-1255a (1/5); hsa-miR-1288-3p (1/5); hsa-miR-296-3p (1/5); hsa-miR-3150a-5p (1/5)
- mirna_only/calibrated_svm/median: hsa-miR-200c-3p (5/5); hsa-miR-671-5p (4/5); hsa-miR-30d-5p (3/5); hsa-miR-32-3p (3/5); hsa-miR-200b-3p (2/5); hsa-miR-10a-3p (1/5); hsa-miR-1255a (1/5); hsa-miR-1288-3p (1/5); hsa-miR-296-3p (1/5); hsa-miR-3150a-5p (1/5)
- mirna_only/calibrated_svm/knn_ca125: hsa-miR-200c-3p (5/5); hsa-miR-671-5p (4/5); hsa-miR-30d-5p (3/5); hsa-miR-32-3p (3/5); hsa-miR-200b-3p (2/5); hsa-miR-10a-3p (1/5); hsa-miR-1255a (1/5); hsa-miR-1288-3p (1/5); hsa-miR-296-3p (1/5); hsa-miR-3150a-5p (1/5)
- clinical_only/regularized_logistic/median: hsa-miR-200c-3p (5/5); hsa-miR-671-5p (4/5); hsa-miR-30d-5p (3/5); hsa-miR-32-3p (3/5); hsa-miR-200b-3p (2/5); hsa-miR-10a-3p (1/5); hsa-miR-1255a (1/5); hsa-miR-1288-3p (1/5); hsa-miR-296-3p (1/5); hsa-miR-3150a-5p (1/5)
- clinical_only/regularized_logistic/knn_ca125: hsa-miR-200c-3p (5/5); hsa-miR-671-5p (4/5); hsa-miR-30d-5p (3/5); hsa-miR-32-3p (3/5); hsa-miR-200b-3p (2/5); hsa-miR-10a-3p (1/5); hsa-miR-1255a (1/5); hsa-miR-1288-3p (1/5); hsa-miR-296-3p (1/5); hsa-miR-3150a-5p (1/5)
- clinical_only/calibrated_svm/median: hsa-miR-200c-3p (5/5); hsa-miR-671-5p (4/5); hsa-miR-30d-5p (3/5); hsa-miR-32-3p (3/5); hsa-miR-200b-3p (2/5); hsa-miR-10a-3p (1/5); hsa-miR-1255a (1/5); hsa-miR-1288-3p (1/5); hsa-miR-296-3p (1/5); hsa-miR-3150a-5p (1/5)
- clinical_only/calibrated_svm/knn_ca125: hsa-miR-200c-3p (5/5); hsa-miR-671-5p (4/5); hsa-miR-30d-5p (3/5); hsa-miR-32-3p (3/5); hsa-miR-200b-3p (2/5); hsa-miR-10a-3p (1/5); hsa-miR-1255a (1/5); hsa-miR-1288-3p (1/5); hsa-miR-296-3p (1/5); hsa-miR-3150a-5p (1/5)
- combined/regularized_logistic/median: hsa-miR-200c-3p (5/5); hsa-miR-671-5p (4/5); hsa-miR-30d-5p (3/5); hsa-miR-32-3p (3/5); hsa-miR-200b-3p (2/5); hsa-miR-10a-3p (1/5); hsa-miR-1255a (1/5); hsa-miR-1288-3p (1/5); hsa-miR-296-3p (1/5); hsa-miR-3150a-5p (1/5)
- combined/regularized_logistic/knn_ca125: hsa-miR-200c-3p (5/5); hsa-miR-671-5p (4/5); hsa-miR-30d-5p (3/5); hsa-miR-32-3p (3/5); hsa-miR-200b-3p (2/5); hsa-miR-10a-3p (1/5); hsa-miR-1255a (1/5); hsa-miR-1288-3p (1/5); hsa-miR-296-3p (1/5); hsa-miR-3150a-5p (1/5)
- combined/calibrated_svm/median: hsa-miR-200c-3p (5/5); hsa-miR-671-5p (4/5); hsa-miR-30d-5p (3/5); hsa-miR-32-3p (3/5); hsa-miR-200b-3p (2/5); hsa-miR-10a-3p (1/5); hsa-miR-1255a (1/5); hsa-miR-1288-3p (1/5); hsa-miR-296-3p (1/5); hsa-miR-3150a-5p (1/5)
- combined/calibrated_svm/knn_ca125: hsa-miR-200c-3p (5/5); hsa-miR-671-5p (4/5); hsa-miR-30d-5p (3/5); hsa-miR-32-3p (3/5); hsa-miR-200b-3p (2/5); hsa-miR-10a-3p (1/5); hsa-miR-1255a (1/5); hsa-miR-1288-3p (1/5); hsa-miR-296-3p (1/5); hsa-miR-3150a-5p (1/5)

### early_stage_vs_benign_control

- mirna_only/regularized_logistic/median: hsa-miR-200b-3p (5/5); hsa-miR-296-3p (5/5); hsa-miR-6826-5p (4/5); hsa-miR-1255a (3/5); hsa-let-7a-2-3p (2/5); hsa-let-7a-3p (1/5); hsa-miR-1303 (1/5); hsa-miR-200c-3p (1/5); hsa-miR-3149 (1/5); hsa-miR-32-3p (1/5)
- mirna_only/regularized_logistic/knn_ca125: hsa-miR-200b-3p (5/5); hsa-miR-296-3p (5/5); hsa-miR-6826-5p (4/5); hsa-miR-1255a (3/5); hsa-let-7a-2-3p (2/5); hsa-let-7a-3p (1/5); hsa-miR-1303 (1/5); hsa-miR-200c-3p (1/5); hsa-miR-3149 (1/5); hsa-miR-32-3p (1/5)
- mirna_only/calibrated_svm/median: hsa-miR-200b-3p (5/5); hsa-miR-296-3p (5/5); hsa-miR-6826-5p (4/5); hsa-miR-1255a (3/5); hsa-let-7a-2-3p (2/5); hsa-let-7a-3p (1/5); hsa-miR-1303 (1/5); hsa-miR-200c-3p (1/5); hsa-miR-3149 (1/5); hsa-miR-32-3p (1/5)
- mirna_only/calibrated_svm/knn_ca125: hsa-miR-200b-3p (5/5); hsa-miR-296-3p (5/5); hsa-miR-6826-5p (4/5); hsa-miR-1255a (3/5); hsa-let-7a-2-3p (2/5); hsa-let-7a-3p (1/5); hsa-miR-1303 (1/5); hsa-miR-200c-3p (1/5); hsa-miR-3149 (1/5); hsa-miR-32-3p (1/5)
- clinical_only/regularized_logistic/median: hsa-miR-200b-3p (5/5); hsa-miR-296-3p (5/5); hsa-miR-6826-5p (4/5); hsa-miR-1255a (3/5); hsa-let-7a-2-3p (2/5); hsa-let-7a-3p (1/5); hsa-miR-1303 (1/5); hsa-miR-200c-3p (1/5); hsa-miR-3149 (1/5); hsa-miR-32-3p (1/5)
- clinical_only/regularized_logistic/knn_ca125: hsa-miR-200b-3p (5/5); hsa-miR-296-3p (5/5); hsa-miR-6826-5p (4/5); hsa-miR-1255a (3/5); hsa-let-7a-2-3p (2/5); hsa-let-7a-3p (1/5); hsa-miR-1303 (1/5); hsa-miR-200c-3p (1/5); hsa-miR-3149 (1/5); hsa-miR-32-3p (1/5)
- clinical_only/calibrated_svm/median: hsa-miR-200b-3p (5/5); hsa-miR-296-3p (5/5); hsa-miR-6826-5p (4/5); hsa-miR-1255a (3/5); hsa-let-7a-2-3p (2/5); hsa-let-7a-3p (1/5); hsa-miR-1303 (1/5); hsa-miR-200c-3p (1/5); hsa-miR-3149 (1/5); hsa-miR-32-3p (1/5)
- clinical_only/calibrated_svm/knn_ca125: hsa-miR-200b-3p (5/5); hsa-miR-296-3p (5/5); hsa-miR-6826-5p (4/5); hsa-miR-1255a (3/5); hsa-let-7a-2-3p (2/5); hsa-let-7a-3p (1/5); hsa-miR-1303 (1/5); hsa-miR-200c-3p (1/5); hsa-miR-3149 (1/5); hsa-miR-32-3p (1/5)
- combined/regularized_logistic/median: hsa-miR-200b-3p (5/5); hsa-miR-296-3p (5/5); hsa-miR-6826-5p (4/5); hsa-miR-1255a (3/5); hsa-let-7a-2-3p (2/5); hsa-let-7a-3p (1/5); hsa-miR-1303 (1/5); hsa-miR-200c-3p (1/5); hsa-miR-3149 (1/5); hsa-miR-32-3p (1/5)
- combined/regularized_logistic/knn_ca125: hsa-miR-200b-3p (5/5); hsa-miR-296-3p (5/5); hsa-miR-6826-5p (4/5); hsa-miR-1255a (3/5); hsa-let-7a-2-3p (2/5); hsa-let-7a-3p (1/5); hsa-miR-1303 (1/5); hsa-miR-200c-3p (1/5); hsa-miR-3149 (1/5); hsa-miR-32-3p (1/5)
- combined/calibrated_svm/median: hsa-miR-200b-3p (5/5); hsa-miR-296-3p (5/5); hsa-miR-6826-5p (4/5); hsa-miR-1255a (3/5); hsa-let-7a-2-3p (2/5); hsa-let-7a-3p (1/5); hsa-miR-1303 (1/5); hsa-miR-200c-3p (1/5); hsa-miR-3149 (1/5); hsa-miR-32-3p (1/5)
- combined/calibrated_svm/knn_ca125: hsa-miR-200b-3p (5/5); hsa-miR-296-3p (5/5); hsa-miR-6826-5p (4/5); hsa-miR-1255a (3/5); hsa-let-7a-2-3p (2/5); hsa-let-7a-3p (1/5); hsa-miR-1303 (1/5); hsa-miR-200c-3p (1/5); hsa-miR-3149 (1/5); hsa-miR-32-3p (1/5)

## Interpretation Guardrails

- A clinical-only result is not a molecular diagnostic result; it is an ablation comparator.
- Any combined-model improvement is an internal OOF estimate and requires independent validation before being considered robust.
- CA-125 has 60 missing values of 179; KNN imputation is therefore a sensitivity analysis, not a substitute for measured values.
- HE4 is unavailable and was excluded explicitly.

Machine-readable aggregate results: `reports/clinical_ablation_results.csv`.
Machine-readable Top-5 stability: `reports/clinical_ablation_stability.csv`.
