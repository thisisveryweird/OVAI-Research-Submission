# Phase 6 Feature-Selection Stability

This report was generated exclusively from the verified local GSE94533 cohort. No external dataset was loaded, and no feature was selected using held-out observations.

## Cohort and Protocol

The study contains 179 patients and 2578 measured miRNA features: 96 cancer, 46 benign, 22 borderline, and 15 controls.

Each task uses shuffled stratified 5-fold cross-validation with a fixed seed. In every training fold, median imputation, standardization, L1 or Elastic-Net selection, and the final classifier are fit using training patients only. The requested signature sizes are 5, 10, 14, 20, 30, and 50.

The early-stage endpoint is Stage I/II cancer versus benign/control. Borderline tumors are excluded from that endpoint and are not silently relabeled.

## Cross-Validated Metrics

Metrics below are means across held-out folds. They are descriptive OOF research estimates, not independent validation results.

| Task | Selector | Features | ROC-AUC | Balanced accuracy | Sensitivity | Specificity | Accuracy |
|---|---|---:|---:|---:|---:|---:|---:|
| all_cancer_vs_non_cancer | elastic_net | 5 | 0.626 | 0.617 | 0.561 | 0.673 | 0.614 |
| all_cancer_vs_non_cancer | elastic_net | 10 | 0.680 | 0.641 | 0.624 | 0.659 | 0.642 |
| all_cancer_vs_non_cancer | elastic_net | 14 | 0.698 | 0.641 | 0.634 | 0.649 | 0.642 |
| all_cancer_vs_non_cancer | elastic_net | 20 | 0.691 | 0.648 | 0.623 | 0.674 | 0.648 |
| all_cancer_vs_non_cancer | elastic_net | 30 | 0.656 | 0.629 | 0.707 | 0.551 | 0.637 |
| all_cancer_vs_non_cancer | elastic_net | 50 | 0.623 | 0.616 | 0.666 | 0.565 | 0.620 |
| all_cancer_vs_non_cancer | l1 | 5 | 0.706 | 0.685 | 0.614 | 0.757 | 0.681 |
| all_cancer_vs_non_cancer | l1 | 10 | 0.716 | 0.682 | 0.644 | 0.720 | 0.681 |
| all_cancer_vs_non_cancer | l1 | 14 | 0.687 | 0.676 | 0.655 | 0.696 | 0.676 |
| all_cancer_vs_non_cancer | l1 | 20 | 0.642 | 0.596 | 0.603 | 0.588 | 0.598 |
| all_cancer_vs_non_cancer | l1 | 30 | 0.664 | 0.635 | 0.623 | 0.648 | 0.637 |
| all_cancer_vs_non_cancer | l1 | 50 | 0.659 | 0.645 | 0.654 | 0.635 | 0.648 |
| early_stage_vs_benign_control | elastic_net | 5 | 0.639 | 0.648 | 0.476 | 0.819 | 0.660 |
| early_stage_vs_benign_control | elastic_net | 10 | 0.670 | 0.676 | 0.633 | 0.719 | 0.679 |
| early_stage_vs_benign_control | elastic_net | 14 | 0.669 | 0.647 | 0.575 | 0.719 | 0.653 |
| early_stage_vs_benign_control | elastic_net | 20 | 0.698 | 0.642 | 0.613 | 0.672 | 0.645 |
| early_stage_vs_benign_control | elastic_net | 30 | 0.682 | 0.650 | 0.629 | 0.672 | 0.654 |
| early_stage_vs_benign_control | elastic_net | 50 | 0.694 | 0.622 | 0.556 | 0.688 | 0.628 |
| early_stage_vs_benign_control | l1 | 5 | 0.749 | 0.723 | 0.595 | 0.851 | 0.733 |
| early_stage_vs_benign_control | l1 | 10 | 0.726 | 0.709 | 0.616 | 0.801 | 0.716 |
| early_stage_vs_benign_control | l1 | 14 | 0.728 | 0.695 | 0.638 | 0.753 | 0.699 |
| early_stage_vs_benign_control | l1 | 20 | 0.724 | 0.686 | 0.638 | 0.735 | 0.691 |
| early_stage_vs_benign_control | l1 | 30 | 0.739 | 0.634 | 0.598 | 0.669 | 0.638 |
| early_stage_vs_benign_control | l1 | 50 | 0.724 | 0.650 | 0.598 | 0.703 | 0.655 |

## Feature Stability

Selection frequency is the number of outer folds in which a feature was selected divided by five. High frequency indicates reproducibility under these folds; it is not proof of biological causality or external generalization.

### all_cancer_vs_non_cancer

#### l1

- Top 5: hsa-miR-200c-3p (5/5); hsa-miR-32-3p (5/5); hsa-miR-708-3p (3/5); hsa-miR-1538 (2/5); hsa-miR-1255a (1/5); hsa-miR-1288-3p (1/5); hsa-miR-30d-5p (1/5); hsa-miR-3150a-5p (1/5); hsa-miR-380-5p (1/5); hsa-miR-4638-5p (1/5); hsa-miR-495-3p (1/5); hsa-miR-548ah-3p (1/5); hsa-miR-671-5p (1/5); hsa-miR-6891-5p (1/5)
- Top 10: hsa-miR-200c-3p (5/5); hsa-miR-32-3p (5/5); hsa-miR-1538 (4/5); hsa-miR-671-5p (4/5); hsa-miR-380-5p (3/5); hsa-miR-708-3p (3/5); hsa-miR-933 (3/5); hsa-miR-200a-5p (2/5); hsa-miR-30d-5p (2/5); hsa-miR-4638-5p (2/5); hsa-miR-495-3p (2/5); hsa-miR-10a-3p (1/5); hsa-miR-1255a (1/5); hsa-miR-1288-3p (1/5); hsa-miR-1291 (1/5)
- Top 14: hsa-miR-200c-3p (5/5); hsa-miR-32-3p (5/5); hsa-miR-671-5p (5/5); hsa-miR-1538 (4/5); hsa-miR-380-5p (4/5); hsa-miR-200a-5p (3/5); hsa-miR-708-3p (3/5); hsa-miR-933 (3/5); hsa-miR-10a-3p (2/5); hsa-miR-1291 (2/5); hsa-miR-30d-5p (2/5); hsa-miR-4638-5p (2/5); hsa-miR-495-3p (2/5); hsa-miR-1255a (1/5); hsa-miR-1288-3p (1/5)
- Top 20: hsa-miR-1538 (5/5); hsa-miR-200c-3p (5/5); hsa-miR-32-3p (5/5); hsa-miR-671-5p (5/5); hsa-miR-380-5p (4/5); hsa-miR-10a-3p (3/5); hsa-miR-200a-5p (3/5); hsa-miR-3150a-5p (3/5); hsa-miR-4638-5p (3/5); hsa-miR-708-3p (3/5); hsa-miR-933 (3/5); hsa-let-7a-2-3p (2/5); hsa-miR-1291 (2/5); hsa-miR-30d-5p (2/5); hsa-miR-320c (2/5)
- Top 30: hsa-miR-1538 (5/5); hsa-miR-200c-3p (5/5); hsa-miR-32-3p (5/5); hsa-miR-671-5p (5/5); hsa-let-7a-2-3p (4/5); hsa-let-7a-3p (4/5); hsa-let-7a-5p (4/5); hsa-let-7b-3p (4/5); hsa-let-7b-5p (4/5); hsa-let-7c-3p (4/5); hsa-miR-380-5p (4/5); hsa-miR-708-3p (4/5); hsa-miR-933 (4/5); hsa-let-7c-5p (3/5); hsa-let-7d-3p (3/5)
- Top 50: hsa-let-7a-2-3p (5/5); hsa-let-7a-3p (5/5); hsa-let-7a-5p (5/5); hsa-let-7b-3p (5/5); hsa-let-7b-5p (5/5); hsa-let-7c-3p (5/5); hsa-let-7c-5p (5/5); hsa-let-7d-3p (5/5); hsa-let-7d-5p (5/5); hsa-let-7e-3p (5/5); hsa-let-7e-5p (5/5); hsa-let-7f-1-3p (5/5); hsa-let-7f-2-3p (5/5); hsa-let-7f-5p (5/5); hsa-let-7g-3p (5/5)

#### elastic_net

- Top 5: hsa-miR-32-3p (5/5); hsa-miR-200c-3p (4/5); hsa-miR-708-3p (3/5); hsa-miR-4638-5p (2/5); hsa-miR-933 (2/5); hsa-miR-1255a (1/5); hsa-miR-1288-3p (1/5); hsa-miR-1538 (1/5); hsa-miR-4661-5p (1/5); hsa-miR-4685-3p (1/5); hsa-miR-4754 (1/5); hsa-miR-548ah-3p (1/5); hsa-miR-605-3p (1/5); hsa-miR-6891-5p (1/5)
- Top 10: hsa-miR-200c-3p (5/5); hsa-miR-32-3p (5/5); hsa-miR-1255a (3/5); hsa-miR-1538 (3/5); hsa-miR-708-3p (3/5); hsa-miR-4638-5p (2/5); hsa-miR-4661-5p (2/5); hsa-miR-4754 (2/5); hsa-miR-671-5p (2/5); hsa-miR-6826-5p (2/5); hsa-miR-933 (2/5); hsa-miR-1288-3p (1/5); hsa-miR-1291 (1/5); hsa-miR-1307-5p (1/5); hsa-miR-200a-5p (1/5)
- Top 14: hsa-miR-200c-3p (5/5); hsa-miR-32-3p (5/5); hsa-miR-1255a (4/5); hsa-miR-1538 (3/5); hsa-miR-380-5p (3/5); hsa-miR-671-5p (3/5); hsa-miR-6826-5p (3/5); hsa-miR-708-3p (3/5); hsa-miR-933 (3/5); hsa-miR-1307-5p (2/5); hsa-miR-200a-5p (2/5); hsa-miR-4638-5p (2/5); hsa-miR-4661-5p (2/5); hsa-miR-4754 (2/5); hsa-miR-495-3p (2/5)
- Top 20: hsa-miR-1255a (5/5); hsa-miR-1538 (5/5); hsa-miR-200c-3p (5/5); hsa-miR-32-3p (5/5); hsa-miR-380-5p (5/5); hsa-miR-708-3p (5/5); hsa-miR-933 (4/5); hsa-miR-1307-5p (3/5); hsa-miR-200a-5p (3/5); hsa-miR-4638-5p (3/5); hsa-miR-4754 (3/5); hsa-miR-671-5p (3/5); hsa-miR-6826-5p (3/5); hsa-miR-6853-3p (3/5); hsa-miR-378i (2/5)
- Top 30: hsa-miR-1255a (5/5); hsa-miR-1538 (5/5); hsa-miR-200c-3p (5/5); hsa-miR-32-3p (5/5); hsa-miR-380-5p (5/5); hsa-miR-708-3p (5/5); hsa-miR-671-5p (4/5); hsa-miR-933 (4/5); hsa-miR-1307-5p (3/5); hsa-miR-200a-5p (3/5); hsa-miR-21-3p (3/5); hsa-miR-378i (3/5); hsa-miR-4638-5p (3/5); hsa-miR-4661-5p (3/5); hsa-miR-4754 (3/5)
- Top 50: hsa-miR-1255a (5/5); hsa-miR-1538 (5/5); hsa-miR-200c-3p (5/5); hsa-miR-21-3p (5/5); hsa-miR-32-3p (5/5); hsa-miR-380-5p (5/5); hsa-miR-671-5p (5/5); hsa-miR-708-3p (5/5); hsa-miR-933 (5/5); hsa-miR-200a-5p (4/5); hsa-miR-378i (4/5); hsa-miR-495-3p (4/5); hsa-miR-6826-5p (4/5); hsa-let-7b-3p (3/5); hsa-miR-10a-3p (3/5)

### early_stage_vs_benign_control

#### l1

- Top 5: hsa-miR-296-3p (5/5); hsa-miR-1255a (4/5); hsa-miR-200b-3p (4/5); hsa-miR-200a-3p (3/5); hsa-miR-32-3p (3/5); hsa-miR-196a-5p (1/5); hsa-miR-200a-5p (1/5); hsa-miR-3149 (1/5); hsa-miR-4668-5p (1/5); hsa-miR-495-3p (1/5); hsa-miR-6826-5p (1/5)
- Top 10: hsa-miR-1255a (5/5); hsa-miR-200b-3p (5/5); hsa-miR-296-3p (5/5); hsa-let-7a-2-3p (4/5); hsa-let-7a-3p (3/5); hsa-miR-200a-3p (3/5); hsa-miR-32-3p (3/5); hsa-miR-6826-5p (3/5); hsa-let-7a-5p (2/5); hsa-miR-1538 (2/5); hsa-miR-92a-1-5p (2/5); hsa-let-7b-3p (1/5); hsa-let-7f-1-3p (1/5); hsa-miR-1299 (1/5); hsa-miR-196a-5p (1/5)
- Top 14: hsa-let-7a-2-3p (5/5); hsa-let-7a-3p (5/5); hsa-let-7a-5p (5/5); hsa-miR-1255a (5/5); hsa-miR-200b-3p (5/5); hsa-miR-296-3p (5/5); hsa-let-7b-3p (4/5); hsa-let-7b-5p (4/5); hsa-let-7c-3p (3/5); hsa-miR-200a-3p (3/5); hsa-miR-32-3p (3/5); hsa-miR-6826-5p (3/5); hsa-let-7c-5p (2/5); hsa-miR-1538 (2/5); hsa-miR-92a-1-5p (2/5)
- Top 20: hsa-let-7a-2-3p (5/5); hsa-let-7a-3p (5/5); hsa-let-7a-5p (5/5); hsa-let-7b-3p (5/5); hsa-let-7b-5p (5/5); hsa-let-7c-3p (5/5); hsa-let-7c-5p (5/5); hsa-let-7d-3p (5/5); hsa-let-7d-5p (5/5); hsa-miR-1255a (5/5); hsa-miR-200b-3p (5/5); hsa-miR-296-3p (5/5); hsa-let-7e-3p (4/5); hsa-let-7e-5p (4/5); hsa-let-7f-1-3p (4/5)
- Top 30: hsa-let-7a-2-3p (5/5); hsa-let-7a-3p (5/5); hsa-let-7a-5p (5/5); hsa-let-7b-3p (5/5); hsa-let-7b-5p (5/5); hsa-let-7c-3p (5/5); hsa-let-7c-5p (5/5); hsa-let-7d-3p (5/5); hsa-let-7d-5p (5/5); hsa-let-7e-3p (5/5); hsa-let-7e-5p (5/5); hsa-let-7f-1-3p (5/5); hsa-let-7f-2-3p (5/5); hsa-let-7f-5p (5/5); hsa-let-7g-3p (5/5)
- Top 50: hsa-let-7a-2-3p (5/5); hsa-let-7a-3p (5/5); hsa-let-7a-5p (5/5); hsa-let-7b-3p (5/5); hsa-let-7b-5p (5/5); hsa-let-7c-3p (5/5); hsa-let-7c-5p (5/5); hsa-let-7d-3p (5/5); hsa-let-7d-5p (5/5); hsa-let-7e-3p (5/5); hsa-let-7e-5p (5/5); hsa-let-7f-1-3p (5/5); hsa-let-7f-2-3p (5/5); hsa-let-7f-5p (5/5); hsa-let-7g-3p (5/5)

#### elastic_net

- Top 5: hsa-miR-296-3p (5/5); hsa-miR-200b-3p (3/5); hsa-miR-6826-5p (3/5); hsa-miR-1255a (2/5); hsa-miR-1538 (2/5); hsa-miR-1291 (1/5); hsa-miR-1299 (1/5); hsa-miR-1303 (1/5); hsa-miR-196a-5p (1/5); hsa-miR-200a-3p (1/5); hsa-miR-200a-5p (1/5); hsa-miR-3149 (1/5); hsa-miR-32-3p (1/5); hsa-miR-4742-5p (1/5); hsa-miR-708-3p (1/5)
- Top 10: hsa-miR-296-3p (5/5); hsa-miR-1255a (4/5); hsa-miR-200b-3p (4/5); hsa-miR-32-3p (4/5); hsa-miR-6826-5p (4/5); hsa-miR-1538 (3/5); hsa-miR-200a-3p (3/5); hsa-miR-7847-3p (3/5); hsa-miR-4742-5p (2/5); hsa-miR-92a-1-5p (2/5); hsa-let-7f-1-3p (1/5); hsa-miR-1288-3p (1/5); hsa-miR-1291 (1/5); hsa-miR-1299 (1/5); hsa-miR-1303 (1/5)
- Top 14: hsa-miR-1255a (5/5); hsa-miR-296-3p (5/5); hsa-miR-32-3p (5/5); hsa-miR-6826-5p (5/5); hsa-miR-200b-3p (4/5); hsa-miR-1538 (3/5); hsa-miR-200a-3p (3/5); hsa-miR-4742-5p (3/5); hsa-miR-7847-3p (3/5); hsa-miR-708-3p (2/5); hsa-miR-92a-1-5p (2/5); hsa-let-7f-1-3p (1/5); hsa-miR-1288-3p (1/5); hsa-miR-1291 (1/5); hsa-miR-1299 (1/5)
- Top 20: hsa-miR-1255a (5/5); hsa-miR-1538 (5/5); hsa-miR-200b-3p (5/5); hsa-miR-296-3p (5/5); hsa-miR-32-3p (5/5); hsa-miR-6826-5p (5/5); hsa-miR-200a-3p (3/5); hsa-miR-4742-5p (3/5); hsa-miR-4999-5p (3/5); hsa-miR-6503-3p (3/5); hsa-miR-7847-3p (3/5); hsa-miR-92a-1-5p (3/5); hsa-miR-3149 (2/5); hsa-miR-380-5p (2/5); hsa-miR-548ah-3p (2/5)
- Top 30: hsa-miR-1255a (5/5); hsa-miR-1538 (5/5); hsa-miR-200b-3p (5/5); hsa-miR-296-3p (5/5); hsa-miR-32-3p (5/5); hsa-miR-6826-5p (5/5); hsa-miR-4742-5p (4/5); hsa-miR-4999-5p (4/5); hsa-miR-200a-3p (3/5); hsa-miR-200a-5p (3/5); hsa-miR-200c-3p (3/5); hsa-miR-3149 (3/5); hsa-miR-4754 (3/5); hsa-miR-491-5p (3/5); hsa-miR-6503-3p (3/5)
- Top 50: hsa-let-7a-2-3p (5/5); hsa-let-7a-3p (5/5); hsa-miR-1255a (5/5); hsa-miR-1538 (5/5); hsa-miR-196a-5p (5/5); hsa-miR-200b-3p (5/5); hsa-miR-296-3p (5/5); hsa-miR-3149 (5/5); hsa-miR-32-3p (5/5); hsa-miR-6826-5p (5/5); hsa-let-7a-5p (4/5); hsa-miR-200a-3p (4/5); hsa-miR-200a-5p (4/5); hsa-miR-3135a (4/5); hsa-miR-4742-5p (4/5)

## Interpretation Guardrails

- The 14-feature signature size is included as a sparsity point, not privileged as the historical signature.
- Feature stability is estimated within GSE94533 only. It must not be interpreted as stability across GSE106817, CPTAC, TCGA, or other platforms.
- Model and signature-size comparisons should consider ROC-AUC, balanced accuracy, sensitivity, specificity, and fold variability together; accuracy alone is insufficient.
- Borderline tumors remain represented in the all-cancer-versus-non-cancer endpoint but are explicitly excluded from the early-stage-versus-benign/control endpoint.

Machine-readable fold summaries: `reports/feature_selection_results.csv`.
Machine-readable stability table: `reports/feature_selection_stability.csv`.
