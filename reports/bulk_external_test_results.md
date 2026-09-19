# Bulk External GEO Dataset Test Results

Evaluation uses the complete local GSE94533 Stage I/II cancer versus benign/control cohort. The model is a regularized logistic regression with local median imputation and scaling; external labels are never used for fitting or tuning.

- Candidates tested: 10
- Frozen five-miRNA signature: `hsa-miR-200b-3p, hsa-miR-296-3p, hsa-miR-6826-5p, hsa-miR-32-3p, hsa-miR-1255a`
- Success rule: sensitivity > 50% and balanced accuracy > 60%.

| GSE | Status | Recommendation | N evaluated | ROC-AUC | Balanced Accuracy | Sensitivity | Specificity | Detail |
|---|---|---|---:|---:|---:|---:|---:|---|
| GSE106817 | FAILED | Failed: Platform Shift / 0% Sensitivity | 3108 | 0.838 | 0.500 | 1.000 | 0.000 |  |
| GSE113486 | FAILED | Failed: Platform Shift / 0% Sensitivity | 140 | 0.830 | 0.500 | 1.000 | 0.000 |  |
| GSE113740 | FAILED | Failed: Platform Shift / 0% Sensitivity | 1058 | 0.769 | 0.500 | 1.000 | 0.000 |  |
| GSE19711 | Skipped: Missing Features | Skipped: Missing Features | - | - | - | - | - | Missing: hsa-miR-200b-3p, hsa-miR-296-3p, hsa-miR-6826-5p, hsa-miR-32-3p, hsa-miR-1255a |
| GSE226445 | Skipped: Missing Features | Skipped: Missing Features | - | - | - | - | - | Missing: hsa-miR-200b-3p, hsa-miR-296-3p, hsa-miR-6826-5p, hsa-miR-32-3p, hsa-miR-1255a |
| GSE235525 | Skipped: Missing Features | Skipped: Missing Features | - | - | - | - | - | Missing: hsa-miR-6826-5p, hsa-miR-32-3p |
| GSE300882 | Skipped: Missing Features | Skipped: Missing Features | - | - | - | - | - | Missing: hsa-miR-200b-3p, hsa-miR-296-3p, hsa-miR-6826-5p, hsa-miR-32-3p, hsa-miR-1255a |
| GSE61741 | Skipped: Missing Features | Skipped: Missing Features | - | - | - | - | - | Missing: hsa-miR-200b-3p, hsa-miR-6826-5p, hsa-miR-32-3p |
| GSE94533 | Skipped: Missing Features | Skipped: Missing Features | - | - | - | - | - | Missing: hsa-miR-200b-3p, hsa-miR-296-3p, hsa-miR-6826-5p, hsa-miR-32-3p, hsa-miR-1255a |
| GSE138092 | Failed: Missing Labels | Failed: Missing Labels | - | - | - | - | - | Need both ovarian cancer and non-cancer labeled samples |
