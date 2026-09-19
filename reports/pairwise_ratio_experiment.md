# Pairwise Ratio Experiment

Frozen signature: `hsa-miR-200b-3p, hsa-miR-296-3p, hsa-miR-6826-5p, hsa-miR-32-3p, hsa-miR-1255a`.
Pairwise feature count: 10. A feature is 1 when the left miRNA expression exceeds the right miRNA expression, otherwise 0.
External cohorts are evaluated with a model fitted only on the local GSE94533 cohort.

| Cohort | Evaluation | ROC-AUC | Balanced Accuracy | Sensitivity | Specificity | Detail |
|---|---|---:|---:|---:|---:|---|
| GSE94533 | 5-fold stratified CV | 0.7101 | 0.7325 | 0.5962 | 0.8689 |  |
| GSE106817 | frozen local model | 0.5123 | 0.5025 | 0.8812 | 0.1237 |  |
| GSE113486 | frozen local model | 0.4741 | 0.4775 | 0.8750 | 0.0800 |  |
| GSE113740 | frozen local model | 0.4379 | 0.4548 | 0.8400 | 0.0697 |  |
