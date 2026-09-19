# Within-Patient Rank-Transform Experiment

Frozen signature: `hsa-miR-200b-3p, hsa-miR-296-3p, hsa-miR-6826-5p, hsa-miR-32-3p, hsa-miR-1255a`.
Ranks are computed independently within each sample against that sample's full available hsa-miRNA background. Age and CA-125 are excluded.

| Cohort | Space | ROC-AUC | Balanced Accuracy | Sensitivity | Specificity | Detail |
|---|---|---:|---:|---:|---:|---|
| GSE94533 (5-fold CV) | raw | 0.8507 | 0.7859 | 0.6538 | 0.9180 |  |
| GSE94533 (5-fold CV) | rank | 0.8856 | 0.8287 | 0.7885 | 0.8689 |  |
| GSE106817 | raw | 0.9135 | 0.5000 | 1.0000 | 0.0000 |  |
| GSE106817 | rank | 0.6061 | 0.5116 | 0.0281 | 0.9950 |  |
| GSE113486 | raw | 0.9468 | 0.5000 | 1.0000 | 0.0000 |  |
| GSE113486 | rank | 0.6023 | 0.5000 | 0.0000 | 1.0000 |  |
| GSE113740 | raw | 0.8196 | 0.5000 | 1.0000 | 0.0000 |  |
| GSE113740 | rank | 0.5693 | 0.4952 | 0.0000 | 0.9903 |  |
| GSE138092 | raw/rank | - | - | - | - | Not evaluable: GSE138092 does not contain two evaluable label classes |
