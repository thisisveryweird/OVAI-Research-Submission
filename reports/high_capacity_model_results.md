# High-Capacity Local Model Comparison

All preprocessing, interaction construction, and SVM hyperparameter tuning occur inside each outer training fold. Endpoint: GSE94533 Stage I/II cancer versus benign/control.

| Model | Mean OOF ROC-AUC | Balanced Accuracy | Sensitivity | Specificity | >0.90 ROC-AUC |
|---|---:|---:|---:|---:|---|
| frozen_logistic | 0.8994 | 0.8657 | 0.8462 | 0.8852 | No |
| rbf_svm | 0.8884 | 0.8657 | 0.8462 | 0.8852 | No |
| shallow_random_forest | 0.8354 | 0.7945 | 0.7692 | 0.8197 | No |
| engineered_logistic | 0.8997 | 0.8643 | 0.8269 | 0.9016 | No |
