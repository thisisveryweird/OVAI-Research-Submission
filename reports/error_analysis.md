# GSE94533 OOF Error Analysis

This analysis uses the saved out-of-fold predictions in `reports/tables/gse94533_enhancement_oof_predictions.csv`. No model was refit and no threshold was selected from this analysis.

## Cohort and Limitations

The joined cohort contains 179 unique GSE94533 patients: 96 invasive cancer and 83 non-cancer. Metadata include diagnosis, histology, stage group, age, and CA-125 availability. No batch/platform field is present in the verified metadata; the `e-sample` versus `numeric-sample` field is only a sample-ID provenance proxy, not a validated batch variable.

The saved models are binary cancer-versus-non-cancer classifiers. Borderline tumors are not relabeled as cancer. The three-class section therefore reports one-vs-rest ranking and cancer-vs-benign-vs-borderline discrimination; it does not claim a three-class hard-label classifier or invent benign/borderline predictions.

## Overall Error Pattern

- `median_logistic`: ROC-AUC 0.8016, balanced accuracy 0.7518, false negatives 28, false positives 17.
- `knn_logistic`: ROC-AUC 0.8126, balanced accuracy 0.7329, false negatives 27, false positives 21.
- `calibrated_svm`: ROC-AUC 0.8084, balanced accuracy 0.7513, false negatives 20, false positives 24.
- `random_forest`: ROC-AUC 0.7784, balanced accuracy 0.7381, false negatives 26, false positives 21.
- `xgboost`: ROC-AUC 0.7938, balanced accuracy 0.7312, false negatives 25, false positives 23.
- `ratio_logistic`: ROC-AUC 0.8020, balanced accuracy 0.7438, false negatives 33, false positives 14.
- `stacked_ensemble`: ROC-AUC 0.7999, balanced accuracy 0.7321, false negatives 26, false positives 22.

## False Negatives

### median_logistic
- False negatives: 28; median probability 0.3137; median age 53.5000; missing clinical values in 6.
- Stage distribution: {"Stage I/II": 17, "Stage III/IV": 11}.
- Histology distribution: {"Endo/Clear Cell": 11, "Endometrial Adenocarcinoma": 2, "Granulosa Cell": 1, "Myxoid Leiomyosarcoma": 1, "Serous": 13}.
- Largest descriptive miRNA contrasts versus correctly classified cancers (log1p scale, standardized mean difference): hsa-miR-200c-3p delta=-1.39 (error mean=2.53, reference mean=3.96); hsa-miR-200a-3p delta=-1.17 (error mean=2.38, reference mean=3.67); hsa-miR-29a-3p delta=-1.09 (error mean=5.25, reference mean=6.04); hsa-miR-320c delta=-0.97 (error mean=4.27, reference mean=5.05); hsa-miR-1307-5p delta=-0.94 (error mean=3.39, reference mean=3.81).
- Per-sample details are in `reports/error_analysis_predictions.csv`.
### knn_logistic
- False negatives: 27; median probability 0.3694; median age 56.0000; missing clinical values in 5.
- Stage distribution: {"Stage I/II": 15, "Stage III/IV": 12}.
- Histology distribution: {"Endo/Clear Cell": 9, "Endometrial Adenocarcinoma": 2, "Granulosa Cell": 1, "Myxoid Leiomyosarcoma": 1, "Serous": 14}.
- Largest descriptive miRNA contrasts versus correctly classified cancers (log1p scale, standardized mean difference): hsa-miR-200c-3p delta=-1.35 (error mean=2.54, reference mean=3.94); hsa-miR-200a-3p delta=-1.11 (error mean=2.41, reference mean=3.64); hsa-miR-320c delta=-1.08 (error mean=4.20, reference mean=5.06); hsa-miR-320d delta=-0.98 (error mean=3.14, reference mean=3.93); hsa-miR-1307-5p delta=-0.71 (error mean=3.45, reference mean=3.78).
- Per-sample details are in `reports/error_analysis_predictions.csv`.
### calibrated_svm
- False negatives: 20; median probability 0.3492; median age 52.0000; missing clinical values in 5.
- Stage distribution: {"Stage I/II": 14, "Stage III/IV": 6}.
- Histology distribution: {"Endo/Clear Cell": 9, "Endometrial Adenocarcinoma": 1, "Granulosa Cell": 1, "Serous": 9}.
- Largest descriptive miRNA contrasts versus correctly classified cancers (log1p scale, standardized mean difference): hsa-miR-200c-3p delta=-1.18 (error mean=2.54, reference mean=3.81); hsa-miR-320c delta=-0.99 (error mean=4.20, reference mean=4.98); hsa-miR-29a-3p delta=-0.95 (error mean=5.23, reference mean=5.96); hsa-miR-1307-5p delta=-0.94 (error mean=3.36, reference mean=3.77); hsa-miR-200a-3p delta=-0.91 (error mean=2.47, reference mean=3.51).
- Per-sample details are in `reports/error_analysis_predictions.csv`.
### random_forest
- False negatives: 26; median probability 0.3915; median age 57.0000; missing clinical values in 9.
- Stage distribution: {"Stage I/II": 15, "Stage III/IV": 11}.
- Histology distribution: {"Endo/Clear Cell": 9, "Endometrial Adenocarcinoma": 1, "Granulosa Cell": 1, "Myxoid Leiomyosarcoma": 1, "Serous": 14}.
- Largest descriptive miRNA contrasts versus correctly classified cancers (log1p scale, standardized mean difference): hsa-miR-200c-3p delta=-1.57 (error mean=2.42, reference mean=3.96); hsa-miR-200a-3p delta=-1.20 (error mean=2.33, reference mean=3.65); hsa-miR-203a delta=-1.03 (error mean=3.07, reference mean=3.93); hsa-miR-29a-3p delta=-0.93 (error mean=5.30, reference mean=6.00); hsa-miR-450b-5p delta=-0.85 (error mean=2.67, reference mean=3.31).
- Per-sample details are in `reports/error_analysis_predictions.csv`.
### xgboost
- False negatives: 25; median probability 0.3218; median age 53.0000; missing clinical values in 7.
- Stage distribution: {"Stage I/II": 14, "Stage III/IV": 11}.
- Histology distribution: {"Endo/Clear Cell": 9, "Endometrial Adenocarcinoma": 1, "Granulosa Cell": 1, "Myxoid Leiomyosarcoma": 1, "Serous": 13}.
- Largest descriptive miRNA contrasts versus correctly classified cancers (log1p scale, standardized mean difference): hsa-miR-200c-3p delta=-1.43 (error mean=2.51, reference mean=3.91); hsa-miR-200a-3p delta=-1.06 (error mean=2.44, reference mean=3.60); hsa-miR-320c delta=-0.91 (error mean=4.26, reference mean=5.02); hsa-miR-320d delta=-0.89 (error mean=3.18, reference mean=3.89); hsa-miR-1307-5p delta=-0.76 (error mean=3.42, reference mean=3.78).
- Per-sample details are in `reports/error_analysis_predictions.csv`.
### ratio_logistic
- False negatives: 33; median probability 0.4043; median age 56.0000; missing clinical values in 8.
- Stage distribution: {"Stage I/II": 21, "Stage III/IV": 12}.
- Histology distribution: {"Endo/Clear Cell": 11, "Endometrial Adenocarcinoma": 2, "Granulosa Cell": 1, "Myxoid Leiomyosarcoma": 1, "Serous": 18}.
- Largest descriptive miRNA contrasts versus correctly classified cancers (log1p scale, standardized mean difference): hsa-miR-200c-3p delta=-1.40 (error mean=2.60, reference mean=4.04); hsa-miR-200a-3p delta=-1.16 (error mean=2.46, reference mean=3.73); hsa-miR-320c delta=-0.90 (error mean=4.33, reference mean=5.08); hsa-miR-29a-3p delta=-0.75 (error mean=5.42, reference mean=6.01); hsa-miR-203a delta=-0.74 (error mean=3.26, reference mean=3.92).
- Per-sample details are in `reports/error_analysis_predictions.csv`.
### stacked_ensemble
- False negatives: 26; median probability 0.3842; median age 55.0000; missing clinical values in 7.
- Stage distribution: {"Stage I/II": 17, "Stage III/IV": 9}.
- Histology distribution: {"Endo/Clear Cell": 12, "Endometrial Adenocarcinoma": 1, "Granulosa Cell": 1, "Myxoid Leiomyosarcoma": 1, "Serous": 11}.
- Largest descriptive miRNA contrasts versus correctly classified cancers (log1p scale, standardized mean difference): hsa-miR-200c-3p delta=-1.35 (error mean=2.52, reference mean=3.92); hsa-miR-29a-3p delta=-1.09 (error mean=5.23, reference mean=6.02); hsa-miR-200a-3p delta=-1.08 (error mean=2.42, reference mean=3.62); hsa-miR-1307-5p delta=-0.83 (error mean=3.41, reference mean=3.79); hsa-miR-320c delta=-0.77 (error mean=4.35, reference mean=5.00).
- Per-sample details are in `reports/error_analysis_predictions.csv`.

## False Positives

### median_logistic
- False positives: 17; median probability 0.6576; median age 60.0000; diagnosis distribution: {"Benign": 8, "Borderline": 7, "Control": 2}.
- Largest descriptive miRNA contrasts versus correctly classified non-cancer samples (log1p scale, standardized mean difference): hsa-miR-200c-3p delta=0.92 (error mean=3.20, reference mean=2.40); hsa-miR-320d delta=0.84 (error mean=3.61, reference mean=3.12); hsa-miR-203a delta=0.76 (error mean=3.79, reference mean=3.13); hsa-miR-320c delta=0.61 (error mean=4.65, reference mean=4.30); hsa-miR-29a-3p delta=0.58 (error mean=5.67, reference mean=5.27).
### knn_logistic
- False positives: 21; median probability 0.5717; median age 56.0000; diagnosis distribution: {"Benign": 12, "Borderline": 7, "Control": 2}.
- Largest descriptive miRNA contrasts versus correctly classified non-cancer samples (log1p scale, standardized mean difference): hsa-miR-335-5p delta=0.84 (error mean=3.32, reference mean=2.81); hsa-miR-200c-3p delta=0.81 (error mean=3.07, reference mean=2.39); hsa-miR-203a delta=0.80 (error mean=3.78, reference mean=3.09); hsa-miR-1307-5p delta=0.75 (error mean=3.69, reference mean=3.35); hsa-miR-320d delta=0.72 (error mean=3.54, reference mean=3.12).
### calibrated_svm
- False positives: 24; median probability 0.5920; median age 56.0000; diagnosis distribution: {"Benign": 15, "Borderline": 6, "Control": 3}.
- Largest descriptive miRNA contrasts versus correctly classified non-cancer samples (log1p scale, standardized mean difference): hsa-miR-200c-3p delta=0.82 (error mean=3.05, reference mean=2.37); hsa-miR-92a-3p delta=-0.73 (error mean=8.84, reference mean=9.38); hsa-miR-320d delta=0.55 (error mean=3.45, reference mean=3.13); hsa-miR-203a delta=0.55 (error mean=3.62, reference mean=3.12); hsa-miR-335-5p delta=0.51 (error mean=3.17, reference mean=2.85).
### random_forest
- False positives: 21; median probability 0.6071; median age 56.0000; diagnosis distribution: {"Benign": 13, "Borderline": 6, "Control": 2}.
- Largest descriptive miRNA contrasts versus correctly classified non-cancer samples (log1p scale, standardized mean difference): hsa-miR-335-5p delta=1.05 (error mean=3.40, reference mean=2.79); hsa-miR-320c delta=0.97 (error mean=4.77, reference mean=4.24); hsa-miR-203a delta=0.92 (error mean=3.85, reference mean=3.06); hsa-miR-1307-5p delta=0.91 (error mean=3.74, reference mean=3.33); hsa-miR-320d delta=0.87 (error mean=3.60, reference mean=3.09).
### xgboost
- False positives: 23; median probability 0.6557; median age 55.0000; diagnosis distribution: {"Benign": 13, "Borderline": 7, "Control": 3}.
- Largest descriptive miRNA contrasts versus correctly classified non-cancer samples (log1p scale, standardized mean difference): hsa-miR-200c-3p delta=1.28 (error mean=3.23, reference mean=2.31); hsa-miR-203a delta=0.96 (error mean=3.85, reference mean=3.03); hsa-miR-29a-3p delta=0.81 (error mean=5.71, reference mean=5.21); hsa-miR-200a-3p delta=0.80 (error mean=3.14, reference mean=2.50); hsa-miR-320d delta=0.77 (error mean=3.54, reference mean=3.10).
### ratio_logistic
- False positives: 14; median probability 0.5916; median age 55.5000; diagnosis distribution: {"Benign": 9, "Borderline": 3, "Control": 2}.
- Largest descriptive miRNA contrasts versus correctly classified non-cancer samples (log1p scale, standardized mean difference): hsa-miR-200c-3p delta=1.31 (error mean=3.46, reference mean=2.38); hsa-miR-203a delta=0.92 (error mean=3.93, reference mean=3.13); hsa-miR-1307-5p delta=0.77 (error mean=3.74, reference mean=3.38); hsa-miR-320d delta=0.74 (error mean=3.60, reference mean=3.15); hsa-miR-320c delta=0.71 (error mean=4.72, reference mean=4.31).
### stacked_ensemble
- False positives: 22; median probability 0.5852; median age 56.0000; diagnosis distribution: {"Benign": 15, "Borderline": 6, "Control": 1}.
- Largest descriptive miRNA contrasts versus correctly classified non-cancer samples (log1p scale, standardized mean difference): hsa-miR-320d delta=0.75 (error mean=3.55, reference mean=3.11); hsa-miR-320c delta=0.72 (error mean=4.68, reference mean=4.26); hsa-miR-335-5p delta=0.63 (error mean=3.23, reference mean=2.84); hsa-miR-203a delta=0.61 (error mean=3.66, reference mean=3.12); hsa-miR-1307-5p delta=0.59 (error mean=3.64, reference mean=3.36).

## Borderline Tumors

Borderline tumors are a distinct non-invasive subgroup in the source labels. They are not silently merged into cancer.
- `median_logistic`, invasive cancer vs all non-cancer: n=179, ROC-AUC=0.8016, balanced accuracy=0.7518.
- `median_logistic`, invasive cancer vs benign: n=142, ROC-AUC=0.7937, balanced accuracy=0.7672.
- `median_logistic`, invasive cancer vs healthy controls: n=111, ROC-AUC=0.8438, balanced accuracy=0.7875.
- `median_logistic`, invasive cancer vs borderline: n=118, ROC-AUC=0.7893, balanced accuracy=0.6951.
- `median_logistic`, three-class ranking (no hard-label 3-class predictions): Cancer one-vs-rest AUC 0.7923; Benign one-vs-rest AUC 0.2618; Borderline one-vs-rest AUC 0.3031.
- `knn_logistic`, invasive cancer vs all non-cancer: n=179, ROC-AUC=0.8126, balanced accuracy=0.7329.
- `knn_logistic`, invasive cancer vs benign: n=142, ROC-AUC=0.8111, balanced accuracy=0.7289.
- `knn_logistic`, invasive cancer vs healthy controls: n=111, ROC-AUC=0.8563, balanced accuracy=0.7927.
- `knn_logistic`, invasive cancer vs borderline: n=118, ROC-AUC=0.7860, balanced accuracy=0.7003.
- `knn_logistic`, three-class ranking (no hard-label 3-class predictions): Cancer one-vs-rest AUC 0.8030; Benign one-vs-rest AUC 0.2432; Borderline one-vs-rest AUC 0.3131.
- `calibrated_svm`, invasive cancer vs all non-cancer: n=179, ROC-AUC=0.8084, balanced accuracy=0.7513.
- `calibrated_svm`, invasive cancer vs benign: n=142, ROC-AUC=0.8039, balanced accuracy=0.7328.
- `calibrated_svm`, invasive cancer vs healthy controls: n=111, ROC-AUC=0.8340, balanced accuracy=0.7958.
- `calibrated_svm`, invasive cancer vs borderline: n=118, ROC-AUC=0.8002, balanced accuracy=0.7595.
- `calibrated_svm`, three-class ranking (no hard-label 3-class predictions): Cancer one-vs-rest AUC 0.8027; Benign one-vs-rest AUC 0.2518; Borderline one-vs-rest AUC 0.2987.
- `random_forest`, invasive cancer vs all non-cancer: n=179, ROC-AUC=0.7784, balanced accuracy=0.7381.
- `random_forest`, invasive cancer vs benign: n=142, ROC-AUC=0.7846, balanced accuracy=0.7233.
- `random_forest`, invasive cancer vs healthy controls: n=111, ROC-AUC=0.8313, balanced accuracy=0.7979.
- `random_forest`, invasive cancer vs borderline: n=118, ROC-AUC=0.7292, balanced accuracy=0.7282.
- `random_forest`, three-class ranking (no hard-label 3-class predictions): Cancer one-vs-rest AUC 0.7667; Benign one-vs-rest AUC 0.2594; Borderline one-vs-rest AUC 0.3608.
- `xgboost`, invasive cancer vs all non-cancer: n=179, ROC-AUC=0.7938, balanced accuracy=0.7312.
- `xgboost`, invasive cancer vs benign: n=142, ROC-AUC=0.8019, balanced accuracy=0.7285.
- `xgboost`, invasive cancer vs healthy controls: n=111, ROC-AUC=0.8278, balanced accuracy=0.7698.
- `xgboost`, invasive cancer vs borderline: n=118, ROC-AUC=0.7538, balanced accuracy=0.7107.
- `xgboost`, three-class ranking (no hard-label 3-class predictions): Cancer one-vs-rest AUC 0.7863; Benign one-vs-rest AUC 0.2491; Borderline one-vs-rest AUC 0.3377.
- `ratio_logistic`, invasive cancer vs all non-cancer: n=179, ROC-AUC=0.8020, balanced accuracy=0.7438.
- `ratio_logistic`, invasive cancer vs benign: n=142, ROC-AUC=0.7971, balanced accuracy=0.7303.
- `ratio_logistic`, invasive cancer vs healthy controls: n=111, ROC-AUC=0.8056, balanced accuracy=0.7615.
- `ratio_logistic`, invasive cancer vs borderline: n=118, ROC-AUC=0.8097, balanced accuracy=0.7599.
- `ratio_logistic`, three-class ranking (no hard-label 3-class predictions): Cancer one-vs-rest AUC 0.8012; Benign one-vs-rest AUC 0.2684; Borderline one-vs-rest AUC 0.2730.
- `stacked_ensemble`, invasive cancer vs all non-cancer: n=179, ROC-AUC=0.7999, balanced accuracy=0.7321.
- `stacked_ensemble`, invasive cancer vs benign: n=142, ROC-AUC=0.7885, balanced accuracy=0.7015.
- `stacked_ensemble`, invasive cancer vs healthy controls: n=111, ROC-AUC=0.8694, balanced accuracy=0.8313.
- `stacked_ensemble`, invasive cancer vs borderline: n=118, ROC-AUC=0.7765, balanced accuracy=0.7282.
- `stacked_ensemble`, three-class ranking (no hard-label 3-class predictions): Cancer one-vs-rest AUC 0.7846; Benign one-vs-rest AUC 0.2607; Borderline one-vs-rest AUC 0.3211.

## Early-Stage Disease

Stage I-only counts are reported separately; the current verified metadata encode most early cases as `Stage I/II`, so Stage I-only analysis is limited.
- `median_logistic`, Stage I/II vs non-cancer: n=135, positive=52, ROC-AUC=0.7706, balanced accuracy=0.7341, sensitivity=0.6731.
- `median_logistic`, Stage III/IV vs non-cancer: n=127, positive=44, ROC-AUC=0.8382, balanced accuracy=0.7726, sensitivity=0.7500.
- `median_logistic`, Stage I only vs non-cancer: n=83, positive=None, status=insufficient_class_variation, ROC-AUC=NA.
- `knn_logistic`, Stage I/II vs non-cancer: n=135, positive=52, ROC-AUC=0.7817, balanced accuracy=0.7293, sensitivity=0.7115.
- `knn_logistic`, Stage III/IV vs non-cancer: n=127, positive=44, ROC-AUC=0.8491, balanced accuracy=0.7371, sensitivity=0.7273.
- `knn_logistic`, Stage I only vs non-cancer: n=83, positive=None, status=insufficient_class_variation, ROC-AUC=NA.
- `calibrated_svm`, Stage I/II vs non-cancer: n=135, positive=52, ROC-AUC=0.7623, balanced accuracy=0.7208, sensitivity=0.7308.
- `calibrated_svm`, Stage III/IV vs non-cancer: n=127, positive=44, ROC-AUC=0.8628, balanced accuracy=0.7872, sensitivity=0.8636.
- `calibrated_svm`, Stage I only vs non-cancer: n=83, positive=None, status=insufficient_class_variation, ROC-AUC=NA.
- `random_forest`, Stage I/II vs non-cancer: n=135, positive=52, ROC-AUC=0.7523, balanced accuracy=0.7293, sensitivity=0.7115.
- `random_forest`, Stage III/IV vs non-cancer: n=127, positive=44, ROC-AUC=0.8091, balanced accuracy=0.7485, sensitivity=0.7500.
- `random_forest`, Stage I only vs non-cancer: n=83, positive=None, status=insufficient_class_variation, ROC-AUC=NA.
- `xgboost`, Stage I/II vs non-cancer: n=135, positive=52, ROC-AUC=0.7850, balanced accuracy=0.7268, sensitivity=0.7308.
- `xgboost`, Stage III/IV vs non-cancer: n=127, positive=44, ROC-AUC=0.8042, balanced accuracy=0.7364, sensitivity=0.7500.
- `xgboost`, Stage I only vs non-cancer: n=83, positive=None, status=insufficient_class_variation, ROC-AUC=NA.
- `ratio_logistic`, Stage I/II vs non-cancer: n=135, positive=52, ROC-AUC=0.7669, balanced accuracy=0.7137, sensitivity=0.5962.
- `ratio_logistic`, Stage III/IV vs non-cancer: n=127, positive=44, ROC-AUC=0.8434, balanced accuracy=0.7793, sensitivity=0.7273.
- `ratio_logistic`, Stage I only vs non-cancer: n=83, positive=None, status=insufficient_class_variation, ROC-AUC=NA.
- `stacked_ensemble`, Stage I/II vs non-cancer: n=135, positive=52, ROC-AUC=0.7755, balanced accuracy=0.7040, sensitivity=0.6731.
- `stacked_ensemble`, Stage III/IV vs non-cancer: n=127, positive=44, ROC-AUC=0.8289, balanced accuracy=0.7652, sensitivity=0.7955.
- `stacked_ensemble`, Stage I only vs non-cancer: n=83, positive=None, status=insufficient_class_variation, ROC-AUC=NA.

## Missingness and Platform Proxies

Clinical missingness is reported as none, one of age/CA-125 missing, or both missing. These are descriptive subgroup analyses, not causal explanations.
No validated platform or batch variable is available in the current GSE94533 metadata. Results split by e-sample versus numeric-sample ID are included only as a provenance proxy and should not be interpreted as a batch effect.

## Outputs

- Per-sample, per-model errors and margins: `reports/error_analysis_predictions.csv`.
- Machine-readable subgroup metrics: `reports/error_analysis_subgroups.csv`.
- Input: saved OOF predictions only; no refitting performed.
