# Gockley 2022 Workbook and Frozen-Signature Audit

## Source Files and Provenance

- `research/cohort_expansion/data/raw/2024_adnexal_mass/DataSheet_2.xlsx`: Excel OOXML workbook; 488924 bytes; SHA256 `29529a5e2b2f17f11f0748df3c0a9270e82e445053bf8bbb13ba8b5fcc11b74a`; modified `2026-09-26 17:15:06.680407+00:00`.
- `research/cohort_expansion/data/raw/2024_adnexal_mass/DataSheet_3.zip`: ZIP archive of PMML XML networks; 8695 bytes; SHA256 `e61793d923ced1d245498a074b5e441a587a3a35aae1d08de22daf4799218e33`; modified `2026-09-26 17:15:06.638407+00:00`.
- DataSheet_3 ZIP integrity: PASS.
- ZIP members (read in place; no extraction): `7miRTrimmedXenoNet-38.xml`, `NewNetXenograft-1.xml`, `signature-best-10.xml`.
- Both files are physically stored under `research/cohort_expansion/data/raw/2024_adnexal_mass/`; the workbook basename present in the repository is `DataSheet_2.xlsx` (underscore), not `DataSheet 2.xlsx`. They were not moved or renamed.
- Workbook sample labels, exact 14-feature panel, and the paired PMML inputs identify this content as the Gockley human serum cohort despite its directory placement.

## Workbook Structure and Cohort

- Worksheet(s): `Sheet1`; used range `A1:GW276`; 275 patient/sample rows and 205 columns including the header.
- Unique patient IDs: 275; duplicate patient IDs: 0; exact duplicate data rows: 0. Expression is one row per patient/sample ID.
- Diagnosis labels: {"Benign": 100, "Cancer Early": 75, "Control": 100}.
- Cancer stage strings: {"1A": 25, "1B": 3, "1C": 22, "2A": 14, "2B": 8, "3A1": 2, "3A2": 1}; grouped stages: {"Stage I/II": 72, "Stage IIIA": 3}.
- The workbook independently contains 75 cancer records, 100 benign masses, and 100 healthy controls; 72 are Stage I/II and 3 are Stage IIIA. Thus all 275 paper-described human serum samples are represented in the workbook.
- The diagnosis field labels all 75 cancer records `Cancer Early`, including the three Stage IIIA cases; use the separate stage field for the Stage I/II endpoint.
- Data split labels: {"Test": 55, "Training": 165, "Validation": 55}.
- Clinical variables: age is available for 275/275; CA-125 is not a workbook column; histopathology is present for 175/275. Benign mass and healthy control are explicitly distinguishable by diagnosis.
- Histopathology counts for non-missing rows: {"carcinosarcoma": 2, "clear cell": 23, "endometrioid": 29, "endometrioma": 25, "mixed serous and transitional cell": 1, "mucinous": 1, "serous": 18, "serous cystadenoma": 75, "transitional cell": 1}. The final ` Histology` field is populated for 75 rows; `Mutation status`, `Stage`, and `Grade` are repeated in trailing columns and exactly duplicate the leading fields.
- Assay block: 193 unique identifiers in columns 9-201; identifiers use mature-name strings such as `hsa-miR-200c-3p`. All 53075 patient-feature values are present (missing assay values: 0).
- Blank worksheet cells total 1500, mostly clinical fields repeated for non-cancer samples; age and all assay values are complete. The canonical first-occurrence clinical fields have missingness for histopathology=100, mutation status=200, stage=200, and grade=200.
- Normalization: the workbook contains no explicit normalization-method column. Its supplementary-data description states that serum miRNA expression values are log10-normalized; treat this as source-caption metadata, not a transform inferred solely from cell values. There is no per-sample CA-125 measurement to join from this workbook.

## Exact Feature Compatibility

- Current OVAI 5-miRNA signature: 0/5 exact matches (none). No compatible frozen current-OVAI test is possible; no substitute features were used.
- Historical 14-miRNA signature: 14/14 exact matches (`hsa-miR-23b-3p`, `hsa-miR-29a-3p`, `hsa-miR-32-5p`, `hsa-miR-92a-3p`, `hsa-miR-150-5p`, `hsa-miR-200a-3p`, `hsa-miR-200c-3p`, `hsa-miR-203a`, `hsa-miR-320c`, `hsa-miR-320d`, `hsa-miR-335-5p`, `hsa-miR-450b-5p`, `hsa-miR-1246`, `hsa-miR-1307-5p`). These are exact workbook identifier matches, not family-level or arm-level substitutions.
- Alias audit: no valid same-mature-product alias was needed or accepted. `hsa-miR-32-5p` is not `hsa-miR-32-3p`; measured miR-200a/200c members are not miR-200b. Full inventory and requested-feature rows are in `2022_gockley_feature_map.csv`.

## Data Sheet 3 PMML Inspection

`7miRTrimmedXenoNet-38.xml`: `Normalized 2_MLP 7-5-2`, architecture 7-5-2; activations tanh, identity; 45 connection weights and 7 biases; 7 explicit input transforms.
`NewNetXenograft-1.xml`: `Normalized 2_MLP 14-6-2`, architecture 14-6-2; activations logistic, exponential; 96 connection weights and 8 biases; 14 explicit input transforms.
`signature-best-10.xml`: `Ovaries_MLP 14-8-2`, architecture 14-8-2; activations tanh, tanh; 128 connection weights and 10 biases; 14 explicit input transforms.
- `signature-best-10.xml` declares all required model components: model name `Ovaries_MLP 14-8-2`, 14 ordered inputs (`hsa-miR-200c-3p`, `hsa-miR-1246`, `hsa-miR-200a-3p`, `hsa-miR-29a-3p`, `hsa-miR-320d`, `hsa-miR-150-5p`, `hsa-miR-203a`, `hsa-miR-23b-3p`, `hsa-miR-92a-3p`, `hsa-miR-32-5p`, `hsa-miR-335-5p`, `hsa-miR-320c`, `hsa-miR-450b-5p`, `hsa-miR-1307-5p`), per-input `NormContinuous` transforms, 2 neural layers, 128 connection weights, 10 biases, activation functions `tanh, tanh`, and explicit output mapping `{"22": "benign/control", "23": "malignant"}`.
- Every PMML shift/scale pair agrees with its two `LinearNorm` points within tolerance: True. The model has no explicit numeric threshold and no `ModelVerification` examples. PMML class decisions are reproduced by the largest output neuron; ROC-AUC uses malignant-minus-benign/control output margin. No probability calibration or threshold tuning was performed.
- The filename says `best-10`, but this PMML network actually declares 14 inputs and an `Ovaries_MLP 14-8-2` architecture. It can be reconstructed from the XML; it is not a 10-feature network.
- The workbook records Training/Validation/Test assignments. Results for all records are apparent cohort performance; Test/Validation rows are reported separately by their workbook `Set`, but the PMML does not encode a training lineage proving which rows were excluded during fitting.

## Frozen Published-Model Results

Metrics use the frozen PMML network only. ROC-AUC ranks the output margin. Sensitivity/specificity and confusion counts use the PMML largest-output-neuron decision, not a threshold selected on Gockley labels. The requested primary endpoint is shown for the full cohort and the workbook Test split below.
| Endpoint | Subset | Cases | Controls | ROC-AUC (margin) | Sensitivity (PMML argmax) | Specificity (PMML argmax) | TN/FP/FN/TP |
|---|---|---:|---:|---:|---:|---:|---|
| stage_I_II_cancer_vs_benign_mass | All | 72 | 100 | 0.8214 | 0.5278 | 0.9300 | 93/7/34/38 |
| stage_I_II_cancer_vs_benign_mass | Test | 14 | 20 | 0.7143 | 0.5000 | 1.0000 | 20/0/7/7 |

### All Endpoints by Workbook Split

| Endpoint | Subset | Cases | Controls | ROC-AUC | Sensitivity | Specificity | TN/FP/FN/TP |
|---|---|---:|---:|---:|---:|---:|---|
| stage_I_II_cancer_vs_benign_mass | All | 72 | 100 | 0.8214 | 0.5278 | 0.9300 | 93/7/34/38 |
| stage_I_II_cancer_vs_benign_mass | Test | 14 | 20 | 0.7143 | 0.5000 | 1.0000 | 20/0/7/7 |
| stage_I_II_cancer_vs_benign_mass | Training | 43 | 60 | 0.8481 | 0.4884 | 0.9167 | 55/5/22/21 |
| stage_I_II_cancer_vs_benign_mass | Validation | 15 | 20 | 0.8500 | 0.6667 | 0.9000 | 18/2/5/10 |
| stage_I_II_cancer_vs_benign_plus_healthy | All | 72 | 200 | 0.8426 | 0.5278 | 0.9600 | 192/8/34/38 |
| stage_I_II_cancer_vs_benign_plus_healthy | Test | 14 | 40 | 0.7196 | 0.5000 | 1.0000 | 40/0/7/7 |
| stage_I_II_cancer_vs_benign_plus_healthy | Training | 43 | 120 | 0.8731 | 0.4884 | 0.9500 | 114/6/22/21 |
| stage_I_II_cancer_vs_benign_plus_healthy | Validation | 15 | 40 | 0.8867 | 0.6667 | 0.9500 | 38/2/5/10 |
| all_cancer_vs_benign_mass | All | 75 | 100 | 0.8165 | 0.5200 | 0.9300 | 93/7/36/39 |
| all_cancer_vs_benign_mass | Test | 15 | 20 | 0.7267 | 0.4667 | 1.0000 | 20/0/8/7 |
| all_cancer_vs_benign_mass | Training | 45 | 60 | 0.8370 | 0.4889 | 0.9167 | 55/5/23/22 |
| all_cancer_vs_benign_mass | Validation | 15 | 20 | 0.8500 | 0.6667 | 0.9000 | 18/2/5/10 |
| all_cancer_vs_benign_plus_healthy | All | 75 | 200 | 0.8385 | 0.5200 | 0.9600 | 192/8/36/39 |
| all_cancer_vs_benign_plus_healthy | Test | 15 | 40 | 0.7317 | 0.4667 | 1.0000 | 40/0/8/7 |
| all_cancer_vs_benign_plus_healthy | Training | 45 | 120 | 0.8637 | 0.4889 | 0.9500 | 114/6/23/22 |
| all_cancer_vs_benign_plus_healthy | Validation | 15 | 40 | 0.8867 | 0.6667 | 0.9500 | 38/2/5/10 |

All four requested endpoints and all workbook splits are provided in `2022_gockley_signature_best_10_metrics.csv`. Patient-level fixed-model outputs are in `2022_gockley_signature_best_10_predictions.csv`.

## GSE94533 Compatibility

- The Gockley source identifies serum small-RNA sequencing. Local GSE94533 artifacts provide serum miRNA TPM values and 2,582 features, consistent with sequencing-derived expression, but do not identify the exact library chemistry. Gockley’s matrix is described as log10-normalized; local OVAI documentation applies `log1p` to GSE94533 TPM values. These are not directly interchangeable raw scales.
- Exact historical feature availability: Gockley 14/14; GSE94533 prediction table 14/14. Current OVAI production features: Gockley 0/5; the production manifest identifies the five current model inputs for GSE94533, but those five are absent from this Gockley workbook.
- Cohort/endpoint: Gockley has 72 Stage I/II cases, 100 benign masses, and 100 healthy controls, with age but no CA-125. GSE94533 has a locally documented Stage I/II-versus-benign/control endpoint of 52 cases and 61 controls in its saved analysis; the broader joined dataset is 179 patients (96 cancers, 83 non-cancer) and includes age/CA-125 where available.
- Identifier comparison: exact normalized ID matches between the 275 Gockley IDs and 179 unique local GSE94533 sample labels: 0. No deterministic patient linkage was found; different ID namespaces mean zero string matches do not alone prove patient independence.
- Sequencing library chemistry and miRBase release are not documented in the workbook or the inspected local GSE94533 artifacts, so equivalence on those details is unverified.
- Decision: suitable for frozen replay of the historical network and biological replication, but not independent external validation of that network because the PMML was trained on this same cohort. The workbook Test/Validation splits are reported separately, although the XML does not verify split-specific training lineage. It is suitable for a separate cross-cohort benchmark after a prespecified scale/feature protocol; it is not suitable for direct frozen testing of current OVAI (0/5 exact features). Do not pool at this stage. Eventual pooled development would require harmonized quantification, verified patient independence, and a new protocol approved before any training.

## Use Decision

Gockley should replace the 2024 FirePlex cohort as the primary clinically aligned cohort-expansion candidate relative to GSE94533: it has 72 Stage I/II cases, separately identifiable benign masses and healthy controls, and all 14 historical inputs. This makes it stronger for historical-signature benchmarking, but the same-cohort PMML results are not independent external validation. It does not replace the need for a valid independent external cohort for current 5-miRNA OVAI, which has 0/5 exact features here. Keep FirePlex as a separate cross-platform biological replication dataset; do not combine either dataset with GSE94533 in this audit.

No model was trained, recalibrated, or tuned. `production_artifacts/` was read only. Original workbook and ZIP remain unchanged.
