# Gockley Published-Experiment Reproduction Audit

Scope: compare the actual workbook assignments and frozen PMML outputs with the article’s human validation experiment. No fitting, refitting, normalization, threshold selection, or model modification was performed.

## Endpoint Definitions

- Published-style: all 75 workbook cancer cases (FIGO I-IIIA2) versus 100 benign adnexal masses plus 100 healthy controls (n=275). The workbook has 2 Stage IIIA1 and 1 Stage IIIA2 case; all 75 belong in this endpoint.
- OVAI-stricter: only Stage I/II cancer versus benign plus healthy controls (72 vs 200; n=272). The three Stage IIIA cases are excluded, not relabeled.

## Workbook Split Check

The workbook has exactly three labels, no duplicate patient IDs, and every row belongs to one set. Counts are 165 Training, 55 Test, and 55 Validation: an exact 3:1:1 allocation of 275 rows. The diagnosis and stage counts reconcile within each split. This verifies aggregate allocation against the paper; individual sample-to-split assignments cannot be cross-checked against a separate published roster because none is supplied in the article text inspected here.
| Split | n | Cancer | Benign | Healthy | Stage I | Stage II | IIIA1 | IIIA2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| All | 275 | 75 | 100 | 100 | 50 | 22 | 2 | 1 |
| Training | 165 | 45 | 60 | 60 | 30 | 13 | 1 | 1 |
| Test | 55 | 15 | 20 | 20 | 9 | 5 | 1 | 0 |
| Validation | 55 | 15 | 20 | 20 | 11 | 4 | 0 | 0 |

## Which XML Is the Human Model?

| Archive member | PMML model name | Inputs (architecture) | Activations / output normalization; input transforms | Target and mapped output neurons | Interpretation |
|---|---|---:|---|---|---|
| `7miRTrimmedXenoNet-38.xml` | `Normalized 2_MLP 7-5-2` | 7 (7-5-2) | tanh, identity; output normalization `softmax`; 7 input transforms, complete=True | `Group`; 12=Cancer, 13=Control | Murine xenograft model: seven features, Group Cancer/Control target, and explicit softmax output normalization. |
| `NewNetXenograft-1.xml` | `Normalized 2_MLP 14-6-2` | 14 (14-6-2) | logistic, exponential; output normalization `None`; 14 input transforms, complete=True | `Group`; 20=Cancer, 21=Control | Murine xenograft model: 14 features, Group Cancer/Control target, and logistic/exponential layers. |
| `signature-best-10.xml` | `Ovaries_MLP 14-8-2` | 14 (14-8-2) | tanh, tanh; output normalization `None`; 14 input transforms, complete=True | `Malignant vs benign/control`; 22=benign/control, 23=malignant | Human early-stage model: human malignancy-versus-benign/control target and 14 historical signature inputs. |

- `7miRTrimmedXenoNet-38.xml` input order: `hsa-mir-150-5p`, `hsa-mir-200a-3p`, `hsa-mir-200c-3p`, `hsa-mir-203a-3p`, `hsa-mir-335-5p`, `hsa-mir-450b-5p`, `hsa-mir-320d`. Each input has a `NormContinuous` transform with two `LinearNorm` points: 7 transforms; all complete=True. PMML producer: STATISTICA Automated Neural Networks (SANN) 2.0; annotations=0; XML comments present=False. Supplement caption: empirically optimized neural network on sera from three ovarian cancer xenografts using seven miRNAs.
- `NewNetXenograft-1.xml` input order: `hsa-mir-1246`, `hsa-mir-1307-5p`, `hsa-mir-150-5p`, `hsa-mir-200a-3p`, `hsa-mir-200c-3p`, `hsa-mir-203a-3p`, `hsa-mir-23b-3p`, `hsa-mir-29a-3p`, `hsa-mir-320c`, `hsa-mir-320d`, `hsa-mir-32-5p`, `hsa-mir-335-5p`, `hsa-mir-450b-5p`, `hsa-mir-92a-3p`. Each input has a `NormContinuous` transform with two `LinearNorm` points: 14 transforms; all complete=True. PMML producer: STATISTICA Automated Neural Networks (SANN) 2.0; annotations=0; XML comments present=False. Supplement caption: neural network trained on sera from three ovarian cancer xenografts using 14 miRNAs.
- `signature-best-10.xml` input order: `hsa-miR-200c-3p`, `hsa-miR-1246`, `hsa-miR-200a-3p`, `hsa-miR-29a-3p`, `hsa-miR-320d`, `hsa-miR-150-5p`, `hsa-miR-203a`, `hsa-miR-23b-3p`, `hsa-miR-92a-3p`, `hsa-miR-32-5p`, `hsa-miR-335-5p`, `hsa-miR-320c`, `hsa-miR-450b-5p`, `hsa-miR-1307-5p`. Each input has a `NormContinuous` transform with two `LinearNorm` points: 14 transforms; all complete=True. PMML producer: STATISTICA Automated Neural Networks (SANN) 2.0; annotations=0; XML comments present=False. Supplement caption: neural network trained on early-stage cohort samples.

The archive caption identifies `signature-best-10.xml` as trained on early-stage cohort samples. Its target is `Malignant vs benign/control`, and its 14 inputs exactly match the historical human 14-miRNA signature. The other two captions explicitly identify xenograft models; their target is `Group` with Cancer/Control categories. `7miRTrimmedXenoNet-38.xml` has seven inputs and softmax output normalization; `NewNetXenograft-1.xml` has 14 inputs but an animal-model target and logistic/exponential layers. Therefore `signature-best-10.xml` is the supplied human Figure 3 candidate; its filename says best-10 but its model declares 14 inputs.

## Probability and Threshold Semantics

`signature-best-10.xml` has a tanh output layer with two mapped class neurons (22=benign/control, 23=malignant). It has no output `normalizationMethod` (`softmax`/`simplemax`), no probability output field, no calibration transform, and no numeric threshold. Its outputs are bounded neural-network scores, not probabilities. The exact cancer score is the raw activation of the mapped malignant neuron; the binary PMML class rule selects the larger output, equivalently malignant-minus-benign/control margin > 0. That argmax boundary is not the paper’s 50% predicted-probability cutoff.
- The paper’s 50% probability threshold is not recoverable from this PMML alone. Missing information is the mathematical mapping from the tanh class scores to calibrated cancer probability (or an explicitly exported probability output/normalization rule). The class mapping and weights do not specify this conversion. No guessed tanh-to-probability transform was applied.

## Frozen PMML Results

ROC-AUC is calculated first from the explicit malignant output neuron; the output margin is also shown because the two class neurons are not probability-normalized. Sensitivity/specificity and confusion matrices use only the PMML maximum-output-neuron rule; they are not presented as results at the paper’s 50% probability cutoff.
| Endpoint | Split | Cases | Controls | n | AUC cancer output | AUC output margin | Sensitivity (argmax) | Specificity (argmax) | TN/FP/FN/TP |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| published_stage_I_to_IIIA2 | All | 75 | 200 | 275 | 0.8394 | 0.8385 | 0.5200 | 0.9600 | 192/8/36/39 |
| published_stage_I_to_IIIA2 | Training | 45 | 120 | 165 | 0.8607 | 0.8637 | 0.4889 | 0.9500 | 114/6/23/22 |
| published_stage_I_to_IIIA2 | Test | 15 | 40 | 55 | 0.7483 | 0.7317 | 0.4667 | 1.0000 | 40/0/8/7 |
| published_stage_I_to_IIIA2 | Validation | 15 | 40 | 55 | 0.8867 | 0.8867 | 0.6667 | 0.9500 | 38/2/5/10 |
| strict_stage_I_II | All | 72 | 200 | 272 | 0.8427 | 0.8426 | 0.5278 | 0.9600 | 192/8/34/38 |
| strict_stage_I_II | Training | 43 | 120 | 163 | 0.8686 | 0.8731 | 0.4884 | 0.9500 | 114/6/22/21 |
| strict_stage_I_II | Test | 14 | 40 | 54 | 0.7304 | 0.7196 | 0.5000 | 1.0000 | 40/0/7/7 |
| strict_stage_I_II | Validation | 15 | 40 | 55 | 0.8867 | 0.8867 | 0.6667 | 0.9500 | 38/2/5/10 |

Validation, published-style endpoint: n=55 (15 cancer, 40 controls); AUC=0.8867 from the malignant output neuron (margin AUC=0.8867). PMML argmax gives sensitivity=0.6667, specificity=0.9500, confusion TN/FP/FN/TP=38/2/5/10. The stricter endpoint has the same validation n and scores because this split contains no Stage IIIA cases (AUC=0.8867).
- Published Figure 3: AUC 0.87, sensitivity 73%, specificity 91% at predicted probability 50%. Reconstructed AUC 0.8867 is close but not identical (difference +0.0167); the PMML argmax operating point is 66.7%/95.0%, not 73%/91%. The published cutoff cannot be tested without the missing probability conversion.
- Integer-count check: validation has 15 cancer and 40 non-cancer rows. 73% sensitivity is consistent with 11/15 (73.3%); 91% specificity cannot be obtained from 40 controls as an integer confusion count (36/40=90.0%; 37/40=92.5%). Thus the published operating percentages do not exactly correspond to an integer confusion matrix on the workbook Validation rows as counted. The article passage does not separately state another denominator for those percentages, so the discrepancy remains unresolved from the files.
- Diagnostic only, not an accepted substitute: using the inverse benign-node score gives validation AUC 0.8667, which rounds to 0.87, but the XML maps the malignant class to neuron 23 and does not state that the inverse benign score was used. This coincidence is insufficient to claim reproduction.

## Model-Selection Independence

The paper says human samples were randomized 3:1:1 and that more than 100,000 networks were built, retaining the best based on Test-set class assignment. Under that stated procedure, Training is used for fitting; Test is used for model selection and is therefore not an unbiased final evaluation after selecting among 100,000 candidates; Validation is the only set described as reserved for final evaluation. It is an internal held-out split from the same 275-subject cohort, not an external cohort. This independence interpretation assumes Validation was not consulted for model selection or threshold calibration; the paper reports it as the validation set, while the PMML itself has no row-level training lineage.

No model was fit, refit, tuned, recalibrated, or altered. Original workbook and ZIP were read in place. Detailed split counts are in `2022_gockley_published_split_stage_counts.csv`; fixed-score metrics are in `2022_gockley_published_reproduction_metrics.csv`.
