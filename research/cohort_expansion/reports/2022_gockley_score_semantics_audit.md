# Gockley PMML Score-Semantics Audit

Scope: determine which validation score representations are encoded or justified by the frozen `signature-best-10.xml`. Validation only (15 cancers, 40 non-cancer controls); no weights, transforms, thresholds, or model artifacts were fitted or changed.

## Serialized Network Mathematics

- PMML version: `3.0`; top-level children: `Header`, `DataDictionary`, `NeuralNetwork`; producer: `STATISTICA Automated Neural Networks (SANN)` version `2.0`.
- Network: `Ovaries_MLP 14-8-2`, classification; 14 ordered inputs, 8 hidden neurons, 2 output neurons; 128 serialized connection weights and 10 biases. All weights and biases were decoded and used in the frozen computation; both layers use tanh. The 14 input-side `NormContinuous` elements each include shift/scale and two `LinearNorm` points mapping the stored assay range to 0..1.
- Hidden neuron IDs are 14-21; output layer neuron IDs are 22 and 23. Finite tanh outputs lie strictly in (-1,1) (closed range convention [-1,1]). The raw Validation outputs empirically range: malignant `-0.281134` to `0.975767`; benign/control `-0.822951` to `0.980808`.
- Neuron 22 (benign/control): activation `tanh`; bias `8.18883239751340e-001`; incoming weights from hidden neuron IDs 14-21: `14:1.56519495987080e-001, 15:-1.37929619949043e+000, 16:-1.17774015394747e+000, 17:-5.66097288207622e-001, 18:-1.50532303145756e+000, 19:1.42606504606583e+000, 20:4.63638472200952e-002, 21:1.75470836154840e+000`.
- Neuron 23 (malignant): activation `tanh`; bias `4.22284470194467e-002`; incoming weights from hidden neuron IDs 14-21: `14:6.26735611827409e-001, 15:1.06637102984711e+000, 16:2.55792716028695e+000, 17:-3.20594268282178e-001, 18:9.21131991697708e-001, 19:-7.97719723574940e-001, 20:-1.80909016686707e-001, 21:-1.97559053469825e+000`.
- `NeuralOutputs` maps neuron 22 to `benign/control` and neuron 23 to `malignant`. This is a categorical mapping (`NormDiscrete`), not a probability output declaration.
- No `TransformationDictionary`, `LocalTransformations`, `Targets`, `Output`, `ModelVerification`, `Extension`, or `Annotation` elements occur outside the network (observed such elements: []). No probability `OutputField` occurs. The only `DerivedField` transforms are the input `NormContinuous` mappings inside `NeuralInputs`.
- The outputs are separately activated tanh neuron scores. No `normalizationMethod` is set on the output layer, and no equation constrains their sum; in Validation, their observed sum varies from 0.130688 to 1.298486, rather than being fixed at 1. Each output is negative on some rows; the quotient `m/(m+b)` falls outside [0,1] for 15/55 Validation rows (range -2.1512 to 6.3852).
- Neither output is intended by the serialized PMML to be a probability: the XML declares class scores and a maximum-output classification rule. No calibrated probability map, bipolar target coding, sigmoid link, simplemax, or softmax is serialized.

## Validation Score Transformations

AUC uses the same frozen Validation labels for every row. 'Same ranking' compares the ordering of all 55 samples to the PMML-mapped malignant output. Transformations marked unsupported are mathematical diagnostics only; they are not substituted for the published probability or cutoff.

| Representation | Formula | AUC | Same ordering as malignant output? | Probability status | PMML evidence classification |
|---|---|---:|---:|---|---|
| malignant_output | `m = tanh(z_malignant), the raw output of the neuron mapped to malignant` | 0.8867 | True | No: categorical class score, not a probability field. | directly supported by PMML semantics |
| negative_benign_control_output | `-b, where b is the raw output of neuron 22 mapped to benign/control` | 0.8667 | False | No: reverse-coded benign score, not a probability field. | directly supported as an oriented class score |
| malignant_minus_benign_margin | `m - b; positive exactly when the malignant output exceeds the benign/control output` | 0.8867 | False | No: signed decision score, not a probability field. | directly supported by PMML classification semantics |
| tanh_malignant_to_unit_interval | `(m + 1) / 2` | 0.8867 | True | No for the released model: unit-interval rescaling is not calibration and the target coding is absent. | plausible but unsupported |
| inverted_tanh_benign_to_unit_interval | `1 - (b + 1) / 2` | 0.8667 | False | No for the released model: inverse rescaling is not calibrated cancer probability. | plausible but unsupported |
| sigmoid_of_malignant_output | `sigmoid(m) = 1 / (1 + exp(-m))` | 0.8867 | True | No for the released model: sigmoid(m) is monotonic pseudo-probability, not a specified calibration. | plausible but unsupported |
| posthoc_softmax_of_two_tanh_outputs | `exp(m) / (exp(m) + exp(b))` | 0.8867 | False | Would be a normalized pair only if softmax were specified; it is not the released model output. | plausible but unsupported |
| malignant_over_sum_of_outputs | `m / (m + b)` | not calculated | not evaluated | No: outputs are signed and have no positive-sum constraint; the quotient is not a probability definition. | mathematically unjustified |

The inverse benign-node result is AUC 0.8667, close to the published 0.87, but PMML does not designate `-output_benign` as a cancer probability. The designated malignant-node AUC is 0.8867; the exact PMML argmax margin AUC is 0.8867. Similar AUC values do not identify which unsupported probability conversion, if any, the paper’s software used.
- `m/(m+b)` was not scored as an AUC representation: the quotient is outside [0,1] for 15/55 Validation rows (range -2.1512 to 6.3852), and PMML declares no simplemax or positive-sum constraint. A post-hoc two-output softmax is included only as an unsupported diagnostic; it is not in this network's PMML.

## The 50% Cutoff and Published Result

The paper states a 50% predicted-probability cutoff, but the released PMML has no probability transformation. The probability scale is therefore not recoverable from the released model files. Candidate monotone maps (tanh affine rescaling, sigmoid, post-hoc softmax) produce valid-looking numbers but are not selected or calibrated by the XML; no one can be called the paper’s predicted probability based on these files.
- The exact PMML class rule is `output_malignant > output_benign/control`, equivalently margin `m-b > 0`; it is not a 50% probability threshold. No candidate 50% boundary was applied to labels.
- The supplement caption calls `signature-best-10.xml` a neural network trained on early-stage cohort samples; its human target and 14 historical-signature inputs match the Figure 3 task. The paper states that >100,000 networks were compared and the best retained by Test-set assignment, but neither the filename nor PMML records a winner rank, selection flag, or Figure 3 AUC. 'best-10' has no defined technical meaning in the released files and does not mean 10 inputs. The file is the supplied human-model candidate; exact identity as the final selected winner is unverified.
- The PMML producer is STATISTICA Automated Neural Networks (SANN) 2.0, but no software-side probability processing is serialized. Processing outside the PMML is possible but not demonstrated by any supplied file; its formula cannot be recovered here.
- On Validation, the designated malignant-node AUC is 0.8867, versus the published 0.87 (95% CI 0.84-0.94); the value is within the reported interval and approximately reproduces ranking. The inverse benign-node AUC 0.8667 also rounds to 0.87, but closeness does not justify choosing it.
- PMML argmax classification on Validation is 10/15 cancers detected (66.7%) and 38/40 controls correctly negative (95.0%), confusion TN/FP/FN/TP=38/2/5/10. This is not the published 73%/91% point.
- The stated Validation denominator has 15 cancers and 40 controls. 73% is compatible with 11/15=73.3%; 91% specificity is not an attainable integer fraction of 40 controls (36/40=90%, 37/40=92.5%). The exact published specificity therefore does not match this validation denominator under a standard confusion-matrix calculation. The paper passage does not provide a different denominator or exported probability map.

## Final Reproducibility Classification

| Level | Classification | Basis |
|---|---|---|
| Architecture reproducibility | Reproduced for the released XML | All 14 input transforms, neuron counts, activations, weights, biases, and class mappings are serialized and decoded; exact identity as the selected final winner is unverified. |
| Ranking reproducibility | Approximately reproduced | Designated malignant output and PMML margin each give Validation AUC 0.8867 versus reported 0.87; the reported interval 0.84-0.94 contains this value. |
| Operating-point reproducibility | Not reproducible from released files | The 50% predicted-probability scale/conversion is absent; PMML argmax gives 66.7%/95.0%, not 73%/91%. |

## Scientific Classification

- Released human PMML model successfully reconstructed: yes; correspondence to the exact final network selected from >100,000 candidates is unverified.
- Validation ranking consistent with paper: approximately; AUC 0.8867 is close to and within the reported CI, but does not exactly equal 0.87.
- Published 50% probability scale recoverable: no.
- Published operating point recoverable: no.
- Missing: the model’s probability conversion/calibration definition, any selected-network identifier, and a validation-set confusion-count/denominator consistent with the reported 91% specificity. The released artifacts do not establish whether software-side probability processing occurred outside PMML.
- Supported dataset labels: same-cohort model reconstruction; historical-signature replication; approximate internal validation ranking reproduction. Not independent external validation, and not current OVAI validation (0/5 exact current features).

No training, refitting, tuning, recalibration, or PMML modification was performed. The raw workbook and ZIP were read in place.
