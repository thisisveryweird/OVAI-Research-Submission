# Published 2024 Model and Cohort Reproduction Assessment

## Reproduced from the Workbook

- Workbook file: `research/cohort_expansion/data/raw/2024_adnexal_mass/NIHMS2021079-supplement-1.xlsx`; 1,744,787 bytes; SHA256 `e98f5a60fe3117fc69cd9243bc60f146979b831b20dbb6dee4d8e7ce6affc748`.
- Cohort reconstruction: 678 unique IDs, 286 cancer labels, 392 control labels.
- development/training: n=468, cancer=191, control=277.
- external validation: n=110, cancer=51, control=59.
- internal validation: n=100, cancer=44, control=56.
- External validation subset: n=110, 51 cancer and 59 controls; workbook has 20 unambiguous Stage I/II cancer IDs and 31 Stage III/IV IDs.
- The workbook reproduces the published split and label counts exactly. The external stage group counts (20 Stage I/II, 31 Stage III/IV) also agree.
- The workbook demographics reproduce 290 premenopausal subjects; mean age and postmenopausal counts by source group are summarized in the workbook audit.

## Published Combined-Model miRNA Crosswalk

Exact matches to the 10 miRNAs listed in the paper's model-feature table: 9/10.

| Published feature | Exact match | Measured |
|---|---:|---:|
| `hsa-miR-320e` | `yes` | `yes` |
| `hsa-miR-195-5p` | `yes` | `yes` |
| `hsa-miR-103-3p` | `no` | `no` |
| `hsa-miR-320c` | `yes` | `yes` |
| `cel-miR-39-3p` | `yes` | `yes` |
| `hsa-miR-200c-3p` | `yes` | `yes` |
| `hsa-miR-770-5p` | `yes` | `yes` |
| `hsa-miR-122-5p` | `yes` | `yes` |
| `hsa-miR-19a-3p` | `yes` | `yes` |
| `hsa-miR-424-5p` | `yes` | `yes` |

The paper lists `hsa-miR-103-3p`, while the workbook lists `hsa-mir-103a-3p` (and a distinct `hsa-mir-103a-2-5p`). These are not treated as identical without assay/nomenclature documentation. The other listed features, including `cel-miR-39-3p`, are exact workbook identifiers.

## Performance Reproduction Decision

The paper describes a neural-network combined model, forward regression with ten-fold cross-validation to select the number of miRNAs, and a training-derived probability threshold chosen for at least 80% cross-validated sensitivity. It reports combined-model threshold 0.45, external AUC 0.95 (95% CI 0.90-0.98), sensitivity 0.92, and specificity 0.80. The workbook contains the source patient-level matrix but does not contain fitted network weights or a fully specified executable network architecture, selected-feature coefficients/weights, training fold assignments, or all preprocessing details. Further, one published selected-miRNA identifier does not have an exact workbook-column match.

Accordingly, the cohort/split structure and published feature availability were audited, but published validation predictions and performance were **not recreated**. Fitting a new network or reselecting features on this workbook would be a new model-development exercise, not reproduction of the published frozen model. No such fit was performed in this run.
