# Partial Historical-Signature Replication

**partial historical-signature replication** is a descriptive, feature-by-feature comparison of the eight exact shared historical miRNAs. It is not validation of the original 14-feature model. Six of the original 14 inputs are absent from the 2024 FirePlex panel. No composite score, coefficient, or threshold was constructed or tuned.

GSE94533 source direction uses the saved, patient-deduplicated Stage I/II cancer versus benign/control endpoint (52 cancer, 61 benign/control). The 2024 direction uses all 286 malignant versus 392 pooled controls. Effect sizes are Hedges' g on each platform's native values, so magnitudes should not be compared across platforms. Direction agreement means only that the feature-wise case-minus-control mean contrasts have the same sign.

The saved stability report provides selection frequencies, not fitted coefficient signs. GSE abundance direction is reported as a separate source-cohort contrast, not misrepresented as a model coefficient.

| Historical feature | GSE94533 direction (g) | 2024 direction (g) | Direction agrees | Prior stable signature member |
|---|---|---|---:|---:|
| `hsa-miR-23b-3p` | higher (g=0.058) | lower (g=-0.113) | no | no |
| `hsa-miR-29a-3p` | higher (g=0.366) | lower (g=-0.200) | no | no |
| `hsa-miR-32-5p` | higher (g=0.216) | higher (g=0.195) | yes | no |
| `hsa-miR-92a-3p` | lower (g=-0.213) | lower (g=-0.463) | yes | no |
| `hsa-miR-150-5p` | lower (g=-0.035) | lower (g=-0.327) | yes | no |
| `hsa-miR-200c-3p` | higher (g=0.580) | higher (g=0.624) | yes | yes |
| `hsa-miR-320c` | higher (g=0.539) | higher (g=0.192) | yes | no |
| `hsa-miR-320d` | higher (g=0.526) | higher (g=0.096) | yes | no |

Directions agree for 6/8 exact overlaps. The only feature meeting the predeclared stability rule (selection in at least 3/5 folds within an existing 5/10/14-feature signature) is `hsa-miR-200c-3p`; it was stable for the all-cancer endpoint, not the Stage I/II endpoint. The strongest cross-platform directional agreement involving a previously stable feature is miR-200c-3p (GSE g=0.580; FirePlex g=0.624).

Interpretation: miR-200c-3p is the clearest repeatable component in this feature-level comparison. miR-29a-3p reverses direction and miR-23b-3p is near zero/opposite. These descriptive univariate associations do not demonstrate a transportable multivariate assay or clinical performance.

Source artifacts: `feature_selection_stability.csv`, `error_analysis_predictions.csv`, and the audited 2024 workbook. No 2024 labels were used in feature selection or model fitting.
