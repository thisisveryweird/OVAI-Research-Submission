# Cohort expansion research plan

This directory is intentionally separate from the canonical OVAI artifacts in [production_artifacts](../../production_artifacts).

## Scientific rule

The canonical source model and its results remain frozen. This directory is only for cohort-expansion analysis and assay-compatibility auditing.

## Current status of the 2024 serum adnexal-mass study

The study is a promising high-priority replication candidate because it contains a larger adnexal-mass dataset with both serum miRNA and clinical/protein features. However, the raw individual-level supplementary Data File S1 is not present in this workspace and has not been fully audited. The published study-level counts alone are not sufficient for patient-level OVAI modeling.

## Allowed activities here

- cohort acquisition tracking
- accession-level compatibility checks
- assay/platform documentation
- feature mapping against OVAI miRNAs
- external validation planning
- leakage-safe evaluation design

## Not allowed here

- overwriting canonical artifacts
- training on external labels for feature selection or threshold tuning
- creating patient-level data by inference
- pooling incompatible assay systems without explicit harmonization
- treating literature-level counts as patient-level evidence

## Main audit outputs

- [reports/2024_adnexal_mass_data_audit.md](../../reports/2024_adnexal_mass_data_audit.md)
- [reports/2024_adnexal_mass_feature_map.csv](../../reports/2024_adnexal_mass_feature_map.csv)
- [reports/2024_adnexal_mass_compatibility.csv](../../reports/2024_adnexal_mass_compatibility.csv)

## Execution

Run the audit script with:

python3 research/cohort_expansion/audit_2024_adnexal_mass.py

The script writes the audit documents and leaves the canonical production artifacts untouched.
