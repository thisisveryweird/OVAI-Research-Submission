# CPTAC Ovarian Proteomics Acquisition Audit

Audit date: 2026-09-19

## Decision

CPTAC ovarian mass-spectrometry data were **not downloaded into the diagnostic training set**. The requested use as "blood proteins alongside serum miRNA" is not supported by the available CPTAC ovarian resources: the public CPTAC client exposes ovarian tumor-proteomics sources, not a verified blood-protein cohort linked to GSE94533.

Downloading tumor proteomics and appending it to GSE94533 serum rows would create invalid patient-level multimodal records unless the same patients, specimen types, collection times, and endpoint labels were verified. It would also change the biological task from circulating cancer detection to tissue tumor characterization.

## Verified Source Information

The public `cptac` client package (version 1.5.14 available from PyPI) defines ovarian (`ov`) sources including:

- BCM
- Broad
- MSSM
- UMich
- WashU
- Harmonized

These are CPTAC ovarian tumor-proteomics source families. No local CPTAC matrix, patient map, blood specimen table, or GSE94533 crosswalk is present in this repository.

The repository's existing compatibility matrix records CPTAC as conditional only. It has unknown local patient counts, unknown overlap with TCGA/GSE94533, and no verified compatible serum diagnostic endpoint.

## Why No Download Was Performed

1. The requested tissue/source is blood, but the identified CPTAC ovarian resources are tumor tissue proteomics.
2. No patient-level linkage to GSE94533 has been established.
3. No verified cancer-versus-benign/borderline/control phenotype aligned to GSE94533 is available locally.
4. A downloaded tissue matrix would be useful for a separate tissue-proteomics or biological discovery analysis, not for the current serum diagnostic fusion model.
5. The download may include multiple source releases and large files; selecting one without a manifest, access status, assay definition, and patient-level deduplication would reduce reproducibility.

## Safe Next Step

A future CPTAC analysis should first select one named PDC/CPTAC release, save its project/file manifest and checksums, identify specimen type and assay workflow, normalize to one row per patient, audit repeated samples and overlap with TCGA/other CPTAC sources, and evaluate it as a separate tissue endpoint. It should not be merged with GSE94533 unless a validated patient/specimen mapping and compatible clinical endpoint exist.

## Conclusion

CPTAC does not currently provide a scientifically valid way to add "blood proteins" to the GSE94533 serum-miRNA model. The requested accuracy improvement cannot be claimed from this source without changing the endpoint or fabricating multimodal patient rows.
