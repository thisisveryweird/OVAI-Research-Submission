# Data sources draft — not an acquisition/permission registry

Counts and endpoint evidence status are in `ENDPOINTS.csv`; unknown means unavailable or not independently verified in this bounded offline audit. Public access does not itself establish unrestricted reuse, consent or clinical eligibility.

| Source | Primary specimen / platform | Citation evidence | Reuse / IRB | Present evidence |
|---|---|---|---|---|
| GSE94533 | serum; sequencing GPL18573, reported NextSeq500, miRBase20 | Gockley 2022; saved research reports | unknown | source processed inputs absent; source endpoint characterization only |
| GSE106817 | serum / GPL21263 array per saved evaluator | citation not independently verified | unknown | saved phase5 predictions; raw matrix/platform absent |
| GSE113486 | serum / GPL21263 array | local SOFT PMID30382619, bladder detection title | unknown | local SOFT metadata reused and stdlib eligible-label replay; expressions not evaluated |
| GSE113740 | serum array per existing scouting evidence; primary source unverified | unknown | unknown | bulk report only; 1817 PDF total vs1058 evaluated |
| GSE211692 | serum 300uL / GPL21263 Toray 3D-Gene | PMID36426871; DOI10.1093/jncics/pkac080; PMC9825310 in saved research report | unknown | local SOFT metadata reused; raw matrix not freshly verified |
| 2024 adnexal mass | serum / FirePlex 180-miRNA panel plus proteins | DOI10.1016/j.ygyno.2024.08.008 in saved reports | unknown | extracted workbook cells and ZIP inventory reused; no rendered inspection |
| Gockley supplement | serum / historical assay; original PMML14-8-2 | saved Gockley reports; exact publication identity not freshly checked | unknown | PMML score semantics and cohort reconstruction reports; no model run |

Mature-arm/MIMAT mappings must be reconciled against primary assay annotations, not punctuation normalization alone. Five-feature assay availability does not imply availability of nine-input frozen portal validation, particularly without CA125, age, clinical timing and donor linkage. Signed normalized arrays cannot simply be treated as source TPM.

Roadmap: obtain immutable source rosters, checksums, raw/processed assay definitions, missingness and clinical units; secure citation, reuse license and consent/IRB determination; adjudicate early-stage/benign-mass endpoints and non-ovarian disease exclusions; resolve donor/visit overlap across NCC series and training cohort. Then preregister frozen validation separately from harmonization/adaptation or de novo baseline development. No downloads or model development performed.
