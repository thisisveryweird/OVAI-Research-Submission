# Data sources draft — reuse/consent/IRB unknown

Availability is not permission. No rights or exemption inferred from public access. Document original article/accession/version, acquisition date/hash, specimen/units, participant consent, allowed reuse/commercial use, donor/visit crosswalk and institutional determination before use. Author/contact/access routes remain proposals.

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


# Evidence-based data-source addendum (not permission or approval)

| Source | Citation/accession in local evidence | Local evidence and reuse status |
|---|---|---|
| GSE211692 | PMID 36426871; DOI 10.1093/jncics/pkac080; PMC9825310; GPL21263 | Original family SOFT fully parsed offline. Expression matrix cited in report is absent. Article citation recorded, not bibliographically authenticated against a local primary article. Public GEO metadata is not proof of unrestricted redistribution or donor independence. |
| GSE113486 | GEO GSE113486; PMID 30382619; GPL21263 | Original SOFT gives 972 serum sample records. Article/license/consent terms not locally authenticated. |
| Gockley 2022 | Local Gockley reports identify human serum cohort; DataSheet_2.xlsx and DataSheet_3.zip | Workbook and three PMML XML members inspected without scoring. Full article citation, acquisition chain, consent, reuse license and selected-network identity need primary-source verification. Folder placement under 2024 is not source identity. |
| 2024 adnexal mass | DOI 10.1016/j.ygyno.2024.08.008; reported PMCID PMC13283006 | NIHMS2021079-supplement-1.xlsx primary container inspected. Reported serum label conflicts with primary-article serum/plasma wording identified in the review requirement: exact article passage/correction and sample-level specimen map are not supplied locally, so unresolved, not adjudicated. No inference from title or workbook filename resolves specimen. |
| Danish pelvic mass 380 | Oliveira et al.; DOI 10.1371/journal.pone.0225249; PMC6860451; reported Figshare 10322942 | Local 48-assay transcription and aggregate/access reports only. No original patient matrix, supplement sequences or article authenticated here. Publication availability statement is not individual access approval. Author-contact route is source-reported; no request sent. |
| JanusRNA | Stawiski et al.; DOI 10.3389/fonc.2024.1389066; PMC11231195 | Local summary only; reported Federated EGA Norway/DAC access. No approved access, consent scope, data-use agreement or patient roster in reviewed evidence. |
| GSE94533 / historical eLife | Reference training/historical cohort in local research reports | Same underlying historical samples are not an independent cohort. Citation/reuse conditions and clinical donor linkage must be verified separately; previous finalized historical audits are comparisons only. |
| GSE106817 / GSE113740 | GEO accessions in transport/PDF reports | Saved metrics are not original patient evidence in this scope. Counts/endpoints vary; identify exact roster before use. |

**IRB/human-subjects determination, consent compatibility, citation correctness, source reuse/redistribution rights, and organizational approval are UNKNOWN where no documentary proof is supplied.** Public/de-identified does not itself authorize a particular secondary use. No legal or clinical exemption is certified here. PDF pp21–22 recommends documentation but its “very likely fine” wording is not a determination.

## Blocked roadmap — no applications or implementation

1. Obtain authenticated article/supplement citations and license/reuse terms; resolve 2024 specimen contradiction using the exact methods/correction and collection records.
2. Require documented institutional human-subjects determination and data-use permission before controlled-data access or redistribution. No dbGaP/EGA/author applications were submitted.
3. Require donor/visit crosswalks and biobank/site overlap attestations, not string-zero overlap. Preserve original split and distinguish model-selected Test from untouched Validation.
4. Require benign vs healthy/control phenotype adjudication, stage conflicts resolved by source records, timing and assay scale/batch metadata before clinical endpoint claims.
5. Obtain published frozen network/probability map, selected-network ID and consistent validation denominators before claiming reproduced published operating points. For 2024 the workbook supplies inputs, not a frozen fitted model.
6. Any prospective, HE4-linked, harmonization, calibration, threshold, fusion or longitudinal work in PDF pp8–9,20–23 is a proposal requiring separate approval and protocol. Nothing was applied, fitted, downloaded or patched in this audit.


2024 serum/plasma conflict requires authenticated article/correction and collection map; HE4 cells elsewhere do not establish GSE94533 linkage. GSE211692 SOFT metadata is available, expression matrix/control roster not established; Danish/Janus access remains source-reported. Historical primary workbooks/split/processing absent. All owner/IRB/regulatory approvals pending. No new acquisition performed.
