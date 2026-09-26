# GSE211692 OVAI External-Validation Audit

Scope: research-only GEO, publication, platform, sample-level expression, label, stage, feature, clinical-variable, and identifier-overlap audit. No model was trained, refit, tuned, recalibrated, or changed; no cohorts were pooled; `production_artifacts/` was read only.

## Accession and Publication

- GEO accession `GSE211692`, title “Prediction of tissue-of-origin of early-stage cancers using serum miRNomes”; public series SOFT record. PMID `36426871`; DOI `10.1093/jncics/pkac080`; PMC `PMC9825310`.
- GEO series type: non-coding RNA profiling by array; platform `GPL21263`, “3D-Gene Human miRNA V21_1.0.0”, commercial spotted oligonucleotide array manufactured by Toray. Platform annotation has 2,565 mature-miRNA feature rows, including MIMAT accessions. The title's V21_1.0.0 is the product version; an explicit miRBase release is not stated in the GEO platform fields.
- Specimen: individual human serum. GEO extract protocol states total RNA from 300 uL serum per sample, extracted with the 3D-Gene RNA extraction reagent. One-color Cy5 3D-Gene miRNA labeling, 16-hour hybridization at 32 C, Toray 3D-Gene scanner, Extraction software. Publication says samples were preserved in the NCC Biobank; storage duration and pre-treatment status are not established by the inspected GEO records.
- The publication states it used a 3D-Gene v21 microarray and 3 internal-control miRNAs (miR-149-3p, miR-2861, miR-4463) for normalization. GEO describes negative-control background subtraction; negative/undetected signal replacement at 0.1 on a base-2 logarithm scale; and normalized signal intensity by internal controls. The processed matrix contains signed normalized values (both negative and positive), not sequencing counts, TPM, RPM, or CPM. Exact final transform order/units are not specified sufficiently to equate these values to the GSE94533 TPM scale.

## Patient-Level Matrix and Integrity

- Authoritative sources: [NCBI GEO Series GSE211692](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211692), [official family-SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211692/soft/GSE211692_family.soft.gz), [processed expression matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211692/suppl/GSE211692_processed_data_matrix.txt.gz), [GPL21263 platform record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL21263), [PubMed PMID 36426871](https://pubmed.ncbi.nlm.nih.gov/36426871/), and [PMC full text PMC9825310](https://pmc.ncbi.nlm.nih.gov/articles/PMC9825310/). The official GEO family-SOFT and processed matrix passed gzip integrity checks.
- Processed matrix local path: `research/cohort_expansion/data/raw/GSE211692/GSE211692_processed_data_matrix.txt.gz`; 230724363 bytes; SHA256 `5c3445626c426129c6d98eacabaacadd649325925ccc2e8ebd1337e9fef9b70e`.
- Family SOFT local path: `research/cohort_expansion/data/raw/GSE211692/GSE211692_family.soft.gz`; 741790 bytes; SHA256 `46ab13eaa51873e4a0a30ca8c46d9ea3b41c528eed5c2fe6fa99a4fdc926c1ff`.
- GEO family-SOFT has 16190 unique GSM records and 16190 unique sample titles. The matrix has 16190 sample columns plus `ID_REF`; every column title maps one-to-one to a SOFT sample title, with no unmatched columns or records.
- The processed matrix contains 2565 rows/features; all feature IDs are unique and all 41,527,350 expression cells are populated. Exact duplicate sample expression vectors: 0.
- GEO provides one GSM/sample record and one per-sample raw supplementary file for each processed column. No technical-replicate field or duplicate sample title is present. The `Cell1`-`Cell4` strings in raw filenames are acquisition-file labels; GEO does not document them as biological replicates. No separate subject/visit identifier permits definitive exclusion of multiple samples from one individual. Thus expression is sample-level; one-to-one patient independence cannot be proven from public IDs alone.

## Reconstructed Cohort Labels

- Matrix/GEO total: 16,190 serum sample records. GEO design totals reconcile: 9,921 cancer diagnoses across 13 cancer types, 602 clearly labeled benign-disease samples, 24 mixed extraparenchymal-brain-tumor/benign-brain-disease samples, and 5,643 `no cancer` controls, total 16,190. All GEO sample records have a nonempty disease-state field; no explicit borderline-tumor label occurs.
| GEO disease-state label | Count |
|---|---:|
| `benign disease in the bone and soft tissue` | 313 |
| `benign disease in the breast` | 31 |
| `benign disease in the ovary` | 28 |
| `benign disease in the prostate` | 230 |
| `biliary tract cancer` | 402 |
| `bladder cancer` | 399 |
| `bone and soft tissue sarcomas` | 299 |
| `breast cancer` | 675 |
| `colorectal cancer` | 1596 |
| `esophageal squamous cell cancer` | 566 |
| `extraparenchymal brain tumor and benign disease in the brain` | 24 |
| `gastric cancer` | 1418 |
| `hepatocellular cancer` | 348 |
| `intraparenchymal brain tumors` | 241 |
| `lung cancer` | 1699 |
| `no cancer` | 5643 |
| `ovarian cancer` | 400 |
| `pancreatic cancer` | 851 |
| `prostate cancer` | 1027 |

- Ovarian-specific GEO labels: 400 `ovarian cancer`; 28 `benign disease in the ovary`; 5,643 `no cancer` controls (3,175 female, 2,468 male). Other non-cancer disease labels include 574 clearly benign breast/prostate/bone-soft-tissue samples and 24 mixed brain-tumor/benign-brain-disease records; none are silently treated as ovarian benign masses or healthy controls. Other cancer types remain separate.
- The original paper says only 5 benign-disease control samples were available for its cancer-detection analysis and therefore focuses on tissue-of-origin prediction. This differs from the 28 ovarian-benign labels in the broader GEO matrix. The 28 GEO-labeled samples are a metadata-defined candidate comparator, but should not be represented as the paper's evaluated benign-control subset without a sample inclusion crosswalk.
- Stage I/II ovarian cancer vs GEO-labeled benign ovarian disease is reconstructible from the broader matrix as 168 vs 28, but benign subtype/mass details are not supplied. Stage I/II vs benign ovarian disease plus all `no cancer` controls is 168 vs 5,671; these controls include men and should not be silently treated as age/sex-matched women. Miscellaneous cancers and benign diseases are excluded from these endpoints, not relabeled.
- Reported recent subset `168 Stage I/II + 200 non-cancer controls`: 168 cases are exactly reconstructible from ovarian `Stage: 1` (131) + `Stage: 2` (37). The 200 controls are not a distinct GEO label or published split; 5,643 rows are `no cancer`, 28 are benign ovarian disease, and GEO does not flag a particular 200. The paper's aggregate training/validation sets are multicancer and do not identify an ovarian 200-control roster. Therefore the 168/200 subset cannot be independently reconstructed as a uniquely defined cohort.

## Ovarian Stage Audit

| Ovarian cancer stage field | Count |
|---|---:|
| `0` | 0 |
| `1` | 131 |
| `2` | 37 |
| `3` | 214 |
| `4` | 0 |
| `NA` | 18 |
| `missing` | 0 |
- Stage I/II count is 168 (Stage I=131, Stage II=37); this matches 168. Stage III=214, Stage IV=0, unknown=`NA`=18. The GEO phenotype is simply `Stage: 1`, `2`, `3`, `NA`; neither the sample field nor inspected article labels the ovarian stages explicitly as FIGO. A FIGO system cannot be independently certified.
- GEO has no borderline ovarian tumor label. The ovarian record set has no stage 0 or stage 4 entries; do not infer FIGO or cohort-specific inclusion from the general paper's multicancer `stages 0-II` language.

## Exact Feature Compatibility

- Current OVAI signature: 5/5 exact mature identifiers. All five are present for all matrix samples and map to distinct GPL MIMAT accessions. This is analyte availability only; it does not establish assay-scale compatibility.
- Historical 14-miRNA signature: 13/14 exact identifiers. The only absent exact label is `hsa-miR-203a`; GPL/matrix contain `hsa-miR-203a-3p` (MIMAT0000264) and `hsa-miR-203a-5p` (MIMAT0031890). Without a MIMAT accession or arm in the requested legacy identifier, this is recorded as an alias candidate and not counted exact.
- Exact arm separation: current `hsa-miR-32-3p` is MIMAT0004505; historical `hsa-miR-32-5p` is MIMAT0000090. They are different mature products. No arm substitution was made.
| Signature | Requested identifier | Exact found | MIMAT | Negative values / min..max | Alias note (not counted exact) |
|---|---|---:|---|---|---|
| current_5 | `hsa-miR-200b-3p` | True | MIMAT0000318 | 5961 / -5.981344809..12.79332178 | none identified |
| current_5 | `hsa-miR-296-3p` | True | MIMAT0004679 | 3 / -1.341673501..12.98358413 | none identified |
| current_5 | `hsa-miR-6826-5p` | True | MIMAT0027552 | 5 / -2.079342976..11.83318372 | none identified |
| current_5 | `hsa-miR-32-3p` | True | MIMAT0004505 | 6855 / -6.326731167..10.57543835 | hsa-miR-32-5p is separately annotated as MIMAT0000090; opposite arm, not alias. |
| current_5 | `hsa-miR-1255a` | True | MIMAT0005906 | 9494 / -6.326731167..11.99040564 | none identified |
| historical_14 | `hsa-miR-23b-3p` | True | MIMAT0000418 | 1078 / -5.06805668..13.90787827 | none identified |
| historical_14 | `hsa-miR-29a-3p` | True | MIMAT0000086 | 6159 / -6.307353452..11.93653438 | none identified |
| historical_14 | `hsa-miR-32-5p` | True | MIMAT0000090 | 5837 / -6.326731167..11.2352491 | none identified |
| historical_14 | `hsa-miR-92a-3p` | True | MIMAT0000092 | 117 / -3.26072416..12.18131103 | none identified |
| historical_14 | `hsa-miR-150-5p` | True | MIMAT0000451 | 2183 / -5.06805668..11.83845283 | none identified |
| historical_14 | `hsa-miR-200a-3p` | True | MIMAT0000682 | 6934 / -6.326731167..11.18210986 | none identified |
| historical_14 | `hsa-miR-200c-3p` | True | MIMAT0000617 | 5381 / -5.936300843..11.8358106 | none identified |
| historical_14 | `hsa-miR-203a` | False | NA |  / .. | Matrix/platform contains hsa-miR-203a-3p (MIMAT0000264) and hsa-miR-203a-5p (MIMAT0031890); legacy label is not exact and arm identity is not stated in requested name. |
| historical_14 | `hsa-miR-320c` | True | MIMAT0005793 | 949 / -5.06805668..11.9616454 | none identified |
| historical_14 | `hsa-miR-320d` | True | MIMAT0006764 | 2388 / -6.307353452..11.65385964 | none identified |
| historical_14 | `hsa-miR-335-5p` | True | MIMAT0000765 | 7187 / -6.326731167..11.33265775 | none identified |
| historical_14 | `hsa-miR-450b-5p` | True | MIMAT0004909 | 7187 / -6.326731167..10.42375513 | none identified |
| historical_14 | `hsa-miR-1246` | True | MIMAT0005898 | 1822 / -6.307353452..14.35462943 | none identified |
| historical_14 | `hsa-miR-1307-5p` | True | MIMAT0022727 | 2251 / -5.06805668..9.97351544 | none identified |

## Clinical-Variable Compatibility

- Age: present for all 400/400 ovarian cancers, 28/28 ovarian-benign records, and 5641/5643 no-cancer controls. GEO gives numeric ages but does not explicitly state units in the characteristic key; treatment timing is not recorded.
- CA-125: no CA-125/CA125 characteristic or matrix field is present for the GSE211692 samples; availability is 0/16,190. Therefore the full frozen production model cannot be evaluated on this cohort, even though its frozen imputation defaults exist. The local GSE94533 evaluation data have age for 178/179 patients and CA-125 for 119/179; units/timing must be treated per that cohort's source, not transferred to GSE211692.
- The production feature manifest is five miRNAs + age + CA-125 + two interactions (nine predictors). No separate frozen five-miRNA-only production estimator is in the manifest. Therefore exact 5/5 assay coverage does not create a frozen molecular-only model, and deriving one by dropping clinical terms would alter the frozen model. No prediction was attempted.
- Frozen production preprocessing applies `log1p(max(value, 0))`. Negative array values that would be floored to zero occur as follows in the current five features: hsa-miR-200b-3p=5961/16190; hsa-miR-296-3p=3/16190; hsa-miR-6826-5p=5/16190; hsa-miR-32-3p=6855/16190; hsa-miR-1255a=9494/16190. This is direct evidence of an input-distribution/transform mismatch; no values were changed or scored.

## Comparison With GSE94533 and Independence

- Both are serum, but GSE211692 is a 3D-Gene spotted oligonucleotide array with internal-control-normalized signal values; GSE94533 is GEO-annotated as high-throughput sequencing. The GSE94533 paper describes Exiqon serum RNA extraction, NEBNext small-RNA library preparation, Illumina NextSeq 500 single-end 50-nt reads, and mapping to miRBase 20; local OVAI documentation uses TPM followed by `log1p(max(value, 0))`. GPL21263 provides mature names/MIMAT IDs but not an explicit miRBase release. The feature names overlap, but assay chemistry, library method, measurement scale, normalization, and missing-value behavior differ. Applying the production transform to signed array values would clamp all negative signals to zero and is not an equivalence-preserving harmonization.
- GSE94533 GEO series reports GPL18573 and 180 serum samples; the local curated endpoint has 179 patients (52 Stage I/II cancers, 61 benign/control) and 2,582 features. GSE211692 uses GPL21263 and 2,565 features. Age/CA-125 coverage is incomplete in GSE94533 (178/179 age, 119/179 CA-125) and absent for CA-125 in GSE211692.
- Sample-ID audit: no exact or responsibly normalized overlap was detected between GSE211692 sample titles and GSE94533, GSE106817, GSE113486, Gockley 2022, or the 2024 FirePlex workbook identifiers. These cohorts use different identifier namespaces, so **no exact overlap detected, but independence cannot be proven from available identifiers**.
| Comparator | Roster n | Exact overlap | Normalized overlap | Independence proven? |
|---|---:|---:|---:|---:|
| GSE94533 | 179 | 0 | 0 | no |
| GSE106817 | 4046 | 0 | 0 | no |
| GSE113486 | 972 | 0 | 0 | no |
| Gockley 2022 | 275 | 0 | 0 | no |
| 2024 adnexal FirePlex | 678 | 0 | 0 | no |
- GSE106817 and GSE113486 also use GPL21263; GSE94533 uses a different sequencing platform. GSE211692's own publication is a multicancer tissue-of-origin classifier, not a dedicated ovarian cancer vs benign adnexal-mass diagnostic evaluation. Its validation splits are aggregate and not recoverable per ovarian sample in GEO metadata.

## Prevalence and Eligibility Decision

- Study composition prevalence: all cancer types 9,921/16,190 = 61.28%; ovarian cancer 400/16,190 = 2.47%; Stage I/II ovarian cancer 168/16,190 = 1.04%. These are enriched case-control study proportions, not population prevalence.
- Among metadata-defined ovarian Stage I/II plus `no cancer` rows, apparent case fraction is 168/(168+5,643) = 2.89%; adding 28 ovarian-benign labels gives 168/(168+5,643+28) = 2.88%. A deliberately selected 200-control subset would yield a different artificial fraction and cannot be attributed to the publication from GEO.
- Frozen current OVAI external validation: **not suitable**. The matrix is sample-level and includes all five exact feature names, but platform values are normalized microarray signal rather than GSE94533 TPM, signed values interact badly with the frozen `max(x,0)` + `log1p` input transform, the model requires CA-125 which is absent, and no frozen molecular-only model exists. The clinically preferred benign-mass-only endpoint is only 168 vs 28 by broad labels, with the original paper reporting only five benign-disease controls in its analyzed collection and no sample inclusion crosswalk.
- Appropriate present use: cohort characterization and hypothesis-generating biological replication of exact mature miRNA associations, with platform and endpoint kept explicit. Any independent cross-platform model evaluation would require a prespecified, justified harmonization and complete clinical/endpoint audit; no such processing or prediction is performed here.

No model fitting, feature selection, tuning, recalibration, pooling, or canonical artifact modification was performed. The 1-GB raw GEO tar was not downloaded; the official processed matrix and family-SOFT metadata were sufficient for this audit.
