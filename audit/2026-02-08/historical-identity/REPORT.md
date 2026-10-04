# Historical identity / probe / provenance — finalized static audit

All eleven requested scripts completely read; SHA-256/byte/line coverage in COVERAGE.csv. No project execution/import, model fitting, split-search run, network, installation, original modification or deletion. No recursive audit discovery. Counts in ENDPOINTS.csv distinguish source assertions from observations. Source data absent: original 135/44 membership remains UNKNOWN.

## Findings with exact evidence

All source paths below are under scripts/; line ranges inclusive.

### High — tanh score is not probability
validate_historical_reproduction.py::parse_pmml:408–515 requires two tanh layers, seven hidden and two output neurons. forward_pass:518–571 applies x*scale+shift then tanh twice (533,550). It sorts neuron IDs, assumes output 21=cancer/22=BBC, and selects CA only when first activation exceeds second (ties BBC). It does not derive direction from PMML NeuralOutput expressions or calibrate/softmax outputs. Cancer activation is a ranking score, not calibrated cancer risk.

audit_original_pmml_model.py::extract_neural_network:416–665 records output attributes but omits NeuralOutput child expressions; format_report:820–1108 and main:1115–1271 discuss probability without computing it. compare_pmml_to_published_table2.py::published_label:192–204 maps CA=1/BBC=0; main:211–532 maps cancer by substring. validate_pmml_against_table2_by_id.py module uses exact class-name mapping. These compare published NN predictions, not diagnoses.

### High — identity assumptions are not source proof
compare_pmml_to_published_table2.py::main:243 takes prediction-file historical_ids and :364 assigns them directly to publication rows. A one-to-one merge validates its constructed mapping, not independent biological linkage.

validate_pmml_against_table2_by_id.py::table2_id_to_actual_sample_id:163–226 maps numeric publication IDs <=55 to prediction row ID-1, 57–120 to ID-2, and passes e IDs through. Numeric linkage still depends on prediction ordering; set agreement cannot prove that ordering. Module truncates publication records to first 179.

validate_historical_reproduction.py::recover_historical_ids:84–243 requires 120 cohort-2 IDs and 60 Dataset-2 e IDs, excludes numeric position 56 asserted to be 81195, and produces 119+60=179. load_dataset2:247–405 checks e-tail order but uses reconstructed numeric order. parse_table2:574–674 returns publication IDs; main:677–822 joins those directly to actual cohort IDs without publication-number translation. Its conditional SUCCESS text does not establish successful execution or identity. No original training/test membership is recovered.

### Moderate — preprocessing, mature arms and duplicates
Dataset-2 reproduction uses supplementary values, not raw TPM, followed by PMML scale/shift (forward_pass:518–571). Historical batch adjustment/log10 provenance is not reconstructed from raw data by these scripts.

extract_original_14mirna_model.py::extract_pmml_features:139–373 searches all miRNA-like XML text/named fields, not exclusively one network's ordered active inputs; it does not enforce exactly 14 inputs. normalize_mirna_name:78–96 normalizes case/dashes/whitespace while retaining arm suffixes. audit_features:416–510 lists normalized-name collisions but does not resolve or aggregate duplicate probes. validate_historical_reproduction.py::load_dataset2:377 overwrites available[feature] on duplicate exact feature rows; numeric coercion and missing-value rejection do not detect duplicate-row provenance.

map_historical_mirnas_to_mimat.py::norm_name:29–32 removes punctuation other than hyphens, including historical star notation. Module :83–89 reads first header token as name and second as accession and overwrites duplicate normalized keys. Its own example comment has accession first, so actual downloaded header format must be verified, not inferred from the comment. CURRENT URL is mutable; no annotation release/hash/sequence equivalence or historical-arm proof is saved. No mature-arm substitution or duplicate-probe identity is established locally.

### High — compatible split search is exploratory, not holdout validation
search_possible_historical_test_splits.py:60 sets 20,000 trials; :369–375 targets Control4/Benign11/Borderline5/early12/late12, total44 (24 cancer/20 others). :397 seeds RNG94533; :401–459 samples within strata, computes AUC from true_label/cancer_activation and records selected IDs; :463–466 sorts by distance to .90. get_stratum:256–333 assigns missing-stage E1056 to early stage using a publication assertion not independently opened here. The script explicitly disclaims original-seed/split recovery. Counts plus matching AUC cannot prove original membership or independent-test accuracy.

None of the eleven scripts implements training/CV or an independent original holdout evaluation. Published prediction agreement is a different endpoint from disease discrimination; sampling subsets of existing historical predictions is not independent validation.

### Moderate — search coverage and provenance limits
search_original_44_patient_split.py::search_text_file:82–103 and inspect_xlsx:139–200 print keyword hits; module searches project text excluding reports, potentially rediscovering generated assertions. search_original_split_evidence.py::print_context:56–100 caps eight contexts; main:151–199 considers DOCX/XML only. search_split_evidence_clean.py::print_context:61–110 caps five contexts per keyword; XML scope is supplement-directory glob rather than the separate code archive. Absence of printed hits would not prove absence of membership evidence. These scripts warn that counts alone are insufficient, but do not produce a patient membership ledger.

inspect_original_elife_dataset.py::inspect_excel:127–273 previews first25/100 rows; inspect_code:276–386 caps relevant lines at120. download:37–67 and extract_zip:70–94 reuse existing artifacts without checksum/version validation. Nothing was executed in this audit.

## Local availability and saved reporting status

Bounded existence checks found eight candidates absent both repository-relative and configured Path.home()/ovarian-ai roots: data/external/elife-28932-code1/Source-code1.xml; supplementary files 1.docx and 6.xlsx under data/external/elife_28932; data/raw/GSE94533_Processed_file_Cohort1.xlsx and Cohort2.xlsx; reports/tables/original_historical_clean_predictions.csv; data/external/mirbase_mature.fa; data/processed/gse94533_metadata.csv. No patient rows/artifact contents were inspected. Actual counts, duplicate probes, sequences, agreement and split are unknown.

Intended saves only: original_pmml_deep_audit.json/.txt (audit_original_pmml_model.py::main:1115–1271); original_14_mirna_features.json and feature_audit.json (extract_original_14mirna_model.py::main:517–591); historical_mirna_to_mimat.csv (mapping module:127); pmml_vs_published_table2_comparison.csv (comparison main:525); possible_historical_test_splits.csv (search module:473); historical_pmml_vs_published_table2.csv (reproduction main:799); pmml_vs_table2_by_actual_sample_id.csv (ID validator module:568). Writing calls and conditional success messages are not evidence those outputs were generated successfully. Console evidence searches do not establish durable original membership provenance.

## PDF references

Existing audit/2026-02-08/pdf/extracted-pages.txt completely read in two bounded blocks; rendered PDF visuals not reviewed. Page4 identifies Dataset2 as same-patient batch-adjusted log10/192-miRNA historical reconstruction and PMML/Table2 as reference artifacts, not new independent rows. Page18 warns original135/44 membership is not recovered and historical >80% independent-test claims should not be presented as validated: consistent with this review. Page19 requests portable paths and canonical results index; all eleven reviewed scripts use hardcoded home-root paths. Page23 recommends other circulating-miRNA benchmarks, calibration and prospective harmonized collection. Recommendations are not evidence those activities were completed; PDF text cannot replace missing patient-level artifacts.

## Disposition

Keep reconstruction, published-prediction agreement, metric-compatible exploratory splits and independent validation as separate claims. Original split cannot be inferred from matching metrics. No repairs performed. evidence.py::source_evidence is a tiny reusable pure bytes-to-SHA256/numbered-excerpt capability in this area root.
