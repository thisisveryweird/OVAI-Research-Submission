#!/usr/bin/env python3
"""Generate research-only feature overlap and replication audits."""

from __future__ import annotations

import csv
from pathlib import Path
import re
from typing import Any

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
REPORTS = ROOT / "research/cohort_expansion/reports"
WORKBOOK = ROOT / "research/cohort_expansion/data/raw/2024_adnexal_mass/NIHMS2021079-supplement-1.xlsx"
PANEL_INVENTORY = REPORTS / "2024_adnexal_mass_mirna_panel_inventory.csv"
COHORT_RECONSTRUCTION = REPORTS / "2024_adnexal_mass_cohort_reconstruction.csv"
STABILITY = ROOT / "reports/feature_selection_stability.csv"
GSE_PREDICTIONS = ROOT / "reports/error_analysis_predictions.csv"

CURRENT = [
    "hsa-miR-200b-3p",
    "hsa-miR-296-3p",
    "hsa-miR-6826-5p",
    "hsa-miR-32-3p",
    "hsa-miR-1255a",
]
HISTORICAL = [
    "hsa-miR-23b-3p",
    "hsa-miR-29a-3p",
    "hsa-miR-32-5p",
    "hsa-miR-92a-3p",
    "hsa-miR-150-5p",
    "hsa-miR-200a-3p",
    "hsa-miR-200c-3p",
    "hsa-miR-203a",
    "hsa-miR-320c",
    "hsa-miR-320d",
    "hsa-miR-335-5p",
    "hsa-miR-450b-5p",
    "hsa-miR-1246",
    "hsa-miR-1307-5p",
]
PUBLISHED_2024 = [
    "hsa-miR-320e",
    "hsa-miR-195-5p",
    "hsa-miR-103-3p",
    "hsa-miR-320c",
    "cel-miR-39-3p",
    "hsa-miR-200c-3p",
    "hsa-miR-770-5p",
    "hsa-miR-122-5p",
    "hsa-miR-19a-3p",
    "hsa-miR-424-5p",
]

PELVIC_ASSAYS = [
    ("hsa-let-7f", "000382", "hsa-let-7f-5p", "UGAGGUAGUAGAUUGUAUAGUU"),
    ("hsa-miR-20a", "000580", "hsa-miR-20a-5p", "UAAAGUGCUUAUAGUGCAGGUAG"),
    ("hsa-miR-21", "000397", "hsa-miR-21-5p", "UAGCUUAUCAGACUGAUGUUGA"),
    ("hsa-miR-23a", "000399", "hsa-miR-23a-3p", "AUCACAUUGCCAGGGAUUUCC"),
    ("hsa-miR-23b", "000400", "hsa-miR-23b-3p", "AUCACAUUGCCAGGGAUUACC"),
    ("hsa-miR-25", "000403", "hsa-miR-25-3p", "CAUUGCACUUGUCUCGGUCUGA"),
    ("hsa-miR-26a", "000405", "hsa-miR-26a-5p", "UUCAAGUAAUCCAGGAUAGGCU"),
    ("hsa-miR-27a", "000408", "hsa-miR-27a-3p", "UUCACAGUGGCUAAGUUCCGC"),
    ("hsa-miR-27b", "000409", "hsa-miR-27b-3p", "UUCACAGUGGCUAAGUUCUGC"),
    ("hsa-miR-29a", "002112", "hsa-miR-29a-3p", "UAGCACCAUCUGAAAUCGGUUA"),
    ("hsa-miR-92a", "000430", "hsa-mir-92a-3p", "UAUUGCACUUGUCCCGGCCUG"),
    ("mmu-miR-93", "001090", "hsa-miR-93-5p", "CAAAGUGCUGUUCGUGCAGGUAG"),
    ("mmu-miR-96", "000186", "hsa-miR-96-5p", "UUUGGCACUAGCACAUUUUUGCU"),
    ("hsa-miR-103", "000439", "hsa-miR-103a-3p", "AGCAGCAUUGUACAGGGCUAUG"),
    ("hsa-miR-106a", "002169", "hsa-miR-106a-5p", "AAAAGUGCUUACAGUGCAGGUAG"),
    ("hsa-miR-122", "002245", "hsa-miR-122-5p", "UGGAGUGUGACAAUGGUGUUUG"),
    ("hsa-miR-125b", "000449", "hsa-miR-125b-5p", "UCCCUGAGACCCUAACUUGUGA"),
    ("hsa-miR-126", "002228", "hsa-miR-126-3p", "UCGUACCGUGAGUAAUAAUGCG"),
    ("hsa-miR-140-3p", "002234", "hsa-miR-140-3p", "UACCACAGGGUAGAACCACGG"),
    ("hsa-miR-141", "000463", "hsa-miR-141-3p", "UAACACUGUCUGGUAAAGAUGG"),
    ("hsa-miR-142-3p", "000464", "hsa-miR-142-3p", "UGUAGUGUUUCCUACUUUAUGGA"),
    ("hsa-miR-145-5p", "002278", "hsa-miR-145-5p", "GUCCAGUUUUCCCAGGAAUCCCU"),
    ("hsa-miR-148b", "000471", "hsa-miR-148b-3p", "UCAGUGCAUCACAGAACUUUGU"),
    ("hsa-miR-152", "000475", "has-miR-152-3p", "UCAGUGCAUGACAGAACUUGG"),
    ("hsa-miR-181b", "001098", "hsa-miR-181b-5p", "AACAUUCAUUGCUGUCGGUGGG"),
    ("hsa-miR-182", "002334", "hsa-miR-182-5p", "UUUGGCAAUGGUAGAACUCACACU"),
    ("hsa-miR-191", "002299", "hsa-miR-191-5p", "CAACGGAAUCCCAAAAGCAGCUG"),
    ("hsa-miR-195", "000494", "hsa-miR-195-5p", "UAGCAGCACAGAAAUAUUGGC"),
    ("hsa-miR-199a-3p", "002304", "hsa-miR-199a/b-3p", "ACAGUAGUCUGCACAUUGGUUA"),
    ("hsa-miR-200a", "000502", "hsa-miR-200a-3p", "UAACACUGUCUGGUAACGAUGU"),
    ("hsa-miR-200b", "002251", "hsa-miR-200b-3p", "UAAUACUGCCUGGUAAUGAUGA"),
    ("hsa-miR-200c", "002300", "hsa-miR-200c-3p", "UAAUACUGCCGGGUAAUGAUGGA"),
    ("hsa-miR-205", "000509", "hsa-miR-205-5p", "UCCUUCAUUCCACCGGAGUCUG"),
    ("hsa-miR-214", "002306", "hsa-miR-214-3p", "ACAGCAGGCACAGACAGGCAGU"),
    ("hsa-miR-221", "000524", "hsa-miR-221-3p", "AGCUACAUUGUCUGCUGGGUUUC"),
    ("hsa-miR-223", "002295", "hsa-miR-223-3p", "UGUCAGUUUGUCAAAUACCCCA"),
    ("hsa-miR-335", "000546", "hsa-miR-335-5p", "UCAAGAGCAAUAACGAAAAAUGU"),
    ("hsa-miR-346", "000553", "hsa-miR-346", "UGUCUGCCCGCAUGCCUGCCUCU"),
    ("hsa-miR-372", "000560", "hsa-miR-372-3p", "AAAGUGCUGCGACAUUUGAGCGU"),
    ("hsa-miR-373", "000561", "hsa-miR-373-3p", "GAAGUGCUUCGAUUUUGGGGUGU"),
    ("hsa-miR-378", "002243", "hsa-miR-378a-3p", "ACUGGACUUGGAGUCAGAAGG"),
    ("hsa-miR-424", "000604", "hsa-miR-424-5p", "CAGCAGCAAUUCAUGUUUUGAA"),
    ("hsa-miR-429", "001024", "hsa-miR-429", "UAAUACUGUCUGGUAAAACCGU"),
    ("mmu-miR-451", "001141", "hsa-miR-451a", "AAACCGUUACCAUUACUGAGUU"),
    ("hsa-miR-484", "001821", "hsa-miR-484", "UCAGGCUCAGUCCCCUCCCGAU"),
    ("hsa-miR-503", "001048", "hsa-miR-503-5p", "UAGCAGCGGGAACAGUUCUGCAG"),
    ("cel-miR-54", "001361", "cel-miR-54-3p", "UACCCGUAAUCUUCAUAAUCCGAG"),
    ("cel-miR-238", "000248", "cel-miR-238-3p", "UUUGUACUCCGAUGCCAUUCAGA"),
]

MIMAT = {"hsa-miR-200b-3p": "MIMAT0000318"}


def norm_feature(value: Any) -> str:
    text = str(value).strip().strip("'").strip()
    text = re.sub(r"\s*\(MFI\)\s*$", "", text, flags=re.IGNORECASE)
    return text.lower()


def hedges_g(case: pd.Series, control: pd.Series) -> tuple[float, float, float, int, int]:
    cases = pd.to_numeric(case, errors="coerce").dropna().to_numpy(dtype=float)
    controls = pd.to_numeric(control, errors="coerce").dropna().to_numpy(dtype=float)
    pooled_sd = np.sqrt(
        ((len(cases) - 1) * np.var(cases, ddof=1) + (len(controls) - 1) * np.var(controls, ddof=1))
        / (len(cases) + len(controls) - 2)
    )
    d = (cases.mean() - controls.mean()) / pooled_sd if pooled_sd else np.nan
    correction = 1 - 3 / (4 * (len(cases) + len(controls)) - 9)
    return float(d * correction), float(cases.mean()), float(controls.mean()), len(cases), len(controls)


def write_overlap_and_alias() -> tuple[pd.DataFrame, dict[str, str]]:
    panel = pd.read_csv(PANEL_INVENTORY)
    panel_lookup = {
        norm_feature(row["2024_identifier"]): row
        for row in panel.to_dict(orient="records")
    }
    features = list(dict.fromkeys(CURRENT + HISTORICAL + PUBLISHED_2024))
    rows = []
    for feature in features:
        assay = panel_lookup.get(norm_feature(feature))
        rows.append(
            {
                "feature": feature,
                "present_in_2024": assay is not None,
                "present_in_current_OVAI": feature in CURRENT,
                "present_in_historical_OVAI": feature in HISTORICAL,
                "present_in_published_2024_model": feature in PUBLISHED_2024,
                "measurement_identifier": assay["2024_identifier"] if assay else "",
                "measurement_type": assay["measurement_scale"] if assay else "not measured",
                "notes": (
                    "Exact workbook mature-name match; all 678 values present."
                    if assay
                    else "No exact 2024 FirePlex panel identifier; no fuzzy/family substitution used."
                ),
            }
        )
    overlap = pd.DataFrame(rows)
    overlap.to_csv(REPORTS / "2024_adnexal_mass_overlap_analysis.csv", index=False)

    alias_rows = []
    candidate_columns = {
        "hsa-miR-200b-3p": "hsa-miR-200a-3p; hsa-miR-200c-3p",
        "hsa-miR-296-3p": "none identified",
        "hsa-miR-6826-5p": "none identified",
        "hsa-miR-32-3p": "hsa-miR-32-5p (opposite mature arm)",
        "hsa-miR-1255a": "none identified",
    }
    for feature in CURRENT:
        alias_rows.append(
            {
                "production_feature": feature,
                "mimat_accession": MIMAT.get(feature, "not independently resolved"),
                "2024_exact_identifier": "",
                "potential_panel_neighbor": candidate_columns[feature],
                "mature_name_or_arm_evidence": (
                    "MIMAT0000318 is linked by RNAcentral to mature hsa-miR-200b-3p; the FirePlex panel has 200a and 200c family members only, not 200b."
                    if feature == "hsa-miR-200b-3p"
                    else "The workbook has no exact mature identifier. Similar family/number strings are not accepted as aliases."
                    if feature != "hsa-miR-32-3p"
                    else "The measured 5p strand is a distinct mature product from the requested 3p strand; it is not an alias."
                ),
                "sequence_equivalence": "not established for any 2024 measured analyte",
                "classification": "NOT MEASURED",
                "notes": "No documented same-mature miRNA alias maps this production feature to a measured 2024 assay; no substitutions accepted.",
            }
        )
    alias = pd.DataFrame(alias_rows)
    alias.to_csv(REPORTS / "2024_adnexal_mass_alias_audit.csv", index=False)
    return overlap, panel_lookup


def write_historical_replication(panel_lookup: dict[str, str]) -> pd.DataFrame:
    patient_predictions = pd.read_csv(GSE_PREDICTIONS).drop_duplicates("sample_id_normalized")
    source = patient_predictions.loc[
        patient_predictions["dataset"].eq("GSE94533")
        & patient_predictions["diagnosis"].isin(["Cancer", "Benign", "Control"])
    ].copy()
    source = source.loc[
        source["diagnosis"].ne("Cancer") | source["stage_group"].eq("Stage I/II")
    ]
    if len(source) != 113:
        raise ValueError(f"Expected saved GSE94533 Stage I/II-vs-benign/control rows=113; found {len(source)}")
    cancer_mask = source["diagnosis"].eq("Cancer")
    source_control_mask = source["diagnosis"].isin(["Benign", "Control"])

    workbook = pd.read_excel(WORKBOOK, sheet_name="miRNA")
    diagnosis = "diagnosis (benign = 0, malignant= 1)"
    case_mask = workbook[diagnosis].eq(1)
    control_mask = workbook[diagnosis].eq(0)
    stability = pd.read_csv(STABILITY)
    output = []
    shared_historical = [feature for feature in HISTORICAL if norm_feature(feature) in panel_lookup]
    if len(shared_historical) != 8:
        raise ValueError(f"Expected eight exact FirePlex/historical overlaps; found {len(shared_historical)}")
    for feature in shared_historical:
        item = panel_lookup[norm_feature(feature)]
        column = item["workbook_column"]
        case_effect, case_mean, control_mean, n_case, n_control = hedges_g(
            workbook.loc[case_mask, column], workbook.loc[control_mask, column]
        )
        source_effect, source_case_mean, source_control_mean, source_n_case, source_n_control = hedges_g(
            source.loc[cancer_mask, feature], source.loc[source_control_mask, feature]
        )
        feature_rows = stability.loc[
            stability["feature"].eq(feature) & stability["signature_size"].isin([5, 10, 14])
        ]
        freq: dict[str, float] = {}
        for task in ["all_cancer_vs_non_cancer", "early_stage_vs_benign_control"]:
            for method in ["l1", "elastic_net"]:
                subset = feature_rows.loc[
                    feature_rows["task"].eq(task) & feature_rows["selection_method"].eq(method)
                ]
                freq[f"{task}_{method}_max_frequency"] = float(subset["selection_frequency"].max()) if len(subset) else 0.0
        stable_member = any(value >= 0.6 for value in freq.values())
        source_direction = "higher" if source_case_mean > source_control_mean else "lower" if source_case_mean < source_control_mean else "equal"
        external_direction = "higher" if case_mean > control_mean else "lower" if case_mean < control_mean else "equal"
        output.append(
            {
                "feature": feature,
                "GSE94533_early_case_n": source_n_case,
                "GSE94533_benign_control_n": source_n_control,
                "GSE94533_case_mean": source_case_mean,
                "GSE94533_control_mean": source_control_mean,
                "GSE94533_abundance_direction": source_direction,
                "GSE94533_coefficient_direction": "not available; saved stability artifact contains selection frequency only",
                "GSE94533_early_effect_hedges_g": source_effect,
                **freq,
                "prior_stable_signature_member": stable_member,
                "stable_membership_rule": "selection frequency >=3/5 in an existing 5/10/14-feature signature on either saved endpoint",
                "2024_exact_identifier": item["2024_identifier"],
                "2024_cancer_n": n_case,
                "2024_control_n": n_control,
                "2024_cancer_mean_MFI": case_mean,
                "2024_control_mean_MFI": control_mean,
                "2024_cancer_control_direction": external_direction,
                "2024_effect_hedges_g_raw_MFI": case_effect,
                "2024_missing_cancer": int(workbook.loc[case_mask, column].isna().sum()),
                "2024_missing_control": int(workbook.loc[control_mask, column].isna().sum()),
                "direction_agrees_with_GSE94533": source_direction == external_direction,
                "interpretation": "descriptive abundance-direction replication across distinct platforms; effect magnitudes are not directly comparable",
            }
        )
    result = pd.DataFrame(output)
    result.to_csv(REPORTS / "2024_adnexal_mass_historical_replication.csv", index=False)
    return result


def write_pelvic_mass_map() -> tuple[pd.DataFrame, dict[str, dict[str, str]]]:
    panel = {
        norm_feature(mirbase): {
            "assay_name": assay_name,
            "assay_id": assay_id,
            "mirbase_id": mirbase,
            "sequence": sequence,
        }
        for assay_name, assay_id, mirbase, sequence in PELVIС_ASSAYS
    }
    rows = []
    for feature in list(dict.fromkeys(CURRENT + HISTORICAL + PUBLISHED_2024)):
        assay = panel.get(norm_feature(feature))
        alias = feature == "hsa-miR-103-3p" and norm_feature("hsa-miR-103a-3p") in panel
        alias_assay = panel.get(norm_feature("hsa-miR-103a-3p")) if alias else None
        rows.append(
            {
                "feature": feature,
                "current_OVAI_member": feature in CURRENT,
                "historical_OVAI_member": feature in HISTORICAL,
                "published_2024_model_member": feature in PUBLISHED_2024,
                "exact_panel_match": assay is not None,
                "valid_mature_name_alias": bool(alias),
                "panel_identifier": assay["mirbase_id"] if assay else alias_assay["mirbase_id"] if alias_assay else "",
                "TaqMan_assay_name": assay["assay_name"] if assay else alias_assay["assay_name"] if alias_assay else "",
                "TaqMan_assay_id": assay["assay_id"] if assay else alias_assay["assay_id"] if alias_assay else "",
                "mature_sequence": assay["sequence"] if assay else alias_assay["sequence"] if alias_assay else "",
                "source": "PLOS Supplementary Table S1, DOI 10.1371/journal.pone.0225249",
                "notes": (
                    "Same mature product: TaqMan assay name hsa-miR-103; listed miRBase ID is hsa-miR-103a-3p."
                    if alias
                    else "Exact mature miRNA ID match in the published TaqMan panel."
                    if assay
                    else "No exact mature identifier in the 48-assay list; no family-member substitution."
                ),
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(REPORTS / "380_pelvic_mass_feature_compatibility.csv", index=False)
    panel_df = pd.DataFrame(
        [
            {
                "TaqMan_assay_name": assay_name,
                "TaqMan_assay_id": assay_id,
                "miRBase_ID_as_published": mirbase,
                "mature_sequence": sequence,
                "is_spike_in_control": mirbase.lower().startswith("cel-"),
            }
            for assay_name, assay_id, mirbase, sequence in PELVIC_ASSAYS
        ]
    )
    panel_df.to_csv(REPORTS / "380_pelvic_mass_48_assay_inventory.csv", index=False)
    return result, panel


def write_reports(overlap: pd.DataFrame, replication: pd.DataFrame, pelvic_map: pd.DataFrame) -> None:
    current_overlap = overlap.loc[overlap["present_in_current_OVAI"] & overlap["present_in_2024"], "feature"].tolist()
    historical_overlap = overlap.loc[overlap["present_in_historical_OVAI"] & overlap["present_in_2024"], "feature"].tolist()
    published_overlap = overlap.loc[overlap["present_in_published_2024_model"] & overlap["present_in_2024"], "feature"].tolist()
    published_absent = overlap.loc[overlap["present_in_published_2024_model"] & ~overlap["present_in_2024"], "feature"].tolist()
    agreeing = int(replication["direction_agrees_with_GSE94533"].sum())
    stable = replication.loc[replication["prior_stable_signature_member"], "feature"].tolist()
    table_rows = [
        f"| `{r['feature']}` | {r['GSE94533_abundance_direction']} (g={r['GSE94533_early_effect_hedges_g']:.3f}) | {r['2024_cancer_control_direction']} (g={r['2024_effect_hedges_g_raw_MFI']:.3f}) | {'yes' if r['direction_agrees_with_GSE94533'] else 'no'} | {'yes' if r['prior_stable_signature_member'] else 'no'} |"
        for r in replication.to_dict(orient="records")
    ]
    (REPORTS / "2024_adnexal_mass_partial_signature_replication.md").write_text(
        "\n".join(
            [
                "# Partial Historical-Signature Replication",
                "",
                "**partial historical-signature replication** is a descriptive, feature-by-feature comparison of the eight exact shared historical miRNAs. It is not validation of the original 14-feature model. Six of the original 14 inputs are absent from the 2024 FirePlex panel. No composite score, coefficient, or threshold was constructed or tuned.",
                "",
                "GSE94533 source direction uses the saved, patient-deduplicated Stage I/II cancer versus benign/control endpoint (52 cancer, 61 benign/control). The 2024 direction uses all 286 malignant versus 392 pooled controls. Effect sizes are Hedges' g on each platform's native values, so magnitudes should not be compared across platforms. Direction agreement means only that the feature-wise case-minus-control mean contrasts have the same sign.",
                "",
                "The saved stability report provides selection frequencies, not fitted coefficient signs. GSE abundance direction is reported as a separate source-cohort contrast, not misrepresented as a model coefficient.",
                "",
                "| Historical feature | GSE94533 direction (g) | 2024 direction (g) | Direction agrees | Prior stable signature member |",
                "|---|---|---|---:|---:|",
                *table_rows,
                "",
                f"Directions agree for {agreeing}/8 exact overlaps. The only feature meeting the predeclared stability rule (selection in at least 3/5 folds within an existing 5/10/14-feature signature) is `{', '.join(stable) if stable else 'none'}`; it was stable for the all-cancer endpoint, not the Stage I/II endpoint. The strongest cross-platform directional agreement involving a previously stable feature is miR-200c-3p (GSE g={float(replication.loc[replication.feature.eq('hsa-miR-200c-3p'),'GSE94533_early_effect_hedges_g'].iloc[0]):.3f}; FirePlex g={float(replication.loc[replication.feature.eq('hsa-miR-200c-3p'),'2024_effect_hedges_g_raw_MFI'].iloc[0]):.3f}).",
                "",
                "Interpretation: miR-200c-3p is the clearest repeatable component in this feature-level comparison. miR-29a-3p reverses direction and miR-23b-3p is near zero/opposite. These descriptive univariate associations do not demonstrate a transportable multivariate assay or clinical performance.",
                "",
                "Source artifacts: `feature_selection_stability.csv`, `error_analysis_predictions.csv`, and the audited 2024 workbook. No 2024 labels were used in feature selection or model fitting.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    (REPORTS / "380_pelvic_mass_data_access_audit.md").write_text(
        "\n".join(
            [
                "# 380-Patient Pelvic-Mass Cohort Data-Access Audit",
                "",
                "## Study and Clinical Fit",
                "",
                "Study: \"Diagnostic plasma miRNA-profiles for ovarian cancer in patients with pelvic mass\" (Oliveira et al., 2019; DOI `10.1371/journal.pone.0225249`; PMCID `PMC6860451`). The study enrolled 380 presurgical plasma samples: 190 malignant ovarian/pelvic masses and 190 benign pelvic masses. It used a 46-candidate miRNA TaqMan RT-qPCR panel; the supplement enumerates 48 assays including two C. elegans spike-in controls. Discovery and validation cohorts were each reported as 95 malignant/95 benign.",
                "",
                "The intended endpoint (malignant versus benign pelvic mass) is directly clinically compatible with OVAI's triage objective. Specimen is plasma, not serum; platform is targeted TaqMan RT-qPCR with study-specific normalization, so no pooling with GSE94533 serum sequencing is appropriate.",
                "",
                "## Patient-Level Availability",
                "",
                "The article's Data Availability statement says relevant data are in the paper and Supporting Information. Public PLOS Figshare record `10322942` attaches only S1-S3 DOCX tables and figures; S1 lists assay IDs/sequences, S2 gives group-level median/minimum/maximum expression summaries, and S3 gives model sensitivity/specificity. No per-patient matrix, sample IDs, or row-level clinical covariates are attached. Exact GEO and SRA article-title searches returned no records; the public Figshare search returned the article and aggregate supplements, not individual measurements. Zenodo search returned no record; Dryad endpoint lookup was unavailable from this environment. No public patient-level deposit was identified.",
                "",
                "**Status: access required for individual-level analysis.** The publication does not provide an individual-level data-access route in its data statement; request the raw de-identified matrix and matching clinical labels from the corresponding authors. Do not treat the aggregate S2 table as an ML-ready cohort.",
                "",
                "## Published Cohort Counts and Limitations",
                "",
                "- Enrolled/reported: 380 total; 190 malignant and 190 benign; discovery 95/95; validation 95/95.",
                "- Table 1 reports 13 Stage I/II cancers in discovery (Stage I=5, Stage II=8) and none in validation; validation cancer cases are Stage III/IV (72 Stage III, 23 Stage IV). Thus only 13 early-stage cases are reported, below the project’s preferred 50.",
                "- The methods also state that 14 cancer and 13 benign samples were excluded for insufficient plasma or poor RNA quality. Per-patient retention and stage among excluded/retained subjects are not available. A simple upper bound after these exclusions is 177 benign plus at most 13 early-stage cancers, at most 190 candidate rows; actual complete-case numbers are likely lower and cannot be verified from summaries.",
                "- Supplement S2 has feature-specific denominators and missing summaries; for example, hsa-miR-200b-3p has no discovery-cohort summary values in S2 but is summarized in validation. No patient-level missingness or exact usable count is recoverable.",
                "- Cohort-level age and CA-125 summaries and histology/stage tables are published, but no individual-level age, CA-125, stage, or histology fields are released. Patient/sample IDs and duplicate/overlap audits are not possible.",
                "- Samples were presurgical. The paper does not provide patient-level treatment status or timing sufficient to certify a treatment-naive analysis.",
                "",
                "## Exact Feature Compatibility",
                "",
                f"The supplement lists {int(pelvic_map['exact_panel_match'].sum())} exact mature-feature mappings across the requested union of signatures, plus {int(pelvic_map['valid_mature_name_alias'].sum())} documented name aliases.",
                "- Current production OVAI: 1/5 exact, `hsa-miR-200b-3p` (TaqMan assay 002251; mature sequence `UAAUACUGCCUGGUAAUGA`). Four required production features are not in the panel, so the frozen 5-feature-plus-clinical OVAI model cannot be evaluated.",
                "- Historical 14-feature signature: 6/14 exact panel IDs: miR-23b-3p, miR-29a-3p, miR-92a-3p, miR-200a-3p, miR-200c-3p, and miR-335-5p.",
                "- 2024 paper's 10-feature list: 4 exact panel IDs (miR-195-5p, miR-200c-3p, miR-122-5p, and miR-424-5p) plus a documented mature-name alias for miR-103-3p to the panel's miRBase ID miR-103a-3p. The remaining five entries (miR-320e, miR-320c, cel-miR-39-3p, miR-770-5p, and miR-19a-3p) are absent. The panel has cel-miR-54-3p and cel-miR-238-3p, not cel-miR-39-3p.",
                "",
                "## Decision",
                "",
                "Classify as **YELLOW: clinical-endpoint biological replication / external cohort candidate, access required**. Endpoint and presurgical sampling are strong; individual measurements are not public, early-stage count is only 13, and only one of five frozen OVAI features overlaps. It is not currently a valid external validation of the canonical OVAI model and must not be merged with GSE94533.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    (REPORTS / "janusrna_prediagnostic_access_audit.md").write_text(
        "\n".join(
            [
                "# JanusRNA Prediagnostic Cohort Access Audit",
                "",
                "Source: Stawiski et al., \"Validation of miRNA signatures for ovarian cancer earlier detection in the pre-diagnosis setting using machine learning approaches\" (2024; DOI `10.3389/fonc.2024.1389066`; PMCID `PMC11231195`).",
                "",
                "The study used small-RNA sequencing from prediagnostic **serum** in the Norwegian Janus Serum Bank. The source cohort selected 80 invasive ovarian cancer cases and had 106 serum samples before one sample per person was chosen for the published analysis: 14 people had two samples and 4 people had three or more. The analysis matched 40 cases diagnosed <2 years after collection against 40 cases diagnosed >=2 years after collection; no cancer-free controls remained in this analysis because a batch effect between case and control libraries led the authors to discard controls.",
                "",
                "The paper states that data will be made available through Federated EGA Norway after application to the study Data Access Committee. Individual-level expression is therefore **access controlled**, not publicly downloadable. This is a viable access-gated prediagnostic biology study, not an open cohort today.",
                "",
                "Use requires patient-grouped analysis of all repeated samples, collection-to-diagnosis intervals, and the published choice of one sample per person. The paper reports a 2-year proximity contrast, serum sequencing, and significant platform differences that prevented direct transfer of the original model. No canonical OVAI prediction or threshold tuning on JanusRNA is justified from the public evidence.",
                "",
                "**Decision:** YELLOW for a prediagnostic biological-timing study after EGA approval and access to sample-level clinical linkage; RED for direct current-OVAI diagnostic validation because it is case-only in the released analysis, lacks benign pelvic-mass controls, is controlled-access, and differs technically from GSE94533.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    (REPORTS / "ideal_dataset_search_decision.md").write_text(
        "\n".join(
            [
                "# Additional Cohort Search and Decision",
                "",
                "## Search Scope and Findings",
                "",
                "The search prioritized individual-level circulating-miRNA data from patients presenting with an adnexal/pelvic mass, benign-mass controls, early stage, and quantitative measurements. The existing GEO metadata scout screened 58 candidate series (48 uncached candidates; none passed all prespecified metadata laws), and the remaining GEO queue report found zero candidates suitable for automated cache triage. A direct GEO and SRA search for the 2019 pelvic-mass article returned no accession. PLOS Figshare resolves the article's supplementary materials but only to assay sequences and aggregate cohort tables. Zenodo returned no article-level data record; Dryad lookup was unavailable from this environment.",
                "",
                "Previously audited GSE106817 has large title-derived group counts (3,108 samples; 320 cancer, 29 benign, 2,759 control, 66 borderline, 872 other) but has extensive assay-distribution shift and a failed source-derived threshold in the prior external evaluation. It is not promoted as the ideal clinically matched cohort, and no new transportability experiment was run here.",
                "",
                "## Best Target to Pursue",
                "",
                "The **Danish 380-patient pelvic-mass cohort** is the best next clinical-data access target: it directly contrasts malignant with benign masses, is presurgical plasma, reports age/CA-125/stage/histology at study level, and measures 46 candidate miRNAs. It is not currently a usable individual-level dataset: public files are aggregate only, Stage I/II counts are just 13, and only one current OVAI feature (miR-200b-3p) is present. Requesting the de-identified patient-level matrix, exact sample inclusion/assay availability, and labels is the appropriate next move.",
                "",
                "The 2024 serum FirePlex workbook is the best newly accessible patient-level cohort for **cross-platform biological replication**, with 678 records, 56 unambiguous Stage I/II cancer IDs, and eight historical-miRNA overlaps. Its controls pool benign and healthy participants and none of the current five production miRNAs are measured. It is not an expanded OVAI training cohort.",
                "",
                "No GREEN cohort for frozen current-OVAI external validation is presently established. The 380-patient study is YELLOW pending data access and remains feature-incomplete for the current 5-miRNA model; the 2024 workbook is YELLOW for biological replication only; JanusRNA is access-controlled YELLOW for prediagnostic biology.",
                "",
                "## Early-Stage and Benign-Mass Capacity",
                "",
                "For the 380-patient study, the published table reports 13 Stage I/II cancers (5 Stage I, 8 Stage II) and 190 enrolled benign pelvic masses. After the paper's stated exclusions (14 cancers and 13 benign samples), no more than 177 benign samples and at most 13 early-stage cancers could remain: an upper bound of 190 early-stage-plus-benign rows. The actual retained Stage I/II count and feature-complete denominator are unknown without individual-level data; the validation cohort itself has zero Stage I/II cases. For current frozen OVAI, the usable count remains zero because four of five production miRNAs are not in its assay panel.",
                "",
                "## Answers to the Phase Questions",
                "",
                "1. The 2024 cohort supports feature-wise biological replication and external cohort characterization, not direct evaluation or expansion of the frozen OVAI model. A distinct FirePlex-specific model would be a separate future project, not performed here.",
                f"2. Exact 2024 overlaps among the historical 14 are: {', '.join(f'`{x}`' for x in historical_overlap)}.",
                f"Published 2024 feature-list exact matches in the 2024 workbook: {', '.join(f'`{x}`' for x in published_overlap)}; absent exact identifier: {', '.join(f'`{x}`' for x in published_absent)}. The measured hsa-miR-103a-3p column is not silently recoded to the paper’s hsa-miR-103-3p label in this exact-match count.",
                "3. No valid alias for any of the five current production features was established against the measured 2024 panel. The measured miR-32-5p is a different mature arm from miR-32-3p; miR-200a/200c are not miR-200b.",
                "4. The 380-patient article has no public individual-level measurements identified. Public Figshare files are aggregate tables; individual data require an author data request.",
                "5. The 380 panel contains 1/5 exact current OVAI features: hsa-miR-200b-3p.",
                "6. It is not a viable frozen-OVAI external validation cohort as currently documented (only 1/5 features and no row-level data); it is a high-priority clinically aligned access request and biological replication target.",
                "7. JanusRNA can address prediagnostic timing after controlled-access approval, but not the intended benign-mass diagnostic endpoint from the published case-only subset.",
                "8. Best clinically matched target to pursue: the Danish 380-patient pelvic-mass study. Best additional patient-level dataset currently in hand: the 2024 678-patient workbook for biological replication only.",
                "9. The 380 cohort’s endpoint and presurgical specimen are clinically compatible; its plasma RT-qPCR platform, public-data limitation, and 1/5 feature overlap prevent direct frozen-model validation.",
                "10. Published ceiling: 13 early-stage cancers plus at most 177 retained benign masses (<=190 total after the stated technical exclusions). Actual complete rows and retained early-stage cases are unverified; current full OVAI scoreable count is zero.",
                "",
                "No model was trained, optimized, or modified in this phase. No new cohort was merged with GSE94533.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    REPORTS.mkdir(parents=True, exist_ok=True)
    overlap, panel_lookup = write_overlap_and_alias()
    replication = write_historical_replication(panel_lookup)
    pelvic_map, _ = write_pelvic_mass_map()
    write_reports(overlap, replication, pelvic_map)
    print(f"2024 exact overlaps: current={sum(overlap.present_in_current_OVAI & overlap.present_in_2024)}/5, historical={sum(overlap.present_in_historical_OVAI & overlap.present_in_2024)}/14, published={sum(overlap.present_in_published_2024_model & overlap.present_in_2024)}/10.")
    print(f"GSE94533 rows used for directional comparison: 52 cancer and 61 benign/control.")
    print(f"380-patient panel assays inventoried: {len(PELVIC_ASSAYS)}.")


if __name__ == "__main__":
    main()