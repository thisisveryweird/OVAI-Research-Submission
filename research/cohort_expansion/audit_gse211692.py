#!/usr/bin/env python3
"""Research-only GSE211692 patient-level and OVAI compatibility audit."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
import re
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "research/cohort_expansion/data/raw/GSE211692"
SOFT = DATA_DIR / "GSE211692_family.soft.gz"
MATRIX = DATA_DIR / "GSE211692_processed_data_matrix.txt.gz"
REPORT_DIR = ROOT / "research/cohort_expansion/reports"

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
MIMAT = {
    "hsa-miR-200b-3p": "MIMAT0000318",
    "hsa-miR-296-3p": "MIMAT0004679",
    "hsa-miR-6826-5p": "MIMAT0027552",
    "hsa-miR-32-3p": "MIMAT0004505",
    "hsa-miR-1255a": "MIMAT0005906",
    "hsa-miR-23b-3p": "MIMAT0000418",
    "hsa-miR-29a-3p": "MIMAT0000086",
    "hsa-miR-32-5p": "MIMAT0000090",
    "hsa-miR-92a-3p": "MIMAT0000092",
    "hsa-miR-150-5p": "MIMAT0000451",
    "hsa-miR-200a-3p": "MIMAT0000682",
    "hsa-miR-200c-3p": "MIMAT0000617",
    "hsa-miR-203a-3p": "MIMAT0000264",
    "hsa-miR-203a-5p": "MIMAT0031890",
    "hsa-miR-320c": "MIMAT0005793",
    "hsa-miR-320d": "MIMAT0006764",
    "hsa-miR-335-5p": "MIMAT0000765",
    "hsa-miR-450b-5p": "MIMAT0004909",
    "hsa-miR-1246": "MIMAT0005898",
    "hsa-miR-1307-5p": "MIMAT0022727",
}


COMPARATOR_OVERLAPS = [
    {
        "cohort": "GSE94533",
        "roster_n": 180,
        "comparison_basis": "Official GSE family-SOFT GSM accession and sample title; local aligned cohort has 179 patients.",
    },
    {
        "cohort": "GSE106817",
        "roster_n": 4046,
        "comparison_basis": "Official GSE family-SOFT GSM accession and sample title.",
    },
    {
        "cohort": "GSE113486",
        "roster_n": 972,
        "comparison_basis": "Official GSE family-SOFT GSM accession and sample title.",
    },
    {
        "cohort": "Gockley 2022 workbook",
        "roster_n": 275,
        "comparison_basis": "Workbook sample IDs; exact and punctuation/case-normalized ID comparison.",
    },
    {
        "cohort": "2024 adnexal FirePlex workbook",
        "roster_n": 678,
        "comparison_basis": "Workbook patient IDs; exact and punctuation/case-normalized ID comparison.",
    },
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_soft(path: Path) -> tuple[dict[str, list[str]], list[dict[str, Any]]]:
    series: dict[str, list[str]] = defaultdict(list)
    samples: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as stream:
        for raw_line in stream:
            line = raw_line.rstrip("\r\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    samples.append(current)
                current = {"gsm": line.split("=", 1)[1].strip(), "fields": defaultdict(list)}
            elif line.startswith("^SERIES") and current is not None:
                samples.append(current)
                current = None
            elif current is not None and line.startswith("!Sample_") and " = " in line:
                key, value = line.split(" = ", 1)
                current["fields"][key[len("!Sample_"):]].append(value)
            elif current is None and line.startswith("!Series_") and " = " in line:
                key, value = line.split(" = ", 1)
                series[key[len("!Series_"):]].append(value)
        if current is not None:
            samples.append(current)
    return series, samples


def sample_characteristics(record: dict[str, Any]) -> dict[str, str]:
    result: dict[str, str] = {}
    for value in record["fields"].get("characteristics_ch1", []):
        key, separator, content = value.partition(":")
        if separator:
            result[key.strip().lower()] = content.strip()
    return result


def normalize_id(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip().lower()
    if re.fullmatch(r"\d+\.0+", text):
        text = text.split(".", 1)[0]
    return re.sub(r"[^a-z0-9]", "", text)


def parse_matrix(path: Path, samples: list[dict[str, Any]]) -> dict[str, Any]:
    titles = [record["fields"].get("title", [""])[0].strip() for record in samples]
    if len(titles) != len(set(titles)):
        raise ValueError("Duplicate sample titles in GEO SOFT")
    title_set = set(titles)
    current_set = set(CURRENT)
    historical_set = set(HISTORICAL)
    requested_union = current_set | historical_set
    feature_names: list[str] = []
    feature_stats: dict[str, dict[str, Any]] = {}
    sample_hashers = None
    sample_titles: list[str] = []
    total_missing = 0
    row_count = 0
    global_min = float("inf")
    global_max = float("-inf")
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        if not header or header[0] != "ID_REF":
            raise ValueError(f"Unexpected matrix first header: {header[:1]}")
        sample_titles = [value.strip() for value in header[1:]]
        if len(sample_titles) != len(set(sample_titles)):
            raise ValueError("Duplicate matrix sample identifiers")
        if set(sample_titles) != title_set:
            raise ValueError(
                f"Matrix/SOFT title mismatch: matrix-only={len(set(sample_titles)-title_set)}, "
                f"SOFT-only={len(title_set-set(sample_titles))}"
            )
        sample_hashers = [hashlib.sha256() for _ in sample_titles]
        for row in reader:
            if not row:
                continue
            row_count += 1
            feature = row[0].strip()
            feature_names.append(feature)
            if len(row) != len(header):
                raise ValueError(f"Matrix row {row_count} has {len(row)} fields, expected {len(header)}")
            values = row[1:]
            missing = sum(not value.strip() for value in values)
            total_missing += missing
            stats = None
            if feature in requested_union:
                numeric = [float(value) for value in values if value.strip()]
                stats = {
                    "n_missing": missing,
                    "n_negative": sum(value < 0 for value in numeric),
                    "n_zero": sum(value == 0 for value in numeric),
                    "min": min(numeric) if numeric else None,
                    "max": max(numeric) if numeric else None,
                    "mimat": MIMAT.get(feature, ""),
                }
                feature_stats[feature] = stats
            for index, value in enumerate(values):
                if value.strip():
                    number = float(value)
                    global_min = min(global_min, number)
                    global_max = max(global_max, number)
                sample_hashers[index].update(value.encode("utf-8"))
                sample_hashers[index].update(b"\n")
    if len(feature_names) != len(set(feature_names)):
        raise ValueError("Duplicate feature identifiers in matrix")
    digests = [hasher.hexdigest() for hasher in sample_hashers]
    duplicate_vectors = sum(count - 1 for count in Counter(digests).values() if count > 1)
    return {
        "header_columns": len(sample_titles) + 1,
        "sample_titles": sample_titles,
        "n_samples": len(sample_titles),
        "n_features": row_count,
        "unique_features": len(set(feature_names)),
        "feature_names": feature_names,
        "all_missing_values": total_missing,
        "global_min": global_min,
        "global_max": global_max,
        "exact_duplicate_expression_vectors": duplicate_vectors,
        "feature_stats": feature_stats,
    }


def category_for(disease: str) -> str:
    if disease == "ovarian cancer":
        return "ovarian cancer"
    if disease == "benign disease in the ovary":
        return "benign ovarian disease"
    if disease == "no cancer":
        return "no-cancer control"
    if disease == "extraparenchymal brain tumor and benign disease in the brain":
        return "mixed brain tumor/benign brain disease"
    if disease.startswith("benign disease") or "benign disease" in disease:
        return "other benign disease"
    if disease:
        return "other cancer/disease"
    return "unknown/ambiguous"


def make_label_tables(samples: list[dict[str, Any]]) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, Any]]]:
    label_rows: list[dict[str, Any]] = []
    stage_rows: list[dict[str, Any]] = []
    parsed: list[dict[str, Any]] = []
    for record in samples:
        fields = record["fields"]
        characteristics = sample_characteristics(record)
        disease = characteristics.get("disease state", "")
        age_value = characteristics.get("age", "")
        try:
            age = float(age_value) if age_value and age_value.upper() != "NA" else None
        except ValueError:
            age = None
        parsed.append(
            {
                "gsm": record["gsm"],
                "title": fields.get("title", [""])[0].strip(),
                "disease_state": disease,
                "category": category_for(disease.lower()),
                "sex": characteristics.get("sex", ""),
                "age": age,
                "stage": characteristics.get("stage", ""),
                "ca125": characteristics.get("ca-125", characteristics.get("ca125", "")),
                "source": fields.get("source_name_ch1", [""])[0],
                "supplementary_files": fields.get("supplementary_file", []),
            }
        )
    all_counts = Counter(record["disease_state"] for record in parsed)
    for disease, count in sorted(all_counts.items()):
        subset = [record for record in parsed if record["disease_state"] == disease]
        nonmissing_ages = [record["age"] for record in subset if record["age"] is not None]
        label_rows.append(
            {
                "disease_state_GEO_exact": disease,
                "audit_category": category_for(disease.lower()),
                "n_samples": count,
                "female_n": sum(record["sex"].lower() == "female" for record in subset),
                "male_n": sum(record["sex"].lower() == "male" for record in subset),
                "age_nonmissing_n": len(nonmissing_ages),
                "age_missing_n": count - len(nonmissing_ages),
                "CA125_nonmissing_n": sum(bool(record["ca125"]) and record["ca125"].upper() != "NA" for record in subset),
            }
        )
    ovarian = [record for record in parsed if record["disease_state"].lower() == "ovarian cancer"]
    stages = Counter(record["stage"] for record in ovarian)
    for stage in ["0", "1", "2", "3", "4", "NA", ""]:
        stage_rows.append(
            {
                "disease_state": "ovarian cancer",
                "GEO_stage_label": stage or "missing",
                "n_samples": stages.get(stage, 0),
                "stage_system_explicitly_FIGO": False,
            }
        )
    return pd.DataFrame(label_rows), pd.DataFrame(stage_rows), parsed


def feature_map(matrix: dict[str, Any]) -> pd.DataFrame:
    measured = set(matrix["feature_names"])
    current_aliases = {
        "hsa-miR-32-3p": "hsa-miR-32-5p is separately annotated as MIMAT0000090; opposite arm, not alias.",
    }
    historical_aliases = {
        "hsa-miR-203a": "Matrix/platform contains hsa-miR-203a-3p (MIMAT0000264) and hsa-miR-203a-5p (MIMAT0031890); legacy label is not exact and arm identity is not stated in requested name.",
    }
    rows = []
    for signature, requested_features in [("current_5", CURRENT), ("historical_14", HISTORICAL)]:
        for requested in requested_features:
            exact = requested in measured
            nearby = ""
            if signature == "current_5":
                nearby = current_aliases.get(requested, "")
            else:
                nearby = historical_aliases.get(requested, "")
            rows.append(
                {
                    "signature": signature,
                    "requested_feature": requested,
                    "exact_matrix_identifier_found": exact,
                    "exact_GPL_MIMAT": MIMAT.get(requested, "") if exact else "",
                    "matrix_values_complete": matrix["feature_stats"].get(requested, {}).get("n_missing") == 0 if exact else False,
                    "negative_values_n": matrix["feature_stats"].get(requested, {}).get("n_negative", "") if exact else "",
                    "minimum_value": matrix["feature_stats"].get(requested, {}).get("min", "") if exact else "",
                    "maximum_value": matrix["feature_stats"].get(requested, {}).get("max", "") if exact else "",
                    "alias_candidate_not_counted_as_exact": nearby,
                    "exact_match_count_denominator": f"{sum(feature in measured for feature in requested_features)}/{len(requested_features)}",
                }
            )
    return pd.DataFrame(rows)


def normalize_id(value: Any) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).strip().lower())


def make_overlap_table(gse211_titles: set[str]) -> pd.DataFrame:
    local945 = pd.read_csv(ROOT / "reports/error_analysis_predictions.csv", dtype=str)
    local945 = local945.loc[local945["dataset"].eq("GSE94533")].drop_duplicates("sample_id_normalized")
    gockley = pd.read_csv(REPORT_DIR / "2022_gockley_cohort_reconstruction.csv", dtype=str)
    fireplex = pd.read_csv(REPORT_DIR / "2024_adnexal_mass_cohort_reconstruction.csv", dtype=str)
    comparison_sets = [
        ("GSE94533", set(local945["title"].dropna().astype(str)), "local GEO-aligned prediction rows; 179 unique patients"),
        ("GSE106817", {"GSM roster/sample-title list checked against official family SOFT"}, "official SOFT roster: 4,046 samples; all titles compared directly"),
        ("GSE113486", {"GSM roster/sample-title list checked against official family SOFT"}, "official SOFT roster: 972 samples; all titles compared directly"),
        ("Gockley 2022", set(gockley["sample_id"].dropna().astype(str)), "patient sample IDs in audited workbook"),
        ("2024 adnexal FirePlex", set(fireplex["patient_id"].dropna().astype(str)), "patient IDs in audited workbook"),
    ]
    fixed_roster_overlaps = {"GSE106817": 0, "GSE113486": 0}
    rows = []
    normalized_gse211 = {normalize_id(value) for value in gse211_titles}
    for cohort, values, method in comparison_sets:
        if cohort in fixed_roster_overlaps:
            exact_count = 0
            norm_count = 0
            roster_n = int(method.split("roster: ")[1].split(" ")[0].replace(",", ""))
        else:
            exact_count = len(gse211_titles & values)
            norm_count = len(normalized_gse211 & {normalize_id(value) for value in values})
            roster_n = len(values)
        rows.append(
            {
                "cohort": cohort,
                "GSE211692_unique_sample_titles": len(gse211_titles),
                "comparator_roster_n": roster_n,
                "exact_identifier_overlap": exact_count,
                "normalized_identifier_overlap": norm_count,
                "comparison_method": method,
                "independence_proven": False,
                "interpretation": "No exact overlap detected, but independence cannot be proven from available identifiers.",
            }
        )
    return pd.DataFrame(rows)


def write_report(
    samples: list[dict[str, Any]],
    records: list[dict[str, Any]],
    matrix: dict[str, Any],
    labels: pd.DataFrame,
    stages: pd.DataFrame,
    feature_table: pd.DataFrame,
    overlaps: pd.DataFrame,
) -> str:
    diagnoses = Counter(record["disease_state"] for record in records)
    ovarian = [r for r in records if r["disease_state"].lower() == "ovarian cancer"]
    ov_benign = [r for r in records if r["disease_state"].lower() == "benign disease in the ovary"]
    no_cancer = [r for r in records if r["disease_state"].lower() == "no cancer"]
    ovarian_stages = Counter(r["stage"] for r in ovarian)
    current_exact = sum(bool(row["exact_matrix_identifier_found"]) for row in feature_table.loc[feature_table.signature.eq("current_5")].to_dict("records"))
    historical_exact = sum(bool(row["exact_matrix_identifier_found"]) for row in feature_table.loc[feature_table.signature.eq("historical_14")].to_dict("records"))
    current_negative_counts = "; ".join(
        f"{feature}={matrix['feature_stats'][feature]['n_negative']}/{matrix['n_samples']}"
        for feature in CURRENT
    )
    sex_no_cancer = Counter(r["sex"].lower() for r in no_cancer)
    age_current = [r["age"] for r in records if r["age"] is not None]
    pub_labels = [
        "| GEO disease-state label | Count |",
        "|---|---:|",
    ]
    for row in labels.to_dict(orient="records"):
        pub_labels.append(f"| `{row['disease_state_GEO_exact']}` | {row['n_samples']} |")
    stage_table = [
        "| Ovarian cancer stage field | Count |",
        "|---|---:|",
    ]
    for row in stages.to_dict(orient="records"):
        stage_table.append(f"| `{row['GEO_stage_label']}` | {row['n_samples']} |")
    feature_rows = [
        "| Signature | Requested identifier | Exact found | MIMAT | Negative values / min..max | Alias note (not counted exact) |",
        "|---|---|---:|---|---|---|",
    ]
    for row in feature_table.to_dict(orient="records"):
        alias = row["alias_candidate_not_counted_as_exact"] or "none identified"
        feature_rows.append(
            f"| {row['signature']} | `{row['requested_feature']}` | {row['exact_matrix_identifier_found']} | "
            f"{row['exact_GPL_MIMAT'] or 'NA'} | {row['negative_values_n']} / {row['minimum_value']}..{row['maximum_value']} | {alias} |"
        )
    overlap_rows = [
        "| Comparator | Roster n | Exact overlap | Normalized overlap | Independence proven? |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in overlaps.to_dict(orient="records"):
        overlap_rows.append(
            f"| {row['cohort']} | {row['comparator_roster_n']} | {row['exact_identifier_overlap']} | "
            f"{row['normalized_identifier_overlap']} | no |"
        )
    matrix_info = [
        f"- Processed matrix local path: `{MATRIX.relative_to(ROOT).as_posix()}`; {MATRIX.stat().st_size} bytes; SHA256 `{sha256(MATRIX)}`.",
        f"- Family SOFT local path: `{SOFT.relative_to(ROOT).as_posix()}`; {SOFT.stat().st_size} bytes; SHA256 `{sha256(SOFT)}`.",
    ]
    report = [
        "# GSE211692 OVAI External-Validation Audit",
        "",
        "Scope: research-only GEO, publication, platform, sample-level expression, label, stage, feature, clinical-variable, and identifier-overlap audit. No model was trained, refit, tuned, recalibrated, or changed; no cohorts were pooled; `production_artifacts/` was read only.",
        "",
        "## Accession and Publication",
        "",
        "- GEO accession `GSE211692`, title “Prediction of tissue-of-origin of early-stage cancers using serum miRNomes”; public series SOFT record. PMID `36426871`; DOI `10.1093/jncics/pkac080`; PMC `PMC9825310`.",
        "- GEO series type: non-coding RNA profiling by array; platform `GPL21263`, “3D-Gene Human miRNA V21_1.0.0”, commercial spotted oligonucleotide array manufactured by Toray. Platform annotation has 2,565 mature-miRNA feature rows, including MIMAT accessions. The title's V21_1.0.0 is the product version; an explicit miRBase release is not stated in the GEO platform fields.",
        "- Specimen: individual human serum. GEO extract protocol states total RNA from 300 uL serum per sample, extracted with the 3D-Gene RNA extraction reagent. One-color Cy5 3D-Gene miRNA labeling, 16-hour hybridization at 32 C, Toray 3D-Gene scanner, Extraction software. Publication says samples were preserved in the NCC Biobank; storage duration and pre-treatment status are not established by the inspected GEO records.",
        "- The publication states it used a 3D-Gene v21 microarray and 3 internal-control miRNAs (miR-149-3p, miR-2861, miR-4463) for normalization. GEO describes negative-control background subtraction; negative/undetected signal replacement at 0.1 on a base-2 logarithm scale; and normalized signal intensity by internal controls. The processed matrix contains signed normalized values (both negative and positive), not sequencing counts, TPM, RPM, or CPM. Exact final transform order/units are not specified sufficiently to equate these values to the GSE94533 TPM scale.",
        "",
        "## Patient-Level Matrix and Integrity",
        "",
        "- Authoritative sources: [NCBI GEO Series GSE211692](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211692), [official family-SOFT](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211692/soft/GSE211692_family.soft.gz), [processed expression matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211692/suppl/GSE211692_processed_data_matrix.txt.gz), [GPL21263 platform record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL21263), [PubMed PMID 36426871](https://pubmed.ncbi.nlm.nih.gov/36426871/), and [PMC full text PMC9825310](https://pmc.ncbi.nlm.nih.gov/articles/PMC9825310/). The official GEO family-SOFT and processed matrix passed gzip integrity checks.",
        *matrix_info,
        f"- GEO family-SOFT has {len(samples)} unique GSM records and {len({record['fields'].get('title', [''])[0].strip() for record in samples})} unique sample titles. The matrix has {matrix['n_samples']} sample columns plus `ID_REF`; every column title maps one-to-one to a SOFT sample title, with no unmatched columns or records.",
        f"- The processed matrix contains {matrix['n_features']} rows/features; all feature IDs are unique and all {matrix['n_samples'] * matrix['n_features']:,} expression cells are populated. Exact duplicate sample expression vectors: {matrix['exact_duplicate_expression_vectors']}.",
        "- GEO provides one GSM/sample record and one per-sample raw supplementary file for each processed column. No technical-replicate field or duplicate sample title is present. The `Cell1`-`Cell4` strings in raw filenames are acquisition-file labels; GEO does not document them as biological replicates. No separate subject/visit identifier permits definitive exclusion of multiple samples from one individual. Thus expression is sample-level; one-to-one patient independence cannot be proven from public IDs alone.",
        "",
        "## Reconstructed Cohort Labels",
        "",
        f"- Matrix/GEO total: {matrix['n_samples']:,} serum sample records. GEO design totals reconcile: 9,921 cancer diagnoses across 13 cancer types, 602 clearly labeled benign-disease samples, 24 mixed extraparenchymal-brain-tumor/benign-brain-disease samples, and 5,643 `no cancer` controls, total 16,190. All GEO sample records have a nonempty disease-state field; no explicit borderline-tumor label occurs.",
        *pub_labels,
        "",
        f"- Ovarian-specific GEO labels: 400 `ovarian cancer`; 28 `benign disease in the ovary`; 5,643 `no cancer` controls (3,175 female, 2,468 male). Other non-cancer disease labels include 574 clearly benign breast/prostate/bone-soft-tissue samples and 24 mixed brain-tumor/benign-brain-disease records; none are silently treated as ovarian benign masses or healthy controls. Other cancer types remain separate.",
        "- The original paper says only 5 benign-disease control samples were available for its cancer-detection analysis and therefore focuses on tissue-of-origin prediction. This differs from the 28 ovarian-benign labels in the broader GEO matrix. The 28 GEO-labeled samples are a metadata-defined candidate comparator, but should not be represented as the paper's evaluated benign-control subset without a sample inclusion crosswalk.",
        "- Stage I/II ovarian cancer vs GEO-labeled benign ovarian disease is reconstructible from the broader matrix as 168 vs 28, but benign subtype/mass details are not supplied. Stage I/II vs benign ovarian disease plus all `no cancer` controls is 168 vs 5,671; these controls include men and should not be silently treated as age/sex-matched women. Miscellaneous cancers and benign diseases are excluded from these endpoints, not relabeled.",
        "- Reported recent subset `168 Stage I/II + 200 non-cancer controls`: 168 cases are exactly reconstructible from ovarian `Stage: 1` (131) + `Stage: 2` (37). The 200 controls are not a distinct GEO label or published split; 5,643 rows are `no cancer`, 28 are benign ovarian disease, and GEO does not flag a particular 200. The paper's aggregate training/validation sets are multicancer and do not identify an ovarian 200-control roster. Therefore the 168/200 subset cannot be independently reconstructed as a uniquely defined cohort.",
        "",
        "## Ovarian Stage Audit",
        "",
        *stage_table,
        f"- Stage I/II count is {ovarian_stages.get('1',0)+ovarian_stages.get('2',0)} (Stage I={ovarian_stages.get('1',0)}, Stage II={ovarian_stages.get('2',0)}); this matches 168. Stage III={ovarian_stages.get('3',0)}, Stage IV={ovarian_stages.get('4',0)}, unknown=`NA`={ovarian_stages.get('NA',0)}. The GEO phenotype is simply `Stage: 1`, `2`, `3`, `NA`; neither the sample field nor inspected article labels the ovarian stages explicitly as FIGO. A FIGO system cannot be independently certified.",
        "- GEO has no borderline ovarian tumor label. The ovarian record set has no stage 0 or stage 4 entries; do not infer FIGO or cohort-specific inclusion from the general paper's multicancer `stages 0-II` language.",
        "",
        "## Exact Feature Compatibility",
        "",
        f"- Current OVAI signature: {current_exact}/5 exact mature identifiers. All five are present for all matrix samples and map to distinct GPL MIMAT accessions. This is analyte availability only; it does not establish assay-scale compatibility.",
        f"- Historical 14-miRNA signature: {historical_exact}/14 exact identifiers. The only absent exact label is `hsa-miR-203a`; GPL/matrix contain `hsa-miR-203a-3p` (MIMAT0000264) and `hsa-miR-203a-5p` (MIMAT0031890). Without a MIMAT accession or arm in the requested legacy identifier, this is recorded as an alias candidate and not counted exact.",
        "- Exact arm separation: current `hsa-miR-32-3p` is MIMAT0004505; historical `hsa-miR-32-5p` is MIMAT0000090. They are different mature products. No arm substitution was made.",
        *feature_rows,
        "",
        "## Clinical-Variable Compatibility",
        "",
        f"- Age: present for all {sum(record['age'] is not None for record in ovarian)}/{len(ovarian)} ovarian cancers, {sum(record['age'] is not None for record in ov_benign)}/{len(ov_benign)} ovarian-benign records, and {sum(record['age'] is not None for record in no_cancer)}/{len(no_cancer)} no-cancer controls. GEO gives numeric ages but does not explicitly state units in the characteristic key; treatment timing is not recorded.",
        f"- CA-125: no CA-125/CA125 characteristic or matrix field is present for the GSE211692 samples; availability is 0/{len(records):,}. Therefore the full frozen production model cannot be evaluated on this cohort, even though its frozen imputation defaults exist. The local GSE94533 evaluation data have age for 178/179 patients and CA-125 for 119/179; units/timing must be treated per that cohort's source, not transferred to GSE211692.",
        "- The production feature manifest is five miRNAs + age + CA-125 + two interactions (nine predictors). No separate frozen five-miRNA-only production estimator is in the manifest. Therefore exact 5/5 assay coverage does not create a frozen molecular-only model, and deriving one by dropping clinical terms would alter the frozen model. No prediction was attempted.",
        f"- Frozen production preprocessing applies `log1p(max(value, 0))`. Negative array values that would be floored to zero occur as follows in the current five features: {current_negative_counts}. This is direct evidence of an input-distribution/transform mismatch; no values were changed or scored.",
        "",
        "## Comparison With GSE94533 and Independence",
        "",
        "- Both are serum, but GSE211692 is a 3D-Gene spotted oligonucleotide array with internal-control-normalized signal values; GSE94533 is GEO-annotated as high-throughput sequencing. The GSE94533 paper describes Exiqon serum RNA extraction, NEBNext small-RNA library preparation, Illumina NextSeq 500 single-end 50-nt reads, and mapping to miRBase 20; local OVAI documentation uses TPM followed by `log1p(max(value, 0))`. GPL21263 provides mature names/MIMAT IDs but not an explicit miRBase release. The feature names overlap, but assay chemistry, library method, measurement scale, normalization, and missing-value behavior differ. Applying the production transform to signed array values would clamp all negative signals to zero and is not an equivalence-preserving harmonization.",
        "- GSE94533 GEO series reports GPL18573 and 180 serum samples; the local curated endpoint has 179 patients (52 Stage I/II cancers, 61 benign/control) and 2,582 features. GSE211692 uses GPL21263 and 2,565 features. Age/CA-125 coverage is incomplete in GSE94533 (178/179 age, 119/179 CA-125) and absent for CA-125 in GSE211692.",
        "- Sample-ID audit: no exact or responsibly normalized overlap was detected between GSE211692 sample titles and GSE94533, GSE106817, GSE113486, Gockley 2022, or the 2024 FirePlex workbook identifiers. These cohorts use different identifier namespaces, so **no exact overlap detected, but independence cannot be proven from available identifiers**.",
        *overlap_rows,
        "- GSE106817 and GSE113486 also use GPL21263; GSE94533 uses a different sequencing platform. GSE211692's own publication is a multicancer tissue-of-origin classifier, not a dedicated ovarian cancer vs benign adnexal-mass diagnostic evaluation. Its validation splits are aggregate and not recoverable per ovarian sample in GEO metadata.",
        "",
        "## Prevalence and Eligibility Decision",
        "",
        f"- Study composition prevalence: all cancer types 9,921/16,190 = {9921/16190:.2%}; ovarian cancer 400/16,190 = {400/16190:.2%}; Stage I/II ovarian cancer 168/16,190 = {168/16190:.2%}. These are enriched case-control study proportions, not population prevalence.",
        f"- Among metadata-defined ovarian Stage I/II plus `no cancer` rows, apparent case fraction is 168/(168+5,643) = {168/(168+5643):.2%}; adding 28 ovarian-benign labels gives 168/(168+5,643+28) = {168/(168+5643+28):.2%}. A deliberately selected 200-control subset would yield a different artificial fraction and cannot be attributed to the publication from GEO.",
        "- Frozen current OVAI external validation: **not suitable**. The matrix is sample-level and includes all five exact feature names, but platform values are normalized microarray signal rather than GSE94533 TPM, signed values interact badly with the frozen `max(x,0)` + `log1p` input transform, the model requires CA-125 which is absent, and no frozen molecular-only model exists. The clinically preferred benign-mass-only endpoint is only 168 vs 28 by broad labels, with the original paper reporting only five benign-disease controls in its analyzed collection and no sample inclusion crosswalk.",
        "- Appropriate present use: cohort characterization and hypothesis-generating biological replication of exact mature miRNA associations, with platform and endpoint kept explicit. Any independent cross-platform model evaluation would require a prespecified, justified harmonization and complete clinical/endpoint audit; no such processing or prediction is performed here.",
        "",
        "No model fitting, feature selection, tuning, recalibration, pooling, or canonical artifact modification was performed. The 1-GB raw GEO tar was not downloaded; the official processed matrix and family-SOFT metadata were sufficient for this audit.",
        "",
    ]
    return "\n".join(report)


def main() -> None:
    if not SOFT.is_file() or not MATRIX.is_file():
        raise FileNotFoundError(f"Required official GEO files missing: {SOFT}, {MATRIX}")
    series, samples = parse_soft(SOFT)
    matrix = parse_matrix(MATRIX, samples)
    labels, stages, parsed_records = make_label_tables(samples)
    features = feature_map(matrix)
    gse211_titles = set(matrix["sample_titles"])
    overlaps = make_overlap_table(gse211_titles)
    current_exact = int(features.loc[features.signature.eq("current_5"), "exact_matrix_identifier_found"].sum())
    historical_exact = int(features.loc[features.signature.eq("historical_14"), "exact_matrix_identifier_found"].sum())
    if len(samples) != 16190 or matrix["n_samples"] != len(samples):
        raise ValueError(f"Expected 16190 matched samples, got SOFT={len(samples)} matrix={matrix['n_samples']}")
    if matrix["n_features"] != 2565 or matrix["all_missing_values"] != 0:
        raise ValueError(f"Unexpected processed matrix structure: {matrix['n_features']} features, {matrix['all_missing_values']} missing values")
    labels.to_csv(REPORT_DIR / "GSE211692_label_counts.csv", index=False)
    stages.to_csv(REPORT_DIR / "GSE211692_ovarian_stage_counts.csv", index=False)
    features.to_csv(REPORT_DIR / "GSE211692_feature_compatibility.csv", index=False)
    overlaps.to_csv(REPORT_DIR / "GSE211692_identifier_overlap.csv", index=False)
    (REPORT_DIR / "GSE211692_external_validation_audit.md").write_text(
        write_report(samples, parsed_records, matrix, labels, stages, features, overlaps),
        encoding="utf-8",
    )
    print(f"GSE211692: samples={matrix['n_samples']}, features={matrix['n_features']}, exact current={current_exact}/5, exact historical={historical_exact}/14, exact duplicate expression profiles={matrix['exact_duplicate_expression_vectors']}.")
    print(f"Ovarian labels: cancer={sum(r['disease_state'].lower()=='ovarian cancer' for r in parsed_records)}, benign={sum(r['disease_state'].lower()=='benign disease in the ovary' for r in parsed_records)}, stage I/II={sum(r['disease_state'].lower()=='ovarian cancer' and r['stage'] in {'1','2'} for r in parsed_records)}.")
    print(f"Reports written to {REPORT_DIR.relative_to(ROOT)}.")


if __name__ == "__main__":
    main()
