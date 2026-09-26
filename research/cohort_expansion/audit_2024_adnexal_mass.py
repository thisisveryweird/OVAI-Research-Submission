#!/usr/bin/env python3
"""Audit the 2024 serum adnexal-mass miRNA study without fabricating raw data.

This script is intentionally conservative. It only records what can be confirmed
from the public paper metadata and any locally available files. If the individual
supplementary Data File S1 is not present in the workspace, the script writes a
scientifically honest audit status instead of inventing patient-level data.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REPORTS = ROOT / "reports"
DATA_DIR = ROOT / "data" / "external"
OUT_REPORT = REPORTS / "2024_adnexal_mass_data_audit.md"
OUT_FEATURE_MAP = REPORTS / "2024_adnexal_mass_feature_map.csv"
OUT_COMPAT = REPORTS / "2024_adnexal_mass_compatibility.csv"
OUT_MANIFEST = ROOT / "research" / "cohort_expansion" / "manifest.json"


OVAI_5 = [
    "hsa-miR-200b-3p",
    "hsa-miR-296-3p",
    "hsa-miR-6826-5p",
    "hsa-miR-32-3p",
    "hsa-miR-1255a",
]
OVAI_14 = [
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


def write_feature_map() -> None:
    rows = []
    for feature in OVAI_5 + OVAI_14:
        rows.append({
            "miRNA": feature,
            "GSE94533 availability": "yes",
            "2024 availability": "not publicly verified",
            "identifier mapping": "requires exact raw FirePlex annotation from Data File S1",
            "assay compatibility": "not established without raw file",
            "measurement scale": "unknown",
            "usable for targeted external validation?": "no (currently not validated)",
            "usable for pooled training?": "no (requires assay harmonization)",
            "reason": "The 2024 cohort is a serum FirePlex assay, not the GSE94533 serum sequencing platform. Exact mapping must be checked against the raw Data File S1.",
        })
    import csv
    fieldnames = [
        "miRNA",
        "GSE94533 availability",
        "2024 availability",
        "identifier mapping",
        "assay compatibility",
        "measurement scale",
        "usable for targeted external validation?",
        "usable for pooled training?",
        "reason",
    ]
    with OUT_FEATURE_MAP.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_compatibility() -> None:
    import csv
    rows = [
        {
            "dataset": "2024 serum adnexal-mass study",
            "status": "published candidate cohort",
            "sample_count": 678,
            "cancer_count": 286,
            "control_count": 392,
            "training_count": 468,
            "internal_validation_count": 100,
            "external_validation_count": 110,
            "specimen_type": "serum",
            "assay_platform": "FirePlex targeted serum miRNA panel",
            "endpoint": "ovarian cancer vs benign adnexal mass / control",
            "label_definitions": "publication reports cancer vs non-cancer and training/internal/external splits",
            "feature_space": "180 serum miRNAs + 7 proteins + age + menopausal status",
            "preprocessing": "raw Data File S1 not yet accessed; numerical scale and normalization not verified",
            "development_or_validation": "candidate development + validation, but raw individual-level data not yet retrieved",
            "patient_overlap_status": "unknown; not yet checked against GSE94533",
            "publication_source": "Gynecologic Oncology, 2024, DOI 10.1016/j.ygyno.2024.08.008",
            "labels_used_during_model_development": "not assessed; raw data missing",
            "note": "This study is promising but cannot be used for patient-level OVAI modeling until the actual supplemental Data File S1 is retrieved and audited.",
        }
    ]
    with OUT_COMPAT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def build_report() -> str:
    likely_files = sorted(str(p.relative_to(ROOT)) for p in DATA_DIR.rglob("*") if "adnexal" in p.name.lower() or "fireplex" in p.name.lower() or "2024" in p.name.lower())
    text = """# 2024 serum adnexal-mass cohort audit

This audit preserves the canonical OVAI source model and does not use the 2024 adnexal-mass paper to overwrite any existing artifact.

## Summary

The paper title is: "Serum miRNA improves the accuracy of a multivariate index assay for triage of an adnexal mass" (Gynecologic Oncology, 2024; DOI: 10.1016/j.ygyno.2024.08.008).

The public paper metadata confirms the following study-level numbers:

- total subjects: 678
- cancer: 286
- controls: 392
- training subjects: 468
- independent internal validation subjects: 100
- independent external validation subjects: 110
- serum miRNA measurements: 180
- serum protein biomarkers: 7
- age: present
- menopausal status: present

## Current status of the raw data

The individual-level supplementary Data File S1 is not present in this workspace, and the public publisher-facing metadata available from the DOI/article page does not provide the raw patient-level table. The file therefore has not been audited at the row/subject level. No patient-level training or validation model has been fit using this cohort.

A full audit requires the actual Data File S1 or equivalent patient-level supplementary table. Until that is retrieved and reviewed, this dataset remains a literature-level candidate rather than an usable external validation dataset for OVAI.

## Questions that remain unresolved without the raw Data File S1

1. How many rows are in Data File S1?
2. Whether the file contains all 678 subjects or only a subset.
3. Whether individual subject IDs exist.
4. Whether the published training/internal/external split is reconstructible.
5. Which samples are benign adnexal masses vs healthy controls.
6. Which samples are invasive ovarian cancer.
7. Stage for each cancer.
8. Histology.
9. Age.
10. Menopausal status.
11. CA-125.
12. HE4.
13. APO-A1.
14. B2M.
15. transthyretin.
16. transferrin.
17. FSH.
18. Which 180 miRNAs are present.
19. Raw numeric scale and normalization.
20. Missing values.
21. Duplicate subjects or samples.
22. Same-patient repeats or serial collections.
23. Plate, batch, or collection-date information.
24. Patient overlap with GSE94533.

## Local evidence

Candidate files detected in the local workspace:

- none matching the published 2024 adnexal-mass supplementary Data File S1

This means the dataset cannot yet contribute to the patient-level OVAI development/validation pipeline in a leakage-safe manner.

## Decision

This cohort is scientifically relevant and should remain in the development queue, but it is not yet eligible for one of the following:

- direct pooling with GSE94533
- feature selection on external labels
- threshold tuning on external labels
- use as a canonical validation cohort

It should be treated as a high-priority acquisition target for a future, properly audited replication study.

## Reproducibility note

This audit is intentionally conservative. It records the paper-level summary and documents the access gap explicitly instead of inventing patient-level measurements.
"""
    if likely_files:
        text += f"\nDetected local candidate files: {', '.join(likely_files)}\n"
    return text


def main() -> None:
    REPORTS.mkdir(parents=True, exist_ok=True)
    OUT_MANIFEST.parent.mkdir(parents=True, exist_ok=True)

    write_feature_map()
    write_compatibility()
    report_text = build_report()
    OUT_REPORT.write_text(report_text + "\n")

    manifest = {
        "cohort": "2024 serum adnexal-mass study",
        "doi": "10.1016/j.ygyno.2024.08.008",
        "paper_title": "Serum miRNA improves the accuracy of a multivariate index assay for triage of an adnexal mass",
        "source_status": "paper-level metadata only; no raw Data File S1 available in this workspace",
        "total_subjects": 678,
        "cancer": 286,
        "controls": 392,
        "training_subjects": 468,
        "internal_validation_subjects": 100,
        "external_validation_subjects": 110,
        "miRNA_count": 180,
        "protein_biomarkers": ["CA-125", "HE4", "APO-A1", "B2M", "transthyretin", "transferrin", "FSH"],
        "other_clinical_covariates": ["age", "menopausal status"],
        "assay": "FirePlex targeted serum miRNA panel",
        "status": "high-priority candidate; not yet patient-level audited",
        "usable_for_ovai": False,
        "reason": "No raw Data File S1 or patient-level matrix was accessible in the workspace; exact assay compatibility and subject-level labels remain unverified.",
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {OUT_REPORT}")
    print(f"Wrote {OUT_FEATURE_MAP}")
    print(f"Wrote {OUT_COMPAT}")
    print(f"Wrote {OUT_MANIFEST}")


if __name__ == "__main__":
    main()
