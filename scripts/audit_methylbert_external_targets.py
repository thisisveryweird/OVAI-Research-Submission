#!/usr/bin/env python3

"""Audit proposed external methylation targets before frozen prediction."""

import json
import re
from pathlib import Path

import requests
import xenaPython as xena


PROJECT = Path(__file__).resolve().parents[1]
MODEL_METADATA = PROJECT / "models" / "methylbert_ridge_cancer_vs_healthy.json"
OUT = PROJECT / "reports" / "tables" / "methylbert_external_target_audit.json"


def main():
    model = json.loads(MODEL_METADATA.read_text())
    training_features = set(model["features"])
    gse_text = requests.get(
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144365&targ=self&form=text&view=quick",
        timeout=120,
    ).text
    tcga_features = set(str(value) for value in xena.dataset_field(
        "https://tcga.xenahubs.net",
        "TCGA.OV.sampleMap/HumanMethylation450",
    ))
    audit = {
        "training_features": len(training_features),
        "gse144365": {
            "series_type": next((line for line in gse_text.splitlines() if line.startswith("!Series_type")), ""),
            "organism": next((line for line in gse_text.splitlines() if line.startswith("!Series_platform_organism")), ""),
            "valid_ovarian_methylation_target": bool(re.search("ovarian|methyl", gse_text, re.I)),
            "reason": "GSE144365 is a Pseudomonas aeruginosa transcriptome study, not ovarian methylation.",
        },
        "tcga_ov_humanmethylation450": {
            "features": len(tcga_features),
            "overlap_with_methylbert_coordinates": len(training_features & tcga_features),
            "valid_exact_feature_target": bool(training_features & tcga_features),
            "reason": "TCGA uses cg probe IDs and tissue; MethylBERT uses chr_coordinate features and blood samples.",
        },
        "gsa_hra011451": {
            "frozen_prediction_performed": False,
            "reason": "No independently downloadable, feature-aligned validation matrix was available.",
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()