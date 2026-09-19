#!/usr/bin/env python3

"""
Extract the exact 14-miRNA feature set from the original eLife
GSE94533 neural-network PMML and compare it against our cleaned
179-patient dataset.

This script DOES NOT train a new model.

It establishes:

    original PMML
          |
          v
    exact input features
          |
          v
    cleaned GSE94533 matrix
          |
          v
    availability / name-matching audit

The purpose is to prevent accidental "near reproduction" where
similar-looking miRNA names are substituted for the original inputs.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from xml.etree import ElementTree as ET

import pandas as pd


# ================================================================
# PATHS
# ================================================================

ROOT = Path.home() / "ovarian-ai"

CODE_ROOT = (
    ROOT
    / "data"
    / "external"
    / "elife-28932-code1"
)

MATRIX_PATH = (
    ROOT
    / "data"
    / "processed"
    / "gse94533_mirna_tpm.csv.gz"
)

OUTPUT_DIR = (
    ROOT
    / "reports"
    / "tables"
)

FEATURE_OUTPUT = (
    OUTPUT_DIR
    / "original_14_mirna_features.json"
)

AUDIT_OUTPUT = (
    OUTPUT_DIR
    / "original_14_mirna_feature_audit.json"
)


# ================================================================
# NORMALIZATION
# ================================================================

def normalize_mirna_name(value: object) -> str:

    if value is None:
        return ""

    value = str(value).strip().lower()

    # Normalize unicode dashes to ordinary hyphens.
    value = (
        value
        .replace("–", "-")
        .replace("—", "-")
        .replace("−", "-")
    )

    # Remove whitespace.
    value = re.sub(r"\s+", "", value)

    return value


# ================================================================
# FIND PMML
# ================================================================

def find_pmml() -> Path:

    candidates = sorted(
        CODE_ROOT.rglob("*.xml")
    )

    if not candidates:

        raise RuntimeError(
            "No XML/PMML file was found inside:\n"
            f"{CODE_ROOT}\n\n"
            "The original eLife source-code archive may not have "
            "been extracted correctly."
        )

    print()
    print("=" * 90)
    print("PMML FILES FOUND")
    print("=" * 90)

    for path in candidates:
        print(f"  {path}")

    # Prefer Source-code1.xml.
    for path in candidates:

        if path.name.lower() == "source-code1.xml":
            return path

    return candidates[0]


# ================================================================
# EXTRACT PMML INFORMATION
# ================================================================

def extract_pmml_features(
    pmml_path: Path,
) -> dict:

    print()
    print("=" * 90)
    print("READING ORIGINAL eLIFE PMML")
    print("=" * 90)

    print()
    print(f"PMML: {pmml_path}")

    tree = ET.parse(
        pmml_path
    )

    root = tree.getroot()

    # ------------------------------------------------------------
    # Find all fields referenced by the neural network.
    # ------------------------------------------------------------

    neural_networks = []

    for element in root.iter():

        if element.tag.endswith(
            "NeuralNetwork"
        ):

            neural_networks.append(
                element
            )

    print()
    print(
        f"Neural networks found: "
        f"{len(neural_networks)}"
    )

    if not neural_networks:

        raise RuntimeError(
            "No NeuralNetwork element was found in the PMML."
        )

    for i, network in enumerate(
        neural_networks,
        1,
    ):

        print()
        print(
            f"Network {i}: "
            f"modelName={network.attrib.get('modelName')!r}"
        )

    # The original sequencing network is the one with 14 inputs and
    # seven hidden neurons. We therefore inspect all networks and
    # identify candidate input fields.

    network_info = []

    for network in neural_networks:

        name = network.attrib.get(
            "modelName",
            "",
        )

        fields = []

        for element in network.iter():

            # PMML namespace prefixes vary, so inspect local name.
            local = (
                element.tag
                .split("}")[-1]
            )

            if local == "DerivedField":

                field_name = (
                    element.attrib.get("name")
                )

                if field_name:
                    fields.append(
                        field_name
                    )

            elif local == "MiningField":

                field_name = (
                    element.attrib.get("name")
                )

                if field_name:
                    fields.append(
                        field_name
                    )

        network_info.append(
            {
                "model_name": name,
                "fields": list(
                    dict.fromkeys(fields)
                ),
            }
        )

    # ------------------------------------------------------------
    # Search the entire XML text for all likely miRNA names.
    # ------------------------------------------------------------

    xml_text = pmml_path.read_text(
        encoding="utf-8",
        errors="replace",
    )

    mirna_pattern = re.compile(
        r"(?:hsa[-_]?)?"
        r"miR[-_]"
        r"[0-9A-Za-z]+"
        r"(?:[-_][0-9A-Za-z]+)*"
        r"(?:\*)?",
        re.IGNORECASE,
    )

    raw_matches = mirna_pattern.findall(
        xml_text
    )

    mirna_matches = []

    seen = set()

    for match in raw_matches:

        normalized = normalize_mirna_name(
            match
        )

        if normalized and normalized not in seen:

            seen.add(normalized)

            mirna_matches.append(
                match
            )

    print()
    print("miRNA-like identifiers found in PMML:")

    for feature in mirna_matches:

        print(
            f"  {feature}"
        )

    # ------------------------------------------------------------
    # Extract names from <DerivedField> elements as well.
    # ------------------------------------------------------------

    all_named_fields = []

    for element in root.iter():

        local = (
            element.tag
            .split("}")[-1]
        )

        if local in {
            "DerivedField",
            "MiningField",
            "OutputField",
        }:

            name = (
                element.attrib.get("name")
                or element.attrib.get("field")
            )

            if name:
                all_named_fields.append(
                    name
                )

    # Identify likely input features among PMML field names.
    feature_fields = []

    for field in all_named_fields:

        normalized = normalize_mirna_name(
            field
        )

        if (
            "mir" in normalized
            or "miR".lower() in normalized
        ):

            feature_fields.append(
                field
            )

    feature_fields = list(
        dict.fromkeys(
            feature_fields
        )
    )

    print()
    print(
        "miRNA-containing PMML fields:"
    )

    for field in feature_fields:
        print(
            f"  {field}"
        )

    # ------------------------------------------------------------
    # Save all discovered information.
    # ------------------------------------------------------------

    result = {
        "pmml": str(pmml_path),
        "neural_networks": network_info,
        "mirna_like_identifiers": mirna_matches,
        "mirna_named_fields": feature_fields,
    }

    return result


# ================================================================
# LOAD CLEAN MATRIX
# ================================================================

def load_clean_matrix() -> pd.DataFrame:

    if not MATRIX_PATH.exists():

        raise RuntimeError(
            "Clean matrix not found:\n"
            f"{MATRIX_PATH}\n\n"
            "Run build_gse94533_clean_dataset.py first."
        )

    matrix = pd.read_csv(
        MATRIX_PATH,
        index_col=0,
    )

    print()
    print("=" * 90)
    print("CLEAN MATRIX")
    print("=" * 90)

    print()
    print(
        f"Patients: {matrix.shape[0]}"
    )

    print(
        f"Features: {matrix.shape[1]}"
    )

    return matrix


# ================================================================
# MATCH FEATURES
# ================================================================

def audit_features(
    discovered: dict,
    matrix: pd.DataFrame,
) -> dict:

    print()
    print("=" * 90)
    print("FEATURE MATCHING AUDIT")
    print("=" * 90)

    matrix_lookup = {}

    for feature in matrix.columns:

        normalized = normalize_mirna_name(
            feature
        )

        matrix_lookup.setdefault(
            normalized,
            [],
        ).append(
            feature
        )

    # Use all likely PMML features.
    original_fields = (
        discovered["mirna_named_fields"]
    )

    print()
    print(
        f"PMML miRNA-containing fields: "
        f"{len(original_fields)}"
    )

    matches = {}
    missing = {}

    for original in original_fields:

        normalized = normalize_mirna_name(
            original
        )

        candidates = matrix_lookup.get(
            normalized,
            [],
        )

        if candidates:

            matches[original] = candidates

        else:

            missing[original] = []

    print()
    print("MATCHED FEATURES:")

    for original, candidates in matches.items():

        print(
            f"  {original!r}"
            f" -> {candidates}"
        )

    print()
    print("UNMATCHED FEATURES:")

    if missing:

        for original in missing:
            print(
                f"  {original!r}"
            )

    else:

        print(
            "  None"
        )

    result = {
        "matched": matches,
        "missing": list(
            missing.keys()
        ),
        "matrix_feature_count": int(
            matrix.shape[1]
        ),
    }

    return result


# ================================================================
# MAIN
# ================================================================

def main() -> None:

    print("=" * 90)
    print(
        "ORIGINAL 14-miRNA MODEL FEATURE AUDIT"
    )
    print("=" * 90)

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    pmml_path = find_pmml()

    discovered = extract_pmml_features(
        pmml_path
    )

    matrix = load_clean_matrix()

    audit = audit_features(
        discovered,
        matrix,
    )

    FEATURE_OUTPUT.write_text(
        json.dumps(
            discovered,
            indent=2,
        )
    )

    AUDIT_OUTPUT.write_text(
        json.dumps(
            audit,
            indent=2,
        )
    )

    print()
    print("=" * 90)
    print("AUDIT COMPLETE")
    print("=" * 90)

    print()
    print(
        f"Feature information:"
        f"\n  {FEATURE_OUTPUT}"
    )

    print()
    print(
        f"Matching audit:"
        f"\n  {AUDIT_OUTPUT}"
    )

    print()
    print(
        "IMPORTANT:"
    )

    print(
        """
Do not train anything from this output yet.

The purpose of this stage is to determine the exact 14-miRNA
inputs used by the historical neural network and verify that those
features exist in our cleaned GSE94533 matrix.

After this audit, we will construct a benchmark using the original
feature definitions rather than accidentally creating a different
model and calling it a reproduction.
"""
    )


if __name__ == "__main__":
    main()
