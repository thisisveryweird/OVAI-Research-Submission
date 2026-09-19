#!/usr/bin/env python3

"""
Deep audit of the original eLife GSE94533 neural-network PMML.

Purpose
-------
The previous audit established the exact 14 miRNA feature names.

This script now determines:

    1. exact PMML input fields
    2. input transformations / derived fields
    3. normalization/scaling
    4. neural-network activation functions
    5. hidden-layer size
    6. output-layer structure
    7. every neural-network neuron weight
    8. bias/threshold values where present
    9. output probability definitions
   10. all model metadata relevant to reproducing the historical model

IMPORTANT
---------
This script DOES NOT train a model.

We are trying to reproduce the published model before modifying it.

That matters because the original paper states that the winning
neural network used miRNA expression fold-change variables and that
the data were filtered and batch-adjusted before model construction.
Simply feeding raw TPM values into a new 14-input neural network would
therefore be a NEW model, not a reproduction.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from xml.etree import ElementTree as ET


# ================================================================
# PATHS
# ================================================================

ROOT = Path.home() / "ovarian-ai"

PMML_PATH = (
    ROOT
    / "data"
    / "external"
    / "elife-28932-code1"
    / "Source-code1.xml"
)

OUTPUT_DIR = (
    ROOT
    / "reports"
    / "tables"
)

JSON_OUTPUT = (
    OUTPUT_DIR
    / "original_pmml_deep_audit.json"
)

TXT_OUTPUT = (
    OUTPUT_DIR
    / "original_pmml_deep_audit.txt"
)


# ================================================================
# XML HELPERS
# ================================================================

def local_name(tag: str) -> str:
    """Return XML tag without namespace."""
    return tag.split("}")[-1]


def children(element, name: str):
    """Return direct children matching local XML name."""
    return [
        child
        for child in list(element)
        if local_name(child.tag) == name
    ]


def descendants(element, name: str):
    """Return all descendants matching local XML name."""
    return [
        child
        for child in element.iter()
        if local_name(child.tag) == name
    ]


def child_text(element, name: str):
    """Get direct child's text if present."""
    for child in list(element):
        if local_name(child.tag) == name:
            return child.text
    return None


def safe_float(value):
    """Convert XML numeric values to float where possible."""
    if value is None:
        return None

    try:
        return float(value)
    except Exception:
        return None


# ================================================================
# EXPRESSION PARSER
# ================================================================

def xml_expression_to_dict(element):
    """
    Convert PMML expression trees into a readable nested dictionary.

    PMML preprocessing may contain:
        FieldRef
        Constant
        Apply
        NormContinuous
        NormDiscrete
        Discretize
        MapValues
    """

    if element is None:
        return None

    name = local_name(element.tag)

    result = {
        "type": name,
    }

    if element.attrib:
        result["attributes"] = dict(
            element.attrib
        )

    if element.text and element.text.strip():
        result["text"] = element.text.strip()

    child_expr = []

    for child in list(element):

        child_expr.append(
            xml_expression_to_dict(child)
        )

    if child_expr:
        result["children"] = child_expr

    return result


# ================================================================
# PMML ROOT INFORMATION
# ================================================================

def read_pmml(pmml_path: Path):

    if not pmml_path.exists():

        raise RuntimeError(
            f"PMML file not found:\n{pmml_path}"
        )

    print()
    print("=" * 100)
    print("READING PMML")
    print("=" * 100)

    print()
    print(pmml_path)

    tree = ET.parse(
        pmml_path
    )

    return tree.getroot()


# ================================================================
# HEADER / MODEL INFORMATION
# ================================================================

def extract_header(root):

    header_elements = descendants(
        root,
        "Header",
    )

    result = {}

    if not header_elements:
        return result

    header = header_elements[0]

    result["attributes"] = dict(
        header.attrib
    )

    descriptions = descendants(
        header,
        "InlineTable",
    )

    # Basic direct text descriptions.
    text_values = []

    for element in header.iter():

        if element.text and element.text.strip():

            text = element.text.strip()

            if len(text) < 1000:
                text_values.append(text)

    result["text"] = text_values

    return result


# ================================================================
# DATA DICTIONARY
# ================================================================

def extract_data_dictionary(root):

    result = []

    dictionaries = descendants(
        root,
        "DataDictionary",
    )

    if not dictionaries:
        return result

    dictionary = dictionaries[0]

    for field in children(
        dictionary,
        "DataField",
    ):

        result.append(
            {
                "name": field.attrib.get("name"),
                "optype": field.attrib.get("optype"),
                "dataType": field.attrib.get("dataType"),
            }
        )

    return result


# ================================================================
# TRANSFORMATION DICTIONARY
# ================================================================

def extract_transformations(root):

    result = []

    transformations = descendants(
        root,
        "TransformationDictionary",
    )

    if not transformations:

        return result

    dictionary = transformations[0]

    for field in children(
        dictionary,
        "DerivedField",
    ):

        name = field.attrib.get(
            "name"
        )

        result.append(
            {
                "name": name,
                "optype": field.attrib.get(
                    "optype"
                ),
                "dataType": field.attrib.get(
                    "dataType"
                ),
                "expression": (
                    [
                        xml_expression_to_dict(
                            child
                        )
                        for child in list(field)
                    ]
                    or None
                ),
            }
        )

    return result


# ================================================================
# MINING SCHEMA
# ================================================================

def extract_mining_schema(root):

    result = []

    schemas = descendants(
        root,
        "MiningSchema",
    )

    if not schemas:
        return result

    # The PMML contains one model here.
    schema = schemas[0]

    for field in children(
        schema,
        "MiningField",
    ):

        result.append(
            {
                "name": field.attrib.get(
                    "name"
                ),
                "usageType": field.attrib.get(
                    "usageType"
                ),
                "optype": field.attrib.get(
                    "optype"
                ),
                "missingValueReplacement": field.attrib.get(
                    "missingValueReplacement"
                ),
                "missingValueTreatment": field.attrib.get(
                    "missingValueTreatment"
                ),
            }
        )

    return result


# ================================================================
# LOCAL DERIVED FIELDS
# ================================================================

def extract_all_derived_fields(root):

    result = []

    for field in descendants(
        root,
        "DerivedField",
    ):

        result.append(
            {
                "name": field.attrib.get(
                    "name"
                ),
                "optype": field.attrib.get(
                    "optype"
                ),
                "dataType": field.attrib.get(
                    "dataType"
                ),
                "expression": (
                    [
                        xml_expression_to_dict(
                            child
                        )
                        for child in list(field)
                    ]
                ),
            }
        )

    return result


# ================================================================
# NEURAL NETWORK
# ================================================================

def extract_neural_network(root):

    networks = descendants(
        root,
        "NeuralNetwork",
    )

    if not networks:

        raise RuntimeError(
            "No NeuralNetwork found in PMML."
        )

    if len(networks) > 1:

        print(
            f"WARNING: {len(networks)} neural networks found."
        )

    network = None

    for candidate in networks:

        model_name = candidate.attrib.get(
            "modelName",
            "",
        )

        if "14-7-2" in model_name:
            network = candidate
            break

    if network is None:
        network = networks[0]

    result = {

        "attributes": dict(
            network.attrib
        ),

        "neural_inputs": [],

        "neural_outputs": [],

        "neural_layers": [],

        "neurons": [],

        "connections": [],
    }

    # ------------------------------------------------------------
    # NeuralInputs
    # ------------------------------------------------------------

    input_containers = children(
        network,
        "NeuralInputs",
    )

    for container in input_containers:

        for neural_input in children(
            container,
            "NeuralInput",
        ):

            inp = {
                "id": neural_input.attrib.get(
                    "id"
                ),
                "attributes": dict(
                    neural_input.attrib
                ),
            }

            expressions = []

            for child in list(
                neural_input
            ):

                expressions.append(
                    xml_expression_to_dict(
                        child
                    )
                )

            inp["expressions"] = expressions

            result[
                "neural_inputs"
            ].append(inp)

    # ------------------------------------------------------------
    # NeuralOutputs
    # ------------------------------------------------------------

    output_containers = children(
        network,
        "NeuralOutputs",
    )

    for container in output_containers:

        for output in children(
            container,
            "NeuralOutput",
        ):

            result[
                "neural_outputs"
            ].append(
                {
                    "outputNeuron": output.attrib.get(
                        "outputNeuron"
                    ),
                    "output_value": output.attrib.get(
                        "outputValue"
                    ),
                    "resultFeature": output.attrib.get(
                        "resultFeature"
                    ),
                }
            )

    # ------------------------------------------------------------
    # NeuralLayers
    # ------------------------------------------------------------

    for layer in children(
        network,
        "NeuralLayer",
    ):

        layer_info = {
            "activationFunction": layer.attrib.get(
                "activationFunction"
            ),
            "numberOfNeurons": len(
                children(
                    layer,
                    "Neuron",
                )
            ),
            "neurons": [],
        }

        for neuron in children(
            layer,
            "Neuron",
        ):

            neuron_info = {
                "id": neuron.attrib.get(
                    "id"
                ),
                "bias": safe_float(
                    neuron.attrib.get(
                        "bias"
                    )
                ),
                "width": safe_float(
                    neuron.attrib.get(
                        "width"
                    )
                ),
                "connections": [],
            }

            # Some PMML exporters place threshold as an attribute.
            if "threshold" in neuron.attrib:

                neuron_info[
                    "threshold"
                ] = safe_float(
                    neuron.attrib[
                        "threshold"
                    ]
                )

            for connection in children(
                neuron,
                "Con",
            ):

                conn = {
                    "from": connection.attrib.get(
                        "from"
                    ),
                    "weight": safe_float(
                        connection.attrib.get(
                            "weight"
                        )
                    ),
                    "attributes": dict(
                        connection.attrib
                    ),
                }

                neuron_info[
                    "connections"
                ].append(conn)

                result[
                    "connections"
                ].append(
                    {
                        "layer_activation": layer.attrib.get(
                            "activationFunction"
                        ),
                        "neuron": neuron.attrib.get(
                            "id"
                        ),
                        "from": connection.attrib.get(
                            "from"
                        ),
                        "weight": safe_float(
                            connection.attrib.get(
                                "weight"
                            )
                        ),
                    }
                )

            layer_info[
                "neurons"
            ].append(
                neuron_info
            )

            result[
                "neurons"
            ].append(
                {
                    "layer_activation": layer.attrib.get(
                        "activationFunction"
                    ),
                    **neuron_info,
                }
            )

        result[
            "neural_layers"
        ].append(
            layer_info
        )

    return result


# ================================================================
# OUTPUT FIELDS
# ================================================================

def extract_output_fields(root):

    result = []

    for output in descendants(
        root,
        "Output",
    ):

        for field in children(
            output,
            "OutputField",
        ):

            result.append(
                {
                    "name": field.attrib.get(
                        "name"
                    ),
                    "feature": field.attrib.get(
                        "feature"
                    ),
                    "value": field.attrib.get(
                        "value"
                    ),
                    "targetField": field.attrib.get(
                        "targetField"
                    ),
                    "rank": field.attrib.get(
                        "rank"
                    ),
                }
            )

    return result


# ================================================================
# SCORECARD / REGRESSION / OTHER STRUCTURES
# ================================================================

def extract_other_model_elements(root):

    names = [
        "RegressionModel",
        "Scorecard",
        "TreeModel",
        "SupportVectorMachineModel",
        "NaiveBayesModel",
        "GeneralRegressionModel",
        "ClusteringModel",
    ]

    result = {}

    for name in names:

        count = len(
            descendants(
                root,
                name,
            )
        )

        if count:
            result[name] = count

    return result


# ================================================================
# IDENTIFY 14 INPUTS
# ================================================================

def identify_model_inputs(
    neural_network,
):

    inputs = []

    for inp in neural_network[
        "neural_inputs"
    ]:

        expressions = inp.get(
            "expressions",
            [],
        )

        inputs.append(
            {
                "id": inp.get(
                    "id"
                ),
                "expressions": expressions,
            }
        )

    return inputs


# ================================================================
# WEIGHT SUMMARY
# ================================================================

def summarize_weights(
    neural_network,
):

    weights = []

    for connection in neural_network[
        "connections"
    ]:

        weight = connection.get(
            "weight"
        )

        if weight is not None:

            weights.append(weight)

    if not weights:

        return {}

    return {
        "count": len(weights),
        "min": min(weights),
        "max": max(weights),
        "mean": sum(weights)
        / len(weights),
        "abs_max": max(
            abs(x)
            for x in weights
        ),
        "nonzero": sum(
            abs(x) > 1e-12
            for x in weights
        ),
    }


# ================================================================
# PRINT HUMAN-READABLE REPORT
# ================================================================

def format_report(
    audit,
):

    lines = []

    lines.append(
        "=" * 100
    )

    lines.append(
        "ORIGINAL eLIFE GSE94533 PMML DEEP AUDIT"
    )

    lines.append(
        "=" * 100
    )

    # ------------------------------------------------------------
    # Model information.
    # ------------------------------------------------------------

    nn = audit[
        "neural_network"
    ]

    lines.append("")
    lines.append(
        "MODEL"
    )
    lines.append(
        "-" * 100
    )

    for key, value in nn[
        "attributes"
    ].items():

        lines.append(
            f"{key}: {value}"
        )

    lines.append("")
    lines.append(
        "NEURAL INPUTS"
    )

    for inp in nn[
        "neural_inputs"
    ]:

        lines.append(
            f"\n{inp['id']}"
        )

        for expression in inp[
            "expressions"
        ]:

            lines.append(
                json.dumps(
                    expression,
                    indent=2,
                )
            )

    lines.append("")
    lines.append(
        "NEURAL OUTPUTS"
    )

    for output in nn[
        "neural_outputs"
    ]:

        lines.append(
            str(output)
        )

    lines.append("")
    lines.append(
        "NEURAL LAYERS"
    )

    for i, layer in enumerate(
        nn["neural_layers"],
        1,
    ):

        lines.append(
            f"\nLayer {i}"
        )

        lines.append(
            f"Activation: "
            f"{layer.get('activationFunction')}"
        )

        lines.append(
            f"Neurons: "
            f"{layer.get('numberOfNeurons')}"
        )

        for neuron in layer[
            "neurons"
        ]:

            lines.append(
                f"\n  Neuron "
                f"{neuron.get('id')}"
            )

            lines.append(
                f"    bias = "
                f"{neuron.get('bias')}"
            )

            lines.append(
                f"    width = "
                f"{neuron.get('width')}"
            )

            if "threshold" in neuron:

                lines.append(
                    f"    threshold = "
                    f"{neuron.get('threshold')}"
                )

            for conn in neuron[
                "connections"
            ]:

                lines.append(
                    f"    from "
                    f"{conn.get('from')}: "
                    f"weight="
                    f"{conn.get('weight')}"
                )

    # ------------------------------------------------------------
    # Transformations.
    # ------------------------------------------------------------

    lines.append("")
    lines.append(
        "TRANSFORMATION DICTIONARY"
    )
    lines.append(
        "-" * 100
    )

    for derived in audit[
        "transformation_dictionary"
    ]:

        lines.append(
            f"\n{derived['name']}"
        )

        lines.append(
            f"optype={derived.get('optype')}"
        )

        lines.append(
            f"dataType={derived.get('dataType')}"
        )

        lines.append(
            json.dumps(
                derived.get(
                    "expression"
                ),
                indent=2,
            )
        )

    # ------------------------------------------------------------
    # Local derived fields.
    # ------------------------------------------------------------

    lines.append("")
    lines.append(
        "ALL DERIVED FIELDS"
    )
    lines.append(
        "-" * 100
    )

    for derived in audit[
        "derived_fields"
    ]:

        lines.append(
            f"\n{derived['name']}"
        )

        lines.append(
            json.dumps(
                derived.get(
                    "expression"
                ),
                indent=2,
            )
        )

    # ------------------------------------------------------------
    # Mining schema.
    # ------------------------------------------------------------

    lines.append("")
    lines.append(
        "MINING SCHEMA"
    )
    lines.append(
        "-" * 100
    )

    for field in audit[
        "mining_schema"
    ]:

        lines.append(
            str(field)
        )

    # ------------------------------------------------------------
    # Output fields.
    # ------------------------------------------------------------

    lines.append("")
    lines.append(
        "OUTPUT FIELDS"
    )
    lines.append(
        "-" * 100
    )

    for field in audit[
        "output_fields"
    ]:

        lines.append(
            str(field)
        )

    # ------------------------------------------------------------
    # Weight summary.
    # ------------------------------------------------------------

    lines.append("")
    lines.append(
        "WEIGHT SUMMARY"
    )
    lines.append(
        "-" * 100
    )

    lines.append(
        str(
            audit[
                "weight_summary"
            ]
        )
    )

    # ------------------------------------------------------------
    # Other model elements.
    # ------------------------------------------------------------

    lines.append("")
    lines.append(
        "OTHER MODEL ELEMENTS"
    )
    lines.append(
        "-" * 100
    )

    lines.append(
        str(
            audit[
                "other_model_elements"
            ]
        )
    )

    return "\n".join(
        lines
    )


# ================================================================
# MAIN
# ================================================================

def main():

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    root = read_pmml(
        PMML_PATH
    )

    audit = {

        "pmml_path": str(
            PMML_PATH
        ),

        "header": extract_header(
            root
        ),

        "data_dictionary": (
            extract_data_dictionary(
                root
            )
        ),

        "transformation_dictionary": (
            extract_transformations(
                root
            )
        ),

        "derived_fields": (
            extract_all_derived_fields(
                root
            )
        ),

        "mining_schema": (
            extract_mining_schema(
                root
            )
        ),

        "neural_network": (
            extract_neural_network(
                root
            )
        ),

        "output_fields": (
            extract_output_fields(
                root
            )
        ),

        "other_model_elements": (
            extract_other_model_elements(
                root
            )
        ),
    }

    audit[
        "identified_neural_inputs"
    ] = identify_model_inputs(
        audit[
            "neural_network"
        ]
    )

    audit[
        "weight_summary"
    ] = summarize_weights(
        audit[
            "neural_network"
        ]
    )

    # ------------------------------------------------------------
    # Save JSON.
    # ------------------------------------------------------------

    JSON_OUTPUT.write_text(
        json.dumps(
            audit,
            indent=2,
        )
    )

    # ------------------------------------------------------------
    # Save human-readable text.
    # ------------------------------------------------------------

    report = format_report(
        audit
    )

    TXT_OUTPUT.write_text(
        report
    )

    print()
    print(
        report
    )

    print()
    print("=" * 100)
    print(
        "FILES SAVED"
    )
    print("=" * 100)

    print()
    print(
        JSON_OUTPUT
    )

    print(
        TXT_OUTPUT
    )

    print()
    print("=" * 100)
    print(
        "NO MODEL TRAINING PERFORMED"
    )
    print("=" * 100)

    print(
        """
Next step:

Use this audit to determine whether the PMML contains enough
information to implement an exact forward pass.

Specifically, we want to know:

    raw feature
        ↓
    PMML transformation
        ↓
    neural-network input
        ↓
    14 → 7 → 2 network
        ↓
    probability of ovarian cancer

Only after that should we compare the historical PMML model against
our current dataset.

Then we can construct the modern baseline and, eventually, the
multimodal genomic + symptomatic model.
"""
    )


if __name__ == "__main__":
    main()
