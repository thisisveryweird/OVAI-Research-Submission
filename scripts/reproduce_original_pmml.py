#!/usr/bin/env python3

"""
Reproduce the original eLife GSE94533 14-miRNA neural network.

IMPORTANT
---------
This script DOES NOT train anything.

It implements the already-published PMML model directly from:

    original_pmml_deep_audit.json
    +
    gse94533_14mirna_benchmark.csv.gz

The model is:

    14 inputs
        ->
    7 tanh hidden neurons
        ->
    2 tanh output neurons

Input preprocessing is taken from the PMML:

    x_normalized = x * scale + shift

The PMML's two output neurons are:

    neuron 21 = Ovarian Cancer
    neuron 22 = Controls or Borderline

For classification, PMML uses the output neuron with the largest
activation to determine the predicted class.

For evaluation we additionally retain:

    cancer_score = activation(neuron 21)

This is a continuous model score.

We deliberately DO NOT call:

    (score + 1) / 2

or a softmax transformation an "original probability", because
the supplied PMML does not specify such a transformation for the
two tanh output neurons.

This is a research reproduction only. It is not a clinical
diagnostic system.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


# ================================================================
# PATHS
# ================================================================

ROOT = Path.home() / "ovarian-ai"

DATASET_PATH = (
    ROOT
    / "data"
    / "processed"
    / "gse94533_14mirna_benchmark.csv.gz"
)

AUDIT_PATH = (
    ROOT
    / "reports"
    / "tables"
    / "original_pmml_deep_audit.json"
)

OUTPUT_DIR = (
    ROOT
    / "reports"
    / "tables"
)

PREDICTIONS_PATH = (
    OUTPUT_DIR
    / "original_pmml_predictions.csv"
)

RESULTS_PATH = (
    OUTPUT_DIR
    / "original_pmml_reproduction_results.json"
)

TRACE_PATH = (
    OUTPUT_DIR
    / "original_pmml_forward_pass_trace.csv"
)


# ================================================================
# EXACT HISTORICAL INPUT ORDER
# ================================================================

FEATURES = [
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


# ================================================================
# EXPECTED NETWORK STRUCTURE
# ================================================================

HIDDEN_IDS = [
    "14",
    "15",
    "16",
    "17",
    "18",
    "19",
    "20",
]

OUTPUT_IDS = [
    "21",
    "22",
]


# ================================================================
# ACTIVATION FUNCTIONS
# ================================================================

def tanh(x: np.ndarray) -> np.ndarray:
    return np.tanh(x)


# ================================================================
# LOAD DATA
# ================================================================

def load_dataset() -> pd.DataFrame:

    if not DATASET_PATH.exists():

        raise RuntimeError(
            f"Dataset not found:\n{DATASET_PATH}"
        )

    df = pd.read_csv(
        DATASET_PATH
    )

    required = [
        "sample_id",
        "diagnosis",
        "binary_label",
    ] + FEATURES

    missing = [
        col
        for col in required
        if col not in df.columns
    ]

    if missing:

        raise RuntimeError(
            "Dataset missing required columns:\n"
            + "\n".join(
                f"  {x}"
                for x in missing
            )
        )

    return df


def load_audit():

    if not AUDIT_PATH.exists():

        raise RuntimeError(
            f"PMML audit not found:\n{AUDIT_PATH}"
        )

    return json.loads(
        AUDIT_PATH.read_text()
    )


# ================================================================
# EXTRACT PMML NORMALIZATION
# ================================================================

def extract_normalization(audit):

    """
    Extract the PMML normalization parameters.

    For STATISTICA-style PMML seen here:

        normalized = raw * scale + shift

    The LinearNorm endpoints independently confirm the same mapping.
    """

    neural_inputs = (
        audit[
            "neural_network"
        ][
            "neural_inputs"
        ]
    )

    normalization = {}

    for feature, neural_input in zip(
        FEATURES,
        neural_inputs,
    ):

        expressions = (
            neural_input.get(
                "expressions",
                []
            )
        )

        if len(expressions) != 1:

            raise RuntimeError(
                f"Expected one expression for {feature}, "
                f"found {len(expressions)}."
            )

        expression = expressions[0]

        if expression.get(
            "type"
        ) != "DerivedField":

            raise RuntimeError(
                f"{feature}: expected DerivedField."
            )

        children = expression.get(
            "children",
            []
        )

        if len(children) != 1:

            raise RuntimeError(
                f"{feature}: expected one DerivedField child."
            )

        norm = children[0]

        if norm.get(
            "type"
        ) != "NormContinuous":

            raise RuntimeError(
                f"{feature}: expected NormContinuous."
            )

        attrs = norm.get(
            "attributes",
            {}
        )

        field = attrs.get(
            "field"
        )

        if field != feature:

            raise RuntimeError(
                f"Feature mismatch: "
                f"expected {feature}, "
                f"PMML says {field}."
            )

        shift = float(
            attrs["shift"]
        )

        scale = float(
            attrs["scale"]
        )

        linear_norms = [
            child
            for child in norm.get(
                "children",
                []
            )
            if child.get(
                "type"
            ) == "LinearNorm"
        ]

        if len(linear_norms) != 2:

            raise RuntimeError(
                f"{feature}: expected exactly two LinearNorm elements."
            )

        orig0 = float(
            linear_norms[0][
                "attributes"
            ][
                "orig"
            ]
        )

        norm0 = float(
            linear_norms[0][
                "attributes"
            ][
                "norm"
            ]
        )

        orig1 = float(
            linear_norms[1][
                "attributes"
            ][
                "orig"
            ]
        )

        norm1 = float(
            linear_norms[1][
                "attributes"
            ][
                "norm"
            ]
        )

        # --------------------------------------------------------
        # Independent verification:
        #
        # y = (x - orig0)/(orig1 - orig0)
        #
        # should equal:
        #
        # y = x*scale + shift
        # --------------------------------------------------------

        expected_scale = (
            norm1 - norm0
        ) / (
            orig1 - orig0
        )

        expected_shift = (
            norm0
            - expected_scale * orig0
        )

        if not math.isclose(
            scale,
            expected_scale,
            rel_tol=1e-8,
            abs_tol=1e-10,
        ):

            raise RuntimeError(
                f"{feature}: PMML scale does not agree with LinearNorm."
            )

        if not math.isclose(
            shift,
            expected_shift,
            rel_tol=1e-8,
            abs_tol=1e-10,
        ):

            raise RuntimeError(
                f"{feature}: PMML shift does not agree with LinearNorm."
            )

        normalization[
            feature
        ] = {
            "shift": shift,
            "scale": scale,
            "orig_min": orig0,
            "orig_max": orig1,
            "norm_min": norm0,
            "norm_max": norm1,
        }

    return normalization


# ================================================================
# EXTRACT NEURAL NETWORK
# ================================================================

def extract_network(audit):

    nn = audit[
        "neural_network"
    ]

    # ------------------------------------------------------------
    # Check architecture.
    # ------------------------------------------------------------

    layers = nn[
        "neural_layers"
    ]

    if len(layers) != 2:

        raise RuntimeError(
            f"Expected 2 layers, found {len(layers)}."
        )

    hidden_layer = layers[0]
    output_layer = layers[1]

    if (
        hidden_layer[
            "activationFunction"
        ]
        != "tanh"
    ):

        raise RuntimeError(
            "Historical hidden layer is not tanh."
        )

    if (
        output_layer[
            "activationFunction"
        ]
        != "tanh"
    ):

        raise RuntimeError(
            "Historical output layer is not tanh."
        )

    if (
        hidden_layer[
            "numberOfNeurons"
        ]
        != 7
    ):

        raise RuntimeError(
            "Historical hidden layer does not contain 7 neurons."
        )

    if (
        output_layer[
            "numberOfNeurons"
        ]
        != 2
    ):

        raise RuntimeError(
            "Historical output layer does not contain 2 neurons."
        )

    # ------------------------------------------------------------
    # Convert neurons into numerical parameter matrices.
    # ------------------------------------------------------------

    hidden_weights = np.zeros(
        (
            7,
            14,
        ),
        dtype=float,
    )

    hidden_bias = np.zeros(
        7,
        dtype=float,
    )

    for i, neuron in enumerate(
        hidden_layer[
            "neurons"
        ]
    ):

        expected_id = HIDDEN_IDS[i]

        if neuron[
            "id"
        ] != expected_id:

            raise RuntimeError(
                f"Unexpected hidden neuron order: "
                f"expected {expected_id}, "
                f"found {neuron['id']}."
            )

        hidden_bias[i] = float(
            neuron[
                "bias"
            ]
        )

        connections = neuron[
            "connections"
        ]

        if len(connections) != 14:

            raise RuntimeError(
                f"Neuron {expected_id} has "
                f"{len(connections)} inputs; expected 14."
            )

        for connection in connections:

            source = str(
                connection[
                    "from"
                ]
            )

            source_index = int(
                source
            )

            if not (
                0
                <= source_index
                < 14
            ):

                raise RuntimeError(
                    f"Invalid hidden-layer input index: {source}"
                )

            hidden_weights[
                i,
                source_index,
            ] = float(
                connection[
                    "weight"
                ]
            )

    output_weights = np.zeros(
        (
            2,
            7,
        ),
        dtype=float,
    )

    output_bias = np.zeros(
        2,
        dtype=float,
    )

    for i, neuron in enumerate(
        output_layer[
            "neurons"
        ]
    ):

        expected_id = OUTPUT_IDS[i]

        if neuron[
            "id"
        ] != expected_id:

            raise RuntimeError(
                f"Unexpected output neuron order: "
                f"expected {expected_id}, "
                f"found {neuron['id']}."
            )

        output_bias[i] = float(
            neuron[
                "bias"
            ]
        )

        connections = neuron[
            "connections"
        ]

        if len(connections) != 7:

            raise RuntimeError(
                f"Neuron {expected_id} has "
                f"{len(connections)} hidden inputs; expected 7."
            )

        for connection in connections:

            source = str(
                connection[
                    "from"
                ]
            )

            if source not in HIDDEN_IDS:

                raise RuntimeError(
                    f"Invalid output-layer source neuron: {source}"
                )

            source_index = HIDDEN_IDS.index(
                source
            )

            output_weights[
                i,
                source_index,
            ] = float(
                connection[
                    "weight"
                ]
            )

    return {
        "hidden_weights": hidden_weights,
        "hidden_bias": hidden_bias,
        "output_weights": output_weights,
        "output_bias": output_bias,
    }


# ================================================================
# FORWARD PASS
# ================================================================

def forward_pass(
    df: pd.DataFrame,
    normalization,
    network,
):

    raw = df[
        FEATURES
    ].to_numpy(
        dtype=float
    )

    if not np.isfinite(
        raw
    ).all():

        raise RuntimeError(
            "Raw 14-miRNA data contain non-finite values."
        )

    # ------------------------------------------------------------
    # PMML input normalization.
    # ------------------------------------------------------------

    normalized = np.zeros_like(
        raw,
        dtype=float,
    )

    for j, feature in enumerate(
        FEATURES
    ):

        shift = normalization[
            feature
        ][
            "shift"
        ]

        scale = normalization[
            feature
        ][
            "scale"
        ]

        normalized[
            :,
            j
        ] = (
            raw[
                :,
                j
            ]
            * scale
            + shift
        )

    # ------------------------------------------------------------
    # Hidden layer.
    #
    # z = XW^T + bias
    # h = tanh(z)
    # ------------------------------------------------------------

    hidden_linear = (
        normalized
        @ network[
            "hidden_weights"
        ].T
        + network[
            "hidden_bias"
        ]
    )

    hidden_activation = tanh(
        hidden_linear
    )

    # ------------------------------------------------------------
    # Output layer.
    # ------------------------------------------------------------

    output_linear = (
        hidden_activation
        @ network[
            "output_weights"
        ].T
        + network[
            "output_bias"
        ]
    )

    output_activation = tanh(
        output_linear
    )

    return {
        "raw": raw,
        "normalized": normalized,
        "hidden_linear": hidden_linear,
        "hidden_activation": hidden_activation,
        "output_linear": output_linear,
        "output_activation": output_activation,
    }


# ================================================================
# METRICS
# ================================================================

def compute_confusion(
    y_true,
    y_pred,
):

    y_true = np.asarray(
        y_true,
        dtype=int,
    )

    y_pred = np.asarray(
        y_pred,
        dtype=int,
    )

    tp = int(
        (
            (y_true == 1)
            & (y_pred == 1)
        ).sum()
    )

    tn = int(
        (
            (y_true == 0)
            & (y_pred == 0)
        ).sum()
    )

    fp = int(
        (
            (y_true == 0)
            & (y_pred == 1)
        ).sum()
    )

    fn = int(
        (
            (y_true == 1)
            & (y_pred == 0)
        ).sum()
    )

    sensitivity = (
        tp / (tp + fn)
        if tp + fn
        else None
    )

    specificity = (
        tn / (tn + fp)
        if tn + fp
        else None
    )

    accuracy = (
        (tp + tn)
        / len(y_true)
    )

    return {
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "accuracy": accuracy,
    }


def compute_auc(
    y_true,
    score,
):

    """
    Calculate ROC AUC without requiring sklearn.

    Uses the rank/Mann-Whitney formulation.

    AUC = probability that a randomly chosen positive receives
    a higher score than a randomly chosen negative, with ties
    receiving half credit.
    """

    y_true = np.asarray(
        y_true,
        dtype=int,
    )

    score = np.asarray(
        score,
        dtype=float,
    )

    positive = score[
        y_true == 1
    ]

    negative = score[
        y_true == 0
    ]

    if len(positive) == 0:
        return None

    if len(negative) == 0:
        return None

    wins = 0.0

    for p in positive:

        wins += (
            (p > negative)
        ).sum()

        wins += 0.5 * (
            (p == negative)
        ).sum()

    auc = (
        wins
        / (
            len(positive)
            * len(negative)
        )
    )

    return float(
        auc
    )


# ================================================================
# MAIN
# ================================================================

def main():

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    print("=" * 100)
    print("REPRODUCING ORIGINAL eLIFE GSE94533 PMML MODEL")
    print("=" * 100)

    # ------------------------------------------------------------
    # Load.
    # ------------------------------------------------------------

    df = load_dataset()

    audit = load_audit()

    normalization = extract_normalization(
        audit
    )

    network = extract_network(
        audit
    )

    print()
    print(
        f"Patients: {len(df)}"
    )

    print(
        "Architecture: 14 -> 7 -> 2"
    )

    print(
        "Hidden activation: tanh"
    )

    print(
        "Output activation: tanh"
    )

    # ------------------------------------------------------------
    # Forward pass.
    # ------------------------------------------------------------

    forward = forward_pass(
        df,
        normalization,
        network,
    )

    outputs = forward[
        "output_activation"
    ]

    cancer_score = outputs[
        :,
        0
    ]

    noncancer_score = outputs[
        :,
        1
    ]

    # ------------------------------------------------------------
    # PMML classification rule.
    #
    # Neuron 21 = Ovarian Cancer
    # Neuron 22 = Controls or Borderline
    #
    # Larger activation wins.
    # ------------------------------------------------------------

    predicted_class = (
        cancer_score
        > noncancer_score
    ).astype(
        int
    )

    y_true = df[
        "binary_label"
    ].to_numpy(
        dtype=int
    )

    # ------------------------------------------------------------
    # Metrics.
    # ------------------------------------------------------------

    confusion = compute_confusion(
        y_true,
        predicted_class,
    )

    auc = compute_auc(
        y_true,
        cancer_score,
    )

    # ------------------------------------------------------------
    # Exact equality check.
    #
    # Because the output neurons in the original network are
    # nearly opposite, this should behave as expected.
    # ------------------------------------------------------------

    prediction_ties = (
        cancer_score
        == noncancer_score
    )

    # ------------------------------------------------------------
    # Prediction dataframe.
    # ------------------------------------------------------------

    predictions = pd.DataFrame(
        {
            "sample_id": df[
                "sample_id"
            ],
            "diagnosis": df[
                "diagnosis"
            ],
            "true_binary_label": y_true,
            "cancer_output_neuron_21": cancer_score,
            "noncancer_output_neuron_22": noncancer_score,
            "output_margin_21_minus_22": (
                cancer_score
                - noncancer_score
            ),
            "predicted_binary_label": predicted_class,
            "predicted_class": np.where(
                predicted_class == 1,
                "Ovarian Cancer",
                "Controls or Borderline",
            ),
        }
    )

    predictions.to_csv(
        PREDICTIONS_PATH,
        index=False,
    )

    # ------------------------------------------------------------
    # Forward-pass trace.
    #
    # This creates a fully inspectable table containing every
    # normalized input, every hidden activation, and both outputs.
    # ------------------------------------------------------------

    trace = pd.DataFrame(
        {
            "sample_id": df[
                "sample_id"
            ],
        }
    )

    for j, feature in enumerate(
        FEATURES
    ):

        trace[
            f"raw_{feature}"
        ] = forward[
            "raw"
        ][
            :,
            j
        ]

        trace[
            f"normalized_{feature}"
        ] = forward[
            "normalized"
        ][
            :,
            j
        ]

    for i, neuron_id in enumerate(
        HIDDEN_IDS
    ):

        trace[
            f"hidden_{neuron_id}_linear"
        ] = forward[
            "hidden_linear"
        ][
            :,
            i
        ]

        trace[
            f"hidden_{neuron_id}_activation"
        ] = forward[
            "hidden_activation"
        ][
            :,
            i
        ]

    trace[
        "output_21_linear"
    ] = forward[
        "output_linear"
    ][
        :,
        0
    ]

    trace[
        "output_21_activation"
    ] = cancer_score

    trace[
        "output_22_linear"
    ] = forward[
        "output_linear"
    ][
        :,
        1
    ]

    trace[
        "output_22_activation"
    ] = noncancer_score

    trace[
        "predicted_binary_label"
    ] = predicted_class

    trace.to_csv(
        TRACE_PATH,
        index=False,
    )

    # ------------------------------------------------------------
    # Results.
    # ------------------------------------------------------------

    results = {

        "dataset": "GSE94533",

        "n_patients": int(
            len(df)
        ),

        "historical_model": {
            "model_name": (
                "Dataset complete with classificat_MLP 14-7-2"
            ),
            "architecture": "14-7-2",
            "input_features": 14,
            "hidden_neurons": 7,
            "output_neurons": 2,
            "hidden_activation": "tanh",
            "output_activation": "tanh",
            "cancer_output_neuron": "21",
            "noncancer_output_neuron": "22",
        },

        "input_preprocessing": {
            "formula": "normalized = raw * scale + shift",
            "source": str(
                AUDIT_PATH
            ),
        },

        "classification_rule": (
            "Choose the output neuron with the larger activation."
        ),

        "continuous_score": (
            "activation of cancer output neuron 21"
        ),

        "auc": auc,

        "confusion_matrix_metrics": confusion,

        "number_of_prediction_ties": int(
            prediction_ties.sum()
        ),

        "output_range": {
            "cancer_min": float(
                cancer_score.min()
            ),
            "cancer_max": float(
                cancer_score.max()
            ),
            "noncancer_min": float(
                noncancer_score.min()
            ),
            "noncancer_max": float(
                noncancer_score.max()
            ),
        },

        "WARNING_about_performance": (
            "These metrics are resubstitution metrics on the same "
            "179-patient dataset used in the historical study. "
            "They must NOT be interpreted as independent clinical "
            "generalization performance."
        ),

        "clinical_status": (
            "Research reproduction only; not a clinical diagnostic tool."
        ),
    }

    RESULTS_PATH.write_text(
        json.dumps(
            results,
            indent=2,
        )
    )

    # ------------------------------------------------------------
    # Terminal report.
    # ------------------------------------------------------------

    print()
    print("=" * 100)
    print("FORWARD PASS COMPLETE")
    print("=" * 100)

    print()
    print(
        "Input normalization verified:"
    )

    print(
        "  PMML shift/scale agrees with LinearNorm endpoints."
    )

    print()
    print(
        "Network:"
    )

    print(
        "  14 -> 7 tanh -> 2 tanh"
    )

    print()
    print(
        "Output neuron 21 = Ovarian Cancer"
    )

    print(
        "Output neuron 22 = Controls or Borderline"
    )

    print()
    print(
        "ROC AUC on the 179-patient dataset:"
    )

    print(
        f"  {auc:.6f}"
        if auc is not None
        else "  unavailable"
    )

    print()
    print(
        "Classification metrics:"
    )

    print(
        f"  Accuracy:    {confusion['accuracy']:.6f}"
    )

    print(
        f"  Sensitivity: {confusion['sensitivity']:.6f}"
        if confusion["sensitivity"] is not None
        else "  Sensitivity: unavailable"
    )

    print(
        f"  Specificity: {confusion['specificity']:.6f}"
        if confusion["specificity"] is not None
        else "  Specificity: unavailable"
    )

    print()
    print(
        "Confusion matrix:"
    )

    print(
        f"  TP={confusion['tp']} "
        f"FN={confusion['fn']}"
    )

    print(
        f"  FP={confusion['fp']} "
        f"TN={confusion['tn']}"
    )

    print()
    print(
        f"Prediction ties: "
        f"{int(prediction_ties.sum())}"
    )

    print()
    print("=" * 100)
    print("OUTPUT FILES")
    print("=" * 100)

    print()
    print(
        f"Predictions:\n  {PREDICTIONS_PATH}"
    )

    print()
    print(
        f"Results:\n  {RESULTS_PATH}"
    )

    print()
    print(
        f"Forward-pass trace:\n  {TRACE_PATH}"
    )

    print()
    print("=" * 100)
    print("IMPORTANT INTERPRETATION")
    print("=" * 100)

    print(
        """
This is a reproduction of the published network's forward pass,
not a new trained model.

The reported metrics are NOT independent validation because this
dataset overlaps the historical development cohort.

The next scientifically useful experiment is therefore NOT:

    "train a bigger neural network"

Instead we should establish an honest modern benchmark using
patient-level cross-validation, then compare:

    historical 14-miRNA network
              vs
    regularized modern miRNA model
              vs
    multimodal genomic + clinical/symptom model

with completely separated validation data.
"""
    )


if __name__ == "__main__":
    main()
