#!/usr/bin/env python3
"""Audit score semantics of the frozen Gockley PMML on Validation only."""

from __future__ import annotations

from pathlib import Path
import json
from zipfile import ZipFile
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd
from scipy.stats import kendalltau
from sklearn.metrics import roc_auc_score

from audit_2022_gockley_workbook import (
    MODEL_ARCHIVE,
    MODEL_FILENAME,
    REPORT_DIR,
    build_records,
    load_workbook_data,
    parse_pmml_model,
    predict_pmml,
)


EXPECTED_ROOT_CHILDREN = ["Header", "DataDictionary", "NeuralNetwork"]


def _rank_order(values: np.ndarray) -> np.ndarray:
    return np.argsort(np.argsort(values, kind="mergesort"), kind="mergesort")


def _score_row(
    name: str,
    formula: str,
    score: np.ndarray | None,
    base: np.ndarray,
    truth: np.ndarray,
    status: str,
    probability_note: str,
    monotonicity: str,
    evidence: str,
    hypothetical_50_boundary: str = "not applicable",
    not_evaluated_reason: str = "",
) -> dict[str, object]:
    if score is None:
        return {
            "representation": name,
            "formula": formula,
            "validation_cases": int(truth.sum()),
            "validation_controls": int(len(truth) - truth.sum()),
            "validation_auc": np.nan,
            "score_min": np.nan,
            "score_max": np.nan,
            "monotonicity_vs_malignant_output": monotonicity,
            "same_sample_order_as_malignant_output": "not evaluated",
            "kendall_tau_vs_malignant_output": np.nan,
            "probability_interpretation": probability_note,
            "pmml_evidence_classification": status,
            "hypothetical_50_boundary_not_applied": hypothetical_50_boundary,
            "not_evaluated_reason": not_evaluated_reason,
            "evidence": evidence,
        }
    order_same = bool(np.array_equal(_rank_order(score), _rank_order(base)))
    return {
        "representation": name,
        "formula": formula,
        "validation_cases": int(truth.sum()),
        "validation_controls": int(len(truth) - truth.sum()),
        "validation_auc": float(roc_auc_score(truth, score)),
        "score_min": float(np.min(score)),
        "score_max": float(np.max(score)),
        "monotonicity_vs_malignant_output": monotonicity,
        "same_sample_order_as_malignant_output": order_same,
        "kendall_tau_vs_malignant_output": float(kendalltau(base, score).statistic),
        "probability_interpretation": probability_note,
        "pmml_evidence_classification": status,
        "hypothetical_50_boundary_not_applied": hypothetical_50_boundary,
        "not_evaluated_reason": not_evaluated_reason,
        "evidence": evidence,
    }


def main() -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    headers, raw_rows, features, _ = load_workbook_data()
    records = build_records(headers, raw_rows, features)
    with ZipFile(MODEL_ARCHIVE) as archive:
        xml_bytes = archive.read(MODEL_FILENAME)
        root = ET.fromstring(xml_bytes)
        model = parse_pmml_model(archive, MODEL_FILENAME)
    root_child_names = [child.tag for child in list(root)]
    if root_child_names != EXPECTED_ROOT_CHILDREN:
        raise ValueError(f"Unexpected top-level PMML elements: {root_child_names}")
    network = root.find("NeuralNetwork")
    output_layers = network.findall("NeuralLayer")
    if len(output_layers) != 2:
        raise ValueError(f"Expected one hidden and one output layer; found {len(output_layers)}")
    output_layer = output_layers[-1]
    if output_layer.get("activationFunction") != "tanh":
        raise ValueError(f"Unexpected human output activation: {output_layer.attrib}")
    output_nodes = {node.get("id"): node for node in output_layer.findall("Neuron")}
    output_map = model["output_mapping"]
    class_node = {class_name: output_nodes[node_id] for node_id, class_name in output_map.items()}
    malignant_node = class_node["malignant"]
    benign_node = class_node["benign/control"]
    if len(malignant_node.findall("Con")) != 8 or len(benign_node.findall("Con")) != 8:
        raise ValueError("Output neurons do not have the expected eight hidden-layer inputs")
    if model["weight_count"] != 128 or model["bias_count"] != 10:
        raise ValueError("Unexpected number of serialized weights or biases")
    if any(layer.get("normalizationMethod") for layer in output_layers):
        raise ValueError("Unexpected output normalization method in this audit")
    if root.find(".//TransformationDictionary") is not None or root.find(".//LocalTransformations") is not None:
        raise ValueError("Unexpected PMML transformation outside the network")

    predictions = [predict_pmml(model, record) for record in records]
    selected = [
        index for index, record in enumerate(records)
        if record["set"] == "Validation" and record["diagnosis"] in {"Cancer Early", "Benign", "Control"}
    ]
    validation_records = [records[index] for index in selected]
    validation_predictions = [predictions[index] for index in selected]
    truth = np.array([record["diagnosis"] == "Cancer Early" for record in validation_records], dtype=bool)
    malignant = np.array([float(item["output_malignant"]) for item in validation_predictions])
    benign = np.array([float(item["output_benign_control"]) for item in validation_predictions])
    margin = malignant - benign
    negative_benign = -benign
    sigmoid_malignant = 1.0 / (1.0 + np.exp(-malignant))
    tanh_unit_malignant = (malignant + 1.0) / 2.0
    inverted_tanh_benign = 1.0 - (benign + 1.0) / 2.0
    shifted = np.column_stack([malignant, benign])
    shifted = shifted - shifted.max(axis=1, keepdims=True)
    exp_values = np.exp(shifted)
    softmax_malignant = exp_values[:, 0] / exp_values.sum(axis=1)

    sum_outputs = malignant + benign
    ratio_note = (
        "Not evaluated as an AUC score: the two tanh outputs may be negative, are not constrained to sum to a positive constant, "
        "and PMML does not declare simplemax. Dividing by their sum is not a valid probability transform for this model."
    )
    zero_denominator = np.isclose(sum_outputs, 0.0, atol=1e-12)
    ratio_diagnostic = np.full_like(malignant, np.nan)
    ratio_diagnostic[~zero_denominator] = malignant[~zero_denominator] / sum_outputs[~zero_denominator]
    finite_ratios = ratio_diagnostic[np.isfinite(ratio_diagnostic)]
    ratio_row = _score_row(
        "malignant_over_sum_of_outputs",
        "m / (m + b)",
        None,
        malignant,
        truth,
        "mathematically unjustified",
        "No: outputs are signed and have no positive-sum constraint; the quotient is not a probability definition.",
        "not evaluated",
        "The PMML has independent tanh class outputs, no simplemax, and no invariant that m+b is a positive normalizer.",
        not_evaluated_reason=ratio_note,
    )
    ratio_row["score_min"] = float(finite_ratios.min()) if len(finite_ratios) else np.nan
    ratio_row["score_max"] = float(finite_ratios.max()) if len(finite_ratios) else np.nan
    ratio_row["not_evaluated_reason"] = (
        f"AUC not calculated. As an invalidity diagnostic only, m/(m+b) is outside [0,1] for "
        f"{int(np.sum((ratio_diagnostic < 0) | (ratio_diagnostic > 1)))} of {len(ratio_diagnostic)} Validation rows; "
        f"the finite quotient range is {finite_ratios.min():.4f} to {finite_ratios.max():.4f}. PMML declares no simplemax or positive-sum constraint."
    )
    score_rows = [
        _score_row(
            "malignant_output",
            "m = tanh(z_malignant), the raw output of the neuron mapped to malignant",
            malignant,
            malignant,
            truth,
            "directly supported by PMML semantics",
            "No: categorical class score, not a probability field.",
            "identity; reference ordering",
            "NeuralOutput maps neuron 23 to malignant; output activation is tanh; no output normalization is declared.",
        ),
        _score_row(
            "negative_benign_control_output",
            "-b, where b is the raw output of neuron 22 mapped to benign/control",
            negative_benign,
            malignant,
            truth,
            "directly supported as an oriented class score",
            "No: reverse-coded benign score, not a probability field.",
            "reverse orientation of the PMML-mapped benign score; not a function of m alone",
            "The mapped benign/control output should decrease as evidence for the opposite class increases; PMML does not prescribe its negation as a probability.",
        ),
        _score_row(
            "malignant_minus_benign_margin",
            "m - b; positive exactly when the malignant output exceeds the benign/control output",
            margin,
            malignant,
            truth,
            "directly supported by PMML classification semantics",
            "No: signed decision score, not a probability field.",
            "strictly increasing in the PMML decision margin; not a function of m alone",
            "This is the exact pairwise score underlying maximum-output-neuron classification.",
        ),
        _score_row(
            "tanh_malignant_to_unit_interval",
            "(m + 1) / 2",
            tanh_unit_malignant,
            malignant,
            truth,
            "plausible but unsupported",
            "No for the released model: unit-interval rescaling is not calibration and the target coding is absent.",
            "strictly increasing in m; preserves malignant-output ranking",
            "Tanh has range (-1, 1), but PMML does not declare bipolar target coding or this affine probability conversion.",
            "m = 0",
        ),
        _score_row(
            "inverted_tanh_benign_to_unit_interval",
            "1 - (b + 1) / 2",
            inverted_tanh_benign,
            malignant,
            truth,
            "plausible but unsupported",
            "No for the released model: inverse rescaling is not calibrated cancer probability.",
            "strictly decreasing in b; ranking equals -b, not m",
            "Tanh range permits this unit-interval rescaling, but no output coding or inverse-benign probability rule is encoded.",
            "b = 0",
        ),
        _score_row(
            "sigmoid_of_malignant_output",
            "sigmoid(m) = 1 / (1 + exp(-m))",
            sigmoid_malignant,
            malignant,
            truth,
            "plausible but unsupported",
            "No for the released model: sigmoid(m) is monotonic pseudo-probability, not a specified calibration.",
            "strictly increasing in m; preserves malignant-output ranking",
            "The PMML score is a tanh activation, not an encoded logit; no sigmoid output link or calibration is present.",
            "m = 0",
        ),
        _score_row(
            "posthoc_softmax_of_two_tanh_outputs",
            "exp(m) / (exp(m) + exp(b))",
            softmax_malignant,
            malignant,
            truth,
            "plausible but unsupported",
            "Would be a normalized pair only if softmax were specified; it is not the released model output.",
            "strictly increasing in m - b; ranking equals the PMML margin, not necessarily m",
            "Softmax is a standard PMML NeuralLayer normalization option, but this output layer declares tanh and no normalizationMethod.",
            "m = b",
        ),
        ratio_row,
    ]
    transformations = pd.DataFrame(score_rows)
    transformations.to_csv(REPORT_DIR / "2022_gockley_score_transformations.csv", index=False)

    app = root.find("./Header/Application")
    outputs_mapped = [
        (output.get("outputNeuron"), output.find(".//NormDiscrete").get("value"))
        for output in network.findall("./NeuralOutputs/NeuralOutput")
    ]
    outputs_ordered = {class_name: node_id for node_id, class_name in outputs_mapped}
    output_sum_min = float(sum_outputs.min())
    output_sum_max = float(sum_outputs.max())
    probability_output_fields = [
        element.attrib for element in root.iter()
        if element.tag == "OutputField" and element.get("feature", "").lower() == "probability"
    ]
    outside_elements = [
        element.tag for element in root.iter()
        if element.tag in {"TransformationDictionary", "LocalTransformations", "Targets", "Output", "ModelVerification", "Extension", "Annotation"}
    ]
    output_weights = []
    for class_name in ["benign/control", "malignant"]:
        node = class_node[class_name]
        connections = ", ".join(f"{item.get('from')}:{item.get('weight')}" for item in node.findall("Con"))
        output_weights.append(
            f"- Neuron {node.get('id')} ({class_name}): activation `tanh`; bias `{node.get('bias')}`; "
            f"incoming weights from hidden neuron IDs 14-21: `{connections}`."
        )
    base_auc = float(roc_auc_score(truth, malignant))
    inverse_benign_auc = float(roc_auc_score(truth, negative_benign))
    margin_auc = float(roc_auc_score(truth, margin))
    report = [
        "# Gockley PMML Score-Semantics Audit",
        "",
        "Scope: determine which validation score representations are encoded or justified by the frozen `signature-best-10.xml`. Validation only (15 cancers, 40 non-cancer controls); no weights, transforms, thresholds, or model artifacts were fitted or changed.",
        "",
        "## Serialized Network Mathematics",
        "",
        f"- PMML version: `{root.get('version')}`; top-level children: {', '.join(f'`{tag}`' for tag in [child.tag for child in list(root)])}; producer: `{app.get('name')}` version `{app.get('version')}`.",
        "- Network: `Ovaries_MLP 14-8-2`, classification; 14 ordered inputs, 8 hidden neurons, 2 output neurons; 128 serialized connection weights and 10 biases. All weights and biases were decoded and used in the frozen computation; both layers use tanh. The 14 input-side `NormContinuous` elements each include shift/scale and two `LinearNorm` points mapping the stored assay range to 0..1.",
        f"- Hidden neuron IDs are 14-21; output layer neuron IDs are 22 and 23. Finite tanh outputs lie strictly in (-1,1) (closed range convention [-1,1]). The raw Validation outputs empirically range: malignant `{malignant.min():.6f}` to `{malignant.max():.6f}`; benign/control `{benign.min():.6f}` to `{benign.max():.6f}`.",
        *output_weights,
        f"- `NeuralOutputs` maps neuron {outputs_ordered['benign/control']} to `benign/control` and neuron {outputs_ordered['malignant']} to `malignant`. This is a categorical mapping (`NormDiscrete`), not a probability output declaration.",
        f"- No `TransformationDictionary`, `LocalTransformations`, `Targets`, `Output`, `ModelVerification`, `Extension`, or `Annotation` elements occur outside the network (observed such elements: {outside_elements}). No probability `OutputField` occurs. The only `DerivedField` transforms are the input `NormContinuous` mappings inside `NeuralInputs`.",
        f"- The outputs are separately activated tanh neuron scores. No `normalizationMethod` is set on the output layer, and no equation constrains their sum; in Validation, their observed sum varies from {output_sum_min:.6f} to {output_sum_max:.6f}, rather than being fixed at 1. Each output is negative on some rows; the quotient `m/(m+b)` falls outside [0,1] for {int(np.sum((ratio_diagnostic < 0) | (ratio_diagnostic > 1)))}/55 Validation rows (range {finite_ratios.min():.4f} to {finite_ratios.max():.4f}).",
        "- Neither output is intended by the serialized PMML to be a probability: the XML declares class scores and a maximum-output classification rule. No calibrated probability map, bipolar target coding, sigmoid link, simplemax, or softmax is serialized.",
        "",
        "## Validation Score Transformations",
        "",
        "AUC uses the same frozen Validation labels for every row. 'Same ranking' compares the ordering of all 55 samples to the PMML-mapped malignant output. Transformations marked unsupported are mathematical diagnostics only; they are not substituted for the published probability or cutoff.",
        "",
        "| Representation | Formula | AUC | Same ordering as malignant output? | Probability status | PMML evidence classification |",
        "|---|---|---:|---:|---|---|",
    ]
    for row in transformations.to_dict(orient="records"):
        auc_text = "not calculated" if pd.isna(row["validation_auc"]) else f"{row['validation_auc']:.4f}"
        report.append(
            f"| {row['representation']} | `{row['formula']}` | {auc_text} | "
            f"{row['same_sample_order_as_malignant_output']} | {row['probability_interpretation']} | {row['pmml_evidence_classification']} |"
        )
    report.extend(
        [
            "",
            f"The inverse benign-node result is AUC {inverse_benign_auc:.4f}, close to the published 0.87, but PMML does not designate `-output_benign` as a cancer probability. The designated malignant-node AUC is {base_auc:.4f}; the exact PMML argmax margin AUC is {margin_auc:.4f}. Similar AUC values do not identify which unsupported probability conversion, if any, the paper’s software used.",
            f"- `m/(m+b)` was not scored as an AUC representation: the quotient is outside [0,1] for {int(np.sum((ratio_diagnostic < 0) | (ratio_diagnostic > 1)))}/55 Validation rows (range {finite_ratios.min():.4f} to {finite_ratios.max():.4f}), and PMML declares no simplemax or positive-sum constraint. A post-hoc two-output softmax is included only as an unsupported diagnostic; it is not in this network's PMML.",
            "",
            "## The 50% Cutoff and Published Result",
            "",
            "The paper states a 50% predicted-probability cutoff, but the released PMML has no probability transformation. The probability scale is therefore not recoverable from the released model files. Candidate monotone maps (tanh affine rescaling, sigmoid, post-hoc softmax) produce valid-looking numbers but are not selected or calibrated by the XML; no one can be called the paper’s predicted probability based on these files.",
            "- The exact PMML class rule is `output_malignant > output_benign/control`, equivalently margin `m-b > 0`; it is not a 50% probability threshold. No candidate 50% boundary was applied to labels.",
            "- The supplement caption calls `signature-best-10.xml` a neural network trained on early-stage cohort samples; its human target and 14 historical-signature inputs match the Figure 3 task. The paper states that >100,000 networks were compared and the best retained by Test-set assignment, but neither the filename nor PMML records a winner rank, selection flag, or Figure 3 AUC. 'best-10' has no defined technical meaning in the released files and does not mean 10 inputs. The file is the supplied human-model candidate; exact identity as the final selected winner is unverified.",
            "- The PMML producer is STATISTICA Automated Neural Networks (SANN) 2.0, but no software-side probability processing is serialized. Processing outside the PMML is possible but not demonstrated by any supplied file; its formula cannot be recovered here.",
            f"- On Validation, the designated malignant-node AUC is {base_auc:.4f}, versus the published 0.87 (95% CI 0.84-0.94); the value is within the reported interval and approximately reproduces ranking. The inverse benign-node AUC {inverse_benign_auc:.4f} also rounds to 0.87, but closeness does not justify choosing it.",
            "- PMML argmax classification on Validation is 10/15 cancers detected (66.7%) and 38/40 controls correctly negative (95.0%), confusion TN/FP/FN/TP=38/2/5/10. This is not the published 73%/91% point.",
            "- The stated Validation denominator has 15 cancers and 40 controls. 73% is compatible with 11/15=73.3%; 91% specificity is not an attainable integer fraction of 40 controls (36/40=90%, 37/40=92.5%). The exact published specificity therefore does not match this validation denominator under a standard confusion-matrix calculation. The paper passage does not provide a different denominator or exported probability map.",
            "",
            "## Final Reproducibility Classification",
            "",
            "| Level | Classification | Basis |",
            "|---|---|---|",
            "| Architecture reproducibility | Reproduced for the released XML | All 14 input transforms, neuron counts, activations, weights, biases, and class mappings are serialized and decoded; exact identity as the selected final winner is unverified. |",
            "| Ranking reproducibility | Approximately reproduced | Designated malignant output and PMML margin each give Validation AUC 0.8867 versus reported 0.87; the reported interval 0.84-0.94 contains this value. |",
            "| Operating-point reproducibility | Not reproducible from released files | The 50% predicted-probability scale/conversion is absent; PMML argmax gives 66.7%/95.0%, not 73%/91%. |",
            "",
            "## Scientific Classification",
            "",
            "- Released human PMML model successfully reconstructed: yes; correspondence to the exact final network selected from >100,000 candidates is unverified.",
            "- Validation ranking consistent with paper: approximately; AUC 0.8867 is close to and within the reported CI, but does not exactly equal 0.87.",
            "- Published 50% probability scale recoverable: no.",
            "- Published operating point recoverable: no.",
            "- Missing: the model’s probability conversion/calibration definition, any selected-network identifier, and a validation-set confusion-count/denominator consistent with the reported 91% specificity. The released artifacts do not establish whether software-side probability processing occurred outside PMML.",
            "- Supported dataset labels: same-cohort model reconstruction; historical-signature replication; approximate internal validation ranking reproduction. Not independent external validation, and not current OVAI validation (0/5 exact current features).",
            "",
            "No training, refitting, tuning, recalibration, or PMML modification was performed. The raw workbook and ZIP were read in place.",
            "",
        ]
    )
    (REPORT_DIR / "2022_gockley_score_semantics_audit.md").write_text("\n".join(report), encoding="utf-8")
    print(f"Validation AUC malignant={base_auc:.4f}, inverse benign={inverse_benign_auc:.4f}, margin={margin_auc:.4f}; rows={len(validation_records)}.")
    print(f"Wrote {REPORT_DIR / '2022_gockley_score_semantics_audit.md'}")
    print(f"Wrote {REPORT_DIR / '2022_gockley_score_transformations.csv'}")


if __name__ == "__main__":
    main()
