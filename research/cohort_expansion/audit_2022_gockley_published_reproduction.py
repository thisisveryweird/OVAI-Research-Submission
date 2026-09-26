#!/usr/bin/env python3
"""Audit published Gockley split and PMML reproduction without model fitting."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
import json
from zipfile import ZipFile
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from audit_2022_gockley_workbook import (
    DIAGNOSIS_COLUMN,
    MODEL_ARCHIVE,
    MODEL_FILENAME,
    REPORT_DIR,
    WORKBOOK,
    build_records,
    load_workbook_data,
    parse_pmml_model,
    predict_pmml,
)

SPLITS = ["Training", "Test", "Validation"]
ENDPOINTS = ["published_stage_I_to_IIIA2", "strict_stage_I_II"]
MODEL_COMPARISON = {
    "signature-best-10.xml": {
        "source_description": "Supplement caption: neural network trained on early-stage cohort samples.",
        "interpretation": "Human early-stage model: human malignancy-versus-benign/control target and 14 historical signature inputs.",
    },
    "7miRTrimmedXenoNet-38.xml": {
        "source_description": "Supplement caption: empirically optimized neural network on sera from three ovarian cancer xenografts using seven miRNAs.",
        "interpretation": "Murine xenograft model: seven features, Group Cancer/Control target, and explicit softmax output normalization.",
    },
    "NewNetXenograft-1.xml": {
        "source_description": "Supplement caption: neural network trained on sera from three ovarian cancer xenografts using 14 miRNAs.",
        "interpretation": "Murine xenograft model: 14 features, Group Cancer/Control target, and logistic/exponential layers.",
    },
}


def cancer_stage(stage: object) -> str:
    if stage is None or pd.isna(stage):
        return ""
    value = str(stage).strip().upper()
    if value.startswith("1"):
        return "Stage I"
    if value.startswith("2"):
        return "Stage II"
    if value == "3A1":
        return "Stage IIIA1"
    if value == "3A2":
        return "Stage IIIA2"
    return f"Other {value}"


def split_counts(records: list[dict[str, object]]) -> pd.DataFrame:
    rows = []
    partitions = [("All", records)] + [(split, [r for r in records if r["set"] == split]) for split in SPLITS]
    for split, subset in partitions:
        labels = Counter(r["diagnosis"] for r in subset)
        stages = Counter(cancer_stage(r["stage_raw"]) for r in subset if r["diagnosis"] == "Cancer Early")
        cancer_n = labels.get("Cancer Early", 0)
        benign_n = labels.get("Benign", 0)
        healthy_n = labels.get("Control", 0)
        stage_i = stages.get("Stage I", 0)
        stage_ii = stages.get("Stage II", 0)
        stage_iiia1 = stages.get("Stage IIIA1", 0)
        stage_iiia2 = stages.get("Stage IIIA2", 0)
        rows.append(
            {
                "split": split,
                "total_n": len(subset),
                "cancer_n": cancer_n,
                "benign_n": benign_n,
                "healthy_n": healthy_n,
                "stage_I_n": stage_i,
                "stage_II_n": stage_ii,
                "stage_IIIA1_n": stage_iiia1,
                "stage_IIIA2_n": stage_iiia2,
                "published_endpoint_n": cancer_n + benign_n + healthy_n,
                "strict_stage_I_II_cases_n": stage_i + stage_ii,
                "stage_accounts_for_cancer_rows": stage_i + stage_ii + stage_iiia1 + stage_iiia2 == cancer_n,
            }
        )
    return pd.DataFrame(rows)


def evaluate_splits(records: list[dict[str, object]], predictions: list[dict[str, object]]) -> pd.DataFrame:
    rows = []
    for endpoint in ENDPOINTS:
        for split in ["All"] + SPLITS:
            selected = [
                i for i, record in enumerate(records)
                if split == "All" or record["set"] == split
            ]
            chosen = []
            for i in selected:
                diagnosis = records[i]["diagnosis"]
                stage_group = records[i]["stage_group"]
                is_control = diagnosis in {"Benign", "Control"}
                is_cancer = diagnosis == "Cancer Early"
                if endpoint == "published_stage_I_to_IIIA2":
                    include = is_control or is_cancer
                    positive = is_cancer
                else:
                    positive = is_cancer and stage_group == "Stage I/II"
                    include = is_control or positive
                if include:
                    chosen.append((i, positive))
            index = [item[0] for item in chosen]
            truth = np.array([item[1] for item in chosen], dtype=bool)
            scores = np.array([float(predictions[i]["output_malignant"]) for i in index])
            margins = np.array([float(predictions[i]["score_margin_malignant"]) for i in index])
            predicted = np.array([predictions[i]["predicted_class"] == "malignant" for i in index], dtype=bool)
            tp = int(np.sum(truth & predicted))
            tn = int(np.sum(~truth & ~predicted))
            fp = int(np.sum(~truth & predicted))
            fn = int(np.sum(truth & ~predicted))
            rows.append(
                {
                    "endpoint": endpoint,
                    "split": split,
                    "n_cases": int(truth.sum()),
                    "n_controls": int((~truth).sum()),
                    "n_total": len(truth),
                    "roc_auc_malignant_output": float(roc_auc_score(truth, scores)),
                    "roc_auc_output_margin": float(roc_auc_score(truth, margins)),
                    "roc_auc_negative_benign_output_diagnostic": float(
                        roc_auc_score(truth, [-float(predictions[i]["output_benign_control"]) for i in index])
                    ),
                    "sensitivity_pmml_argmax": tp / (tp + fn) if tp + fn else float("nan"),
                    "specificity_pmml_argmax": tn / (tn + fp) if tn + fp else float("nan"),
                    "tn": tn,
                    "fp": fp,
                    "fn": fn,
                    "tp": tp,
                    "classification_rule": "maximum mapped output neuron; equivalent to malignant-minus-benign/control margin > 0 (no tie observed)",
                    "probability_cutoff": "not encoded in PMML",
                }
            )
    return pd.DataFrame(rows)


def model_metadata(archive: ZipFile, filename: str) -> dict[str, object]:
    xml_bytes = archive.read(filename)
    root = ET.fromstring(xml_bytes)
    network = root.find("NeuralNetwork")
    if network is None:
        raise ValueError(f"NeuralNetwork not found in {filename}")
    inputs = [
        item.find(".//NormContinuous").get("field")
        for item in network.findall("./NeuralInputs/NeuralInput")
    ]
    layers = network.findall("NeuralLayer")
    outputs = []
    for item in network.findall("./NeuralOutputs/NeuralOutput"):
        mapping = item.find(".//NormDiscrete")
        outputs.append(
            {
                "output_neuron": item.get("outputNeuron"),
                "field": mapping.get("field") if mapping is not None else None,
                "class": mapping.get("value") if mapping is not None else None,
            }
        )
    mining_target = next(
        (field.get("name") for field in network.findall("./MiningSchema/MiningField") if field.get("usageType") == "predicted"),
        None,
    )
    last_layer = layers[-1] if layers else None
    application = root.find("./Header/Application")
    input_norms = [item.find(".//NormContinuous") for item in network.findall("./NeuralInputs/NeuralInput")]
    return {
        "filename": filename,
        "model_name": network.get("modelName"),
        "function_name": network.get("functionName"),
        "inputs": inputs,
        "architecture": [len(layer.findall("Neuron")) for layer in layers],
        "activations": [layer.get("activationFunction", "identity") for layer in layers],
        "output_normalization": last_layer.get("normalizationMethod") if last_layer is not None else None,
        "target": mining_target,
        "classes": outputs,
        "normalization_transforms": len(network.findall("./NeuralInputs/NeuralInput")),
        "complete_input_norms": all(
            norm is not None and len(norm.findall("LinearNorm")) == 2
            for norm in input_norms
        ),
        "application_name": application.get("name") if application is not None else None,
        "application_version": application.get("version") if application is not None else None,
        "annotation_count": len(root.findall(".//Annotation")),
        "xml_comments_present": b"<!--" in xml_bytes,
        "has_probability_output_field": any(e.tag == "OutputField" and e.get("feature") == "probability" for e in root.iter()),
        "threshold_attributes": [
            {"element": element.tag, "attributes": element.attrib}
            for element in root.iter()
            if any("threshold" in key.lower() for key in element.attrib)
        ],
    }


def model_summary_rows(models: list[dict[str, object]]) -> list[str]:
    output = []
    for model in models:
        architecture = "-".join([str(len(model["inputs"]))] + [str(n) for n in model["architecture"]])
        classes = ", ".join(f"{item['output_neuron']}={item['class']}" for item in model["classes"])
        source = MODEL_COMPARISON[model["filename"]]
        output.append(
            f"| `{model['filename']}` | `{model['model_name']}` | {len(model['inputs'])} ({architecture}) | "
            f"{', '.join(model['activations'])}; output normalization `{model['output_normalization']}`; "
            f"{model['normalization_transforms']} input transforms, complete={model['complete_input_norms']} | "
            f"`{model['target']}`; {classes} | {source['interpretation']} |"
        )
    return output


def model_input_details(models: list[dict[str, object]]) -> list[str]:
    details = []
    for model in models:
        inputs = ", ".join(f"`{feature}`" for feature in model["inputs"])
        source = MODEL_COMPARISON[model["filename"]]["source_description"]
        details.append(
            f"- `{model['filename']}` input order: {inputs}. Each input has a `NormContinuous` "
            f"transform with two `LinearNorm` points: {model['normalization_transforms']} transforms; "
            f"all complete={model['complete_input_norms']}. PMML producer: "
            f"{model['application_name']} {model['application_version']}; annotations={model['annotation_count']}; "
            f"XML comments present={model['xml_comments_present']}. {source}"
        )
    return details


def fmt(value: object) -> str:
    return "NA" if pd.isna(value) else f"{float(value):.4f}"


def main() -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    headers, _, _, _ = load_workbook_data()
    if DIAGNOSIS_COLUMN not in headers:
        raise ValueError("Expected human diagnosis field not found")
    _, raw_rows, features, _ = load_workbook_data()
    records = build_records(headers, raw_rows, features)
    ids = [record["sample_id"] for record in records]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate patient IDs detected")
    if {record["set"] for record in records} != set(SPLITS):
        raise ValueError("Workbook split labels do not match Training/Test/Validation")
    split_table = split_counts(records)
    expected_sizes = {"Training": 165, "Test": 55, "Validation": 55}
    actual_sizes = dict(zip(split_table["split"], split_table["total_n"]))
    if {name: actual_sizes[name] for name in SPLITS} != expected_sizes:
        raise ValueError(f"Unexpected 3:1:1 split sizes: {actual_sizes}")
    if not split_table["stage_accounts_for_cancer_rows"].all():
        raise ValueError("Cancer stage rows do not reconcile with diagnosis counts")
    if not split_table["published_endpoint_n"].eq(split_table["total_n"]).all():
        raise ValueError("A workbook split contains a diagnosis outside the published endpoint")
    with ZipFile(MODEL_ARCHIVE) as archive:
        archive_names = archive.namelist()
        models = [model_metadata(archive, name) for name in archive_names]
        parsed_model = parse_pmml_model(archive, MODEL_FILENAME)
    if set(parsed_model["output_mapping"].values()) != {"malignant", "benign/control"}:
        raise ValueError("Unexpected class mapping in human PMML")
    if any(feature not in features for feature in [item["feature"] for item in parsed_model["inputs"]]):
        raise ValueError("Human PMML inputs are not all present in Data Sheet 2")
    predictions = [predict_pmml(parsed_model, record) for record in records]
    metrics = evaluate_splits(records, predictions)
    split_table.to_csv(REPORT_DIR / "2022_gockley_published_split_stage_counts.csv", index=False)
    metrics.to_csv(REPORT_DIR / "2022_gockley_published_reproduction_metrics.csv", index=False)
    model_lines = model_summary_rows(models)
    input_details = model_input_details(models)
    validation = metrics.loc[
        (metrics["endpoint"] == "published_stage_I_to_IIIA2") & (metrics["split"] == "Validation")
    ].iloc[0]
    validation_strict = metrics.loc[
        (metrics["endpoint"] == "strict_stage_I_II") & (metrics["split"] == "Validation")
    ].iloc[0]
    stage_table = [
        "| Split | n | Cancer | Benign | Healthy | Stage I | Stage II | IIIA1 | IIIA2 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in split_table.to_dict(orient="records"):
        stage_table.append(
            f"| {row['split']} | {row['total_n']} | {row['cancer_n']} | {row['benign_n']} | {row['healthy_n']} | "
            f"{row['stage_I_n']} | {row['stage_II_n']} | {row['stage_IIIA1_n']} | {row['stage_IIIA2_n']} |"
        )
    published_endpoint_table = [
        "| Endpoint | Split | Cases | Controls | n | AUC cancer output | AUC output margin | Sensitivity (argmax) | Specificity (argmax) | TN/FP/FN/TP |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in metrics.loc[metrics["split"].isin(["All", "Training", "Test", "Validation"])].to_dict(orient="records"):
        published_endpoint_table.append(
            f"| {row['endpoint']} | {row['split']} | {row['n_cases']} | {row['n_controls']} | {row['n_total']} | "
            f"{fmt(row['roc_auc_malignant_output'])} | {fmt(row['roc_auc_output_margin'])} | "
            f"{fmt(row['sensitivity_pmml_argmax'])} | {fmt(row['specificity_pmml_argmax'])} | "
            f"{row['tn']}/{row['fp']}/{row['fn']}/{row['tp']} |"
        )
    report = [
        "# Gockley Published-Experiment Reproduction Audit",
        "",
        "Scope: compare the actual workbook assignments and frozen PMML outputs with the article’s human validation experiment. No fitting, refitting, normalization, threshold selection, or model modification was performed.",
        "",
        "## Endpoint Definitions",
        "",
        "- Published-style: all 75 workbook cancer cases (FIGO I-IIIA2) versus 100 benign adnexal masses plus 100 healthy controls (n=275). The workbook has 2 Stage IIIA1 and 1 Stage IIIA2 case; all 75 belong in this endpoint.",
        "- OVAI-stricter: only Stage I/II cancer versus benign plus healthy controls (72 vs 200; n=272). The three Stage IIIA cases are excluded, not relabeled.",
        "",
        "## Workbook Split Check",
        "",
        "The workbook has exactly three labels, no duplicate patient IDs, and every row belongs to one set. Counts are 165 Training, 55 Test, and 55 Validation: an exact 3:1:1 allocation of 275 rows. The diagnosis and stage counts reconcile within each split. This verifies aggregate allocation against the paper; individual sample-to-split assignments cannot be cross-checked against a separate published roster because none is supplied in the article text inspected here.",
        *stage_table,
        "",
        "## Which XML Is the Human Model?",
        "",
        "| Archive member | PMML model name | Inputs (architecture) | Activations / output normalization; input transforms | Target and mapped output neurons | Interpretation |",
        "|---|---|---:|---|---|---|",
        *model_lines,
        "",
        *input_details,
        "",
        "The archive caption identifies `signature-best-10.xml` as trained on early-stage cohort samples. Its target is `Malignant vs benign/control`, and its 14 inputs exactly match the historical human 14-miRNA signature. The other two captions explicitly identify xenograft models; their target is `Group` with Cancer/Control categories. `7miRTrimmedXenoNet-38.xml` has seven inputs and softmax output normalization; `NewNetXenograft-1.xml` has 14 inputs but an animal-model target and logistic/exponential layers. Therefore `signature-best-10.xml` is the supplied human Figure 3 candidate; its filename says best-10 but its model declares 14 inputs.",
        "",
        "## Probability and Threshold Semantics",
        "",
        "`signature-best-10.xml` has a tanh output layer with two mapped class neurons (22=benign/control, 23=malignant). It has no output `normalizationMethod` (`softmax`/`simplemax`), no probability output field, no calibration transform, and no numeric threshold. Its outputs are bounded neural-network scores, not probabilities. The exact cancer score is the raw activation of the mapped malignant neuron; the binary PMML class rule selects the larger output, equivalently malignant-minus-benign/control margin > 0. That argmax boundary is not the paper’s 50% predicted-probability cutoff.",
        "- The paper’s 50% probability threshold is not recoverable from this PMML alone. Missing information is the mathematical mapping from the tanh class scores to calibrated cancer probability (or an explicitly exported probability output/normalization rule). The class mapping and weights do not specify this conversion. No guessed tanh-to-probability transform was applied.",
        "",
        "## Frozen PMML Results",
        "",
        "ROC-AUC is calculated first from the explicit malignant output neuron; the output margin is also shown because the two class neurons are not probability-normalized. Sensitivity/specificity and confusion matrices use only the PMML maximum-output-neuron rule; they are not presented as results at the paper’s 50% probability cutoff.",
        *published_endpoint_table,
        "",
        f"Validation, published-style endpoint: n={int(validation['n_total'])} ({int(validation['n_cases'])} cancer, {int(validation['n_controls'])} controls); AUC={validation['roc_auc_malignant_output']:.4f} from the malignant output neuron (margin AUC={validation['roc_auc_output_margin']:.4f}). PMML argmax gives sensitivity={validation['sensitivity_pmml_argmax']:.4f}, specificity={validation['specificity_pmml_argmax']:.4f}, confusion TN/FP/FN/TP={validation['tn']}/{validation['fp']}/{validation['fn']}/{validation['tp']}. The stricter endpoint has the same validation n and scores because this split contains no Stage IIIA cases (AUC={validation_strict['roc_auc_malignant_output']:.4f}).",
        "- Published Figure 3: AUC 0.87, sensitivity 73%, specificity 91% at predicted probability 50%. Reconstructed AUC 0.8867 is close but not identical (difference +0.0167); the PMML argmax operating point is 66.7%/95.0%, not 73%/91%. The published cutoff cannot be tested without the missing probability conversion.",
        "- Integer-count check: validation has 15 cancer and 40 non-cancer rows. 73% sensitivity is consistent with 11/15 (73.3%); 91% specificity cannot be obtained from 40 controls as an integer confusion count (36/40=90.0%; 37/40=92.5%). Thus the published operating percentages do not exactly correspond to an integer confusion matrix on the workbook Validation rows as counted. The article passage does not separately state another denominator for those percentages, so the discrepancy remains unresolved from the files.",
        "- Diagnostic only, not an accepted substitute: using the inverse benign-node score gives validation AUC 0.8667, which rounds to 0.87, but the XML maps the malignant class to neuron 23 and does not state that the inverse benign score was used. This coincidence is insufficient to claim reproduction.",
        "",
        "## Model-Selection Independence",
        "",
        "The paper says human samples were randomized 3:1:1 and that more than 100,000 networks were built, retaining the best based on Test-set class assignment. Under that stated procedure, Training is used for fitting; Test is used for model selection and is therefore not an unbiased final evaluation after selecting among 100,000 candidates; Validation is the only set described as reserved for final evaluation. It is an internal held-out split from the same 275-subject cohort, not an external cohort. This independence interpretation assumes Validation was not consulted for model selection or threshold calibration; the paper reports it as the validation set, while the PMML itself has no row-level training lineage.",
        "",
        "No model was fit, refit, tuned, recalibrated, or altered. Original workbook and ZIP were read in place. Detailed split counts are in `2022_gockley_published_split_stage_counts.csv`; fixed-score metrics are in `2022_gockley_published_reproduction_metrics.csv`.",
        "",
    ]
    path = REPORT_DIR / "2022_gockley_published_reproduction_audit.md"
    path.write_text("\n".join(report), encoding="utf-8")
    print(f"Wrote {path.relative_to(Path.cwd())}")
    print(f"Validation published-style: n={int(validation['n_total'])}, cases={int(validation['n_cases'])}, controls={int(validation['n_controls'])}, AUC={validation['roc_auc_malignant_output']:.4f}, PMML argmax sens={validation['sensitivity_pmml_argmax']:.4f}, spec={validation['specificity_pmml_argmax']:.4f}; confusion={validation['tn']}/{validation['fp']}/{validation['fn']}/{validation['tp']}.")
    print("No probability cutoff was applied.")


if __name__ == "__main__":
    main()
