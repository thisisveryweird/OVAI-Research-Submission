#!/usr/bin/env python3
"""Analyze saved GSE94533 out-of-fold errors without refitting models."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)


PROJECT = Path(__file__).resolve().parents[1]
OOF = PROJECT / "reports" / "tables" / "gse94533_enhancement_oof_predictions.csv"
META = PROJECT / "data" / "processed" / "gse94533_metadata.csv"
CLINICAL = PROJECT / "data" / "processed" / "gse94533_14mirna_clinical.csv.gz"
OUT_PREDICTIONS = PROJECT / "reports" / "error_analysis_predictions.csv"
OUT_SUBGROUPS = PROJECT / "reports" / "error_analysis_subgroups.csv"
OUT_REPORT = PROJECT / "reports" / "error_analysis.md"

MIRNAS = [
    "hsa-miR-23b-3p", "hsa-miR-29a-3p", "hsa-miR-32-5p",
    "hsa-miR-92a-3p", "hsa-miR-150-5p", "hsa-miR-200a-3p",
    "hsa-miR-200c-3p", "hsa-miR-203a", "hsa-miR-320c",
    "hsa-miR-320d", "hsa-miR-335-5p", "hsa-miR-450b-5p",
    "hsa-miR-1246", "hsa-miR-1307-5p",
]
CLINICAL_VARS = ["age", "ca125"]
MODEL_NAMES = [
    "median_logistic",
    "knn_logistic",
    "calibrated_svm",
    "random_forest",
    "xgboost",
    "ratio_logistic",
    "stacked_ensemble",
]


def binary_metrics(frame: pd.DataFrame, probability: str, subset: pd.Series) -> dict[str, object]:
    data = frame.loc[subset].copy()
    if len(data) == 0 or data["binary_label"].nunique() < 2:
        return {"n": int(len(data)), "status": "insufficient_class_variation"}
    pred = (data[probability] >= 0.5).astype(int)
    tn, fp, fn, tp = confusion_matrix(data["binary_label"], pred, labels=[0, 1]).ravel()
    return {
        "n": int(len(data)),
        "positive": int(data["binary_label"].sum()),
        "negative": int((data["binary_label"] == 0).sum()),
        "roc_auc": float(roc_auc_score(data["binary_label"], data[probability])),
        "pr_auc": float(average_precision_score(data["binary_label"], data[probability])),
        "accuracy": float(accuracy_score(data["binary_label"], pred)),
        "balanced_accuracy": float(balanced_accuracy_score(data["binary_label"], pred)),
        "sensitivity": float(recall_score(data["binary_label"], pred, zero_division=0)),
        "specificity": float(tn / (tn + fp)) if tn + fp else np.nan,
        "ppv": float(precision_score(data["binary_label"], pred, zero_division=0)),
        "f1": float(f1_score(data["binary_label"], pred, zero_division=0)),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
        "status": "completed",
    }


def subgroup_row(model: str, analysis: str, subgroup: str, frame: pd.DataFrame, probability: str, mask: pd.Series) -> dict[str, object]:
    return {
        "model": model,
        "analysis": analysis,
        "subgroup": subgroup,
        **binary_metrics(frame, probability, mask),
    }


def format_value(value: object) -> str:
    if pd.isna(value):
        return "NA"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def feature_contrasts(frame: pd.DataFrame, probability: str, positive_class: bool) -> list[str]:
    if positive_class:
        errors = frame[(frame["binary_label"] == 1) & (frame[probability] < 0.5)]
        reference = frame[(frame["binary_label"] == 1) & (frame[probability] >= 0.5)]
    else:
        errors = frame[(frame["binary_label"] == 0) & (frame[probability] >= 0.5)]
        reference = frame[(frame["binary_label"] == 0) & (frame[probability] < 0.5)]
    contrasts = []
    for feature in MIRNAS:
        error_values = np.log1p(errors[feature].clip(lower=0)).dropna()
        reference_values = np.log1p(reference[feature].clip(lower=0)).dropna()
        if len(error_values) < 2 or len(reference_values) < 2:
            continue
        pooled = np.sqrt((error_values.var(ddof=1) + reference_values.var(ddof=1)) / 2)
        delta = (error_values.mean() - reference_values.mean()) / pooled if pooled else 0.0
        contrasts.append((abs(delta), feature, delta, error_values.mean(), reference_values.mean()))
    contrasts.sort(reverse=True)
    return [
        f"{feature} delta={delta:.2f} (error mean={error_mean:.2f}, reference mean={reference_mean:.2f})"
        for _, feature, delta, error_mean, reference_mean in contrasts[:5]
    ]


def main() -> None:
    oof = pd.read_csv(OOF)
    metadata = pd.read_csv(META)
    clinical = pd.read_csv(CLINICAL)
    frame = (
        oof.merge(
            metadata[["sample_id_normalized", "geo_accession", "title", "stage", "histology", "diagnosis", "binary_label"]],
            on="sample_id_normalized",
            how="left",
            validate="one_to_one",
            suffixes=("", "_metadata"),
        )
        .drop(columns=["binary_label_metadata"], errors="ignore")
        .merge(
            clinical[["sample_id_normalized"] + MIRNAS + CLINICAL_VARS],
            on="sample_id_normalized",
            how="left",
            validate="one_to_one",
        )
    )
    if len(frame) != len(oof) or frame["sample_id_normalized"].duplicated().any():
        raise RuntimeError("OOF-to-metadata join is not one-to-one")

    frame["diagnosis"] = frame["diagnosis"].astype(str).str.strip()
    frame["stage_group"] = frame["stage"].fillna("Unknown").replace("", "Unknown")
    frame["sample_type"] = np.where(
        frame["sample_id_normalized"].str.startswith("e"),
        "e-sample",
        "numeric-sample",
    )
    frame["clinical_missingness"] = frame[CLINICAL_VARS].isna().sum(axis=1).map({0: "none", 1: "one_of_two", 2: "both"})
    frame["available_clinical_variables"] = frame.apply(
        lambda row: ",".join(
            variable for variable in CLINICAL_VARS if pd.notna(row[variable])
        ),
        axis=1,
    )
    frame["age_group"] = pd.cut(
        frame["age"],
        bins=[-np.inf, 50, 60, 70, np.inf],
        labels=["<=50", "51-60", "61-70", ">70"],
    ).astype("string").fillna("missing")

    prediction_rows = []
    for model in MODEL_NAMES:
        probability = f"{model}_probability"
        prediction = f"{model}_prediction"
        if probability not in frame or prediction not in frame:
            raise RuntimeError(f"Missing saved OOF columns for {model}")
        model_frame = frame.copy()
        model_frame["model"] = model
        model_frame["dataset"] = "GSE94533"
        model_frame["predicted_probability"] = model_frame[probability]
        model_frame["predicted_class"] = model_frame[prediction].astype(int)
        model_frame["true_label"] = model_frame["binary_label"].astype(int)
        model_frame["correct"] = model_frame["predicted_class"] == model_frame["binary_label"]
        model_frame["confidence"] = np.maximum(model_frame["predicted_probability"], 1 - model_frame["predicted_probability"])
        model_frame["signed_margin"] = 2 * model_frame["predicted_probability"] - 1
        model_frame["absolute_margin"] = model_frame["signed_margin"].abs()
        model_frame["error_type"] = np.select(
            [
                (model_frame["binary_label"] == 1) & (model_frame["predicted_class"] == 0),
                (model_frame["binary_label"] == 0) & (model_frame["predicted_class"] == 1),
            ],
            ["false_negative", "false_positive"],
            default="correct",
        )
        keep = [
            "sample_id_normalized", "geo_accession", "title", "dataset", "model",
            "diagnosis", "histology", "stage", "stage_group", "binary_label", "true_label",
            "age", "ca125", "clinical_missingness", "available_clinical_variables", "sample_type",
            "predicted_probability", "predicted_class", "correct", "error_type",
            "confidence", "signed_margin", "absolute_margin",
        ] + MIRNAS
        prediction_rows.append(model_frame[keep])
    predictions = pd.concat(prediction_rows, ignore_index=True)
    OUT_PREDICTIONS.parent.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(OUT_PREDICTIONS, index=False)

    subgroup_rows = []
    for model in MODEL_NAMES:
        probability = f"{model}_probability"
        subgroup_rows.append(subgroup_row(model, "invasive_vs_non_cancer", "all", frame, probability, pd.Series(True, index=frame.index)))
        subgroup_rows.append(subgroup_row(model, "invasive_vs_benign", "cancer_vs_benign", frame, probability, frame["diagnosis"].isin(["Cancer", "Benign"])))
        subgroup_rows.append(subgroup_row(model, "invasive_vs_healthy_control", "cancer_vs_control", frame, probability, frame["diagnosis"].isin(["Cancer", "Control"])))
        subgroup_rows.append(subgroup_row(model, "invasive_vs_borderline_non_invasive", "cancer_vs_borderline", frame, probability, frame["diagnosis"].isin(["Cancer", "Borderline"])))
        subgroup_rows.append(subgroup_row(model, "three_class_ranking", "cancer_vs_benign_vs_borderline", frame, probability, frame["diagnosis"].isin(["Cancer", "Benign", "Borderline"])))
        for group in ["Cancer", "Benign", "Borderline", "Control"]:
            subgroup_rows.append(subgroup_row(model, "diagnosis_subgroup", group, frame, probability, frame["diagnosis"].eq(group)))
        for stage in ["Stage I/II", "Stage III/IV"]:
            subgroup_rows.append(subgroup_row(model, "stage_vs_non_cancer", stage, frame, probability, frame["diagnosis"].eq("Cancer") & frame["stage"].eq(stage) | frame["diagnosis"].isin(["Benign", "Borderline", "Control"])))
        subgroup_rows.append(subgroup_row(model, "stage_I_vs_non_cancer", "Stage I only", frame, probability, frame["diagnosis"].eq("Cancer") & frame["stage"].astype(str).str.fullmatch("Stage I")))
        for missingness in ["none", "one_of_two", "both"]:
            subgroup_rows.append(subgroup_row(model, "clinical_missingness", missingness, frame, probability, frame["clinical_missingness"].eq(missingness)))
        for age_group in ["<=50", "51-60", "61-70", ">70", "missing"]:
            subgroup_rows.append(subgroup_row(model, "age_group", age_group, frame, probability, frame["age_group"].eq(age_group)))
        for sample_type in ["e-sample", "numeric-sample"]:
            subgroup_rows.append(subgroup_row(model, "sample_type_proxy", sample_type, frame, probability, frame["sample_type"].eq(sample_type)))
    subgroup = pd.DataFrame(subgroup_rows)
    subgroup.to_csv(OUT_SUBGROUPS, index=False)

    report_lines = [
        "# GSE94533 OOF Error Analysis",
        "",
        "This analysis uses the saved out-of-fold predictions in `reports/tables/gse94533_enhancement_oof_predictions.csv`. No model was refit and no threshold was selected from this analysis.",
        "",
        "## Cohort and Limitations",
        "",
        f"The joined cohort contains {len(frame)} unique GSE94533 patients: {int((frame['binary_label'] == 1).sum())} invasive cancer and {int((frame['binary_label'] == 0).sum())} non-cancer. Metadata include diagnosis, histology, stage group, age, and CA-125 availability. No batch/platform field is present in the verified metadata; the `e-sample` versus `numeric-sample` field is only a sample-ID provenance proxy, not a validated batch variable.",
        "",
        "The saved models are binary cancer-versus-non-cancer classifiers. Borderline tumors are not relabeled as cancer. The three-class section therefore reports one-vs-rest ranking and cancer-vs-benign-vs-borderline discrimination; it does not claim a three-class hard-label classifier or invent benign/borderline predictions.",
        "",
        "## Overall Error Pattern",
        "",
    ]
    summary = []
    for model in MODEL_NAMES:
        probability = f"{model}_probability"
        metrics = binary_metrics(frame, probability, pd.Series(True, index=frame.index))
        fn = int(((frame["binary_label"] == 1) & (frame[probability] < 0.5)).sum())
        fp = int(((frame["binary_label"] == 0) & (frame[probability] >= 0.5)).sum())
        summary.append((model, metrics, fn, fp))
        report_lines.append(
            f"- `{model}`: ROC-AUC {format_value(metrics.get('roc_auc'))}, balanced accuracy {format_value(metrics.get('balanced_accuracy'))}, false negatives {fn}, false positives {fp}."
        )
    report_lines.extend(["", "## False Negatives", ""])
    for model, _, _, _ in summary:
        p = f"{model}_probability"
        errors = frame[(frame["binary_label"] == 1) & (frame[p] < 0.5)]
        by_stage = errors["stage"].fillna("Unknown").value_counts().to_dict()
        by_histology = errors["histology"].fillna("Unknown").value_counts().head(8).to_dict()
        report_lines.append(f"### {model}")
        report_lines.append(f"- False negatives: {len(errors)}; median probability {format_value(errors[p].median())}; median age {format_value(errors['age'].median())}; missing clinical values in {int((errors['clinical_missingness'] != 'none').sum())}.")
        report_lines.append(f"- Stage distribution: {json.dumps({str(k): int(v) for k, v in by_stage.items()}, sort_keys=True)}.")
        report_lines.append(f"- Histology distribution: {json.dumps({str(k): int(v) for k, v in by_histology.items()}, sort_keys=True)}.")
        report_lines.append(f"- Largest descriptive miRNA contrasts versus correctly classified cancers (log1p scale, standardized mean difference): {'; '.join(feature_contrasts(frame, p, True)) or 'insufficient data'}.")
        report_lines.append("- Per-sample details are in `reports/error_analysis_predictions.csv`.")
    report_lines.extend(["", "## False Positives", ""])
    for model, _, _, _ in summary:
        p = f"{model}_probability"
        errors = frame[(frame["binary_label"] == 0) & (frame[p] >= 0.5)]
        by_diagnosis = errors["diagnosis"].value_counts().to_dict()
        report_lines.append(f"### {model}")
        report_lines.append(f"- False positives: {len(errors)}; median probability {format_value(errors[p].median())}; median age {format_value(errors['age'].median())}; diagnosis distribution: {json.dumps({str(k): int(v) for k, v in by_diagnosis.items()}, sort_keys=True)}.")
        report_lines.append(f"- Largest descriptive miRNA contrasts versus correctly classified non-cancer samples (log1p scale, standardized mean difference): {'; '.join(feature_contrasts(frame, p, False)) or 'insufficient data'}.")
    report_lines.extend(["", "## Borderline Tumors", ""])
    report_lines.append("Borderline tumors are a distinct non-invasive subgroup in the source labels. They are not silently merged into cancer.")
    for model in MODEL_NAMES:
        p = f"{model}_probability"
        for analysis, mask in [
            ("invasive cancer vs all non-cancer", frame["diagnosis"].isin(["Cancer", "Benign", "Borderline", "Control"])),
            ("invasive cancer vs benign", frame["diagnosis"].isin(["Cancer", "Benign"])),
            ("invasive cancer vs healthy controls", frame["diagnosis"].isin(["Cancer", "Control"])),
            ("invasive cancer vs borderline", frame["diagnosis"].isin(["Cancer", "Borderline"])),
        ]:
            metrics = binary_metrics(frame, p, mask)
            report_lines.append(f"- `{model}`, {analysis}: n={metrics.get('n')}, ROC-AUC={format_value(metrics.get('roc_auc'))}, balanced accuracy={format_value(metrics.get('balanced_accuracy'))}.")
        tri = frame[frame["diagnosis"].isin(["Cancer", "Benign", "Borderline"])].copy()
        tri_rows = []
        for class_name in ["Cancer", "Benign", "Borderline"]:
            tri_rows.append(f"{class_name} one-vs-rest AUC {format_value(roc_auc_score((tri['diagnosis'] == class_name).astype(int), tri[p]))}")
        report_lines.append(f"- `{model}`, three-class ranking (no hard-label 3-class predictions): " + "; ".join(tri_rows) + ".")
    report_lines.extend(["", "## Early-Stage Disease", ""])
    report_lines.append("Stage I-only counts are reported separately; the current verified metadata encode most early cases as `Stage I/II`, so Stage I-only analysis is limited.")
    for model in MODEL_NAMES:
        p = f"{model}_probability"
        for stage in ["Stage I/II", "Stage III/IV"]:
            mask = (frame["diagnosis"].eq("Cancer") & frame["stage"].eq(stage)) | frame["diagnosis"].isin(["Benign", "Borderline", "Control"])
            metrics = binary_metrics(frame, p, mask)
            report_lines.append(f"- `{model}`, {stage} vs non-cancer: n={metrics.get('n')}, positive={metrics.get('positive')}, ROC-AUC={format_value(metrics.get('roc_auc'))}, balanced accuracy={format_value(metrics.get('balanced_accuracy'))}, sensitivity={format_value(metrics.get('sensitivity'))}.")
        stage_i = frame["diagnosis"].eq("Cancer") & frame["stage"].astype(str).str.fullmatch("Stage I")
        mask = stage_i | frame["diagnosis"].isin(["Benign", "Borderline", "Control"])
        metrics = binary_metrics(frame, p, mask)
        report_lines.append(f"- `{model}`, Stage I only vs non-cancer: n={metrics.get('n')}, positive={metrics.get('positive')}, status={metrics.get('status')}, ROC-AUC={format_value(metrics.get('roc_auc'))}.")
    report_lines.extend(["", "## Missingness and Platform Proxies", ""])
    report_lines.append("Clinical missingness is reported as none, one of age/CA-125 missing, or both missing. These are descriptive subgroup analyses, not causal explanations.")
    report_lines.append("No validated platform or batch variable is available in the current GSE94533 metadata. Results split by e-sample versus numeric-sample ID are included only as a provenance proxy and should not be interpreted as a batch effect.")
    report_lines.extend(["", "## Outputs", "", "- Per-sample, per-model errors and margins: `reports/error_analysis_predictions.csv`.", "- Machine-readable subgroup metrics: `reports/error_analysis_subgroups.csv`.", "- Input: saved OOF predictions only; no refitting performed."])
    OUT_REPORT.write_text("\n".join(report_lines) + "\n")
    print(f"Saved {OUT_REPORT}")
    print(f"Saved {OUT_PREDICTIONS}")
    print(f"Saved {OUT_SUBGROUPS}")


if __name__ == "__main__":
    main()
