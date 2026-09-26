#!/usr/bin/env python3
"""Characterize transportability using the existing saved OVAI outputs.

This script intentionally does not refit models or use external target labels for
selection. It reads the saved analysis artifacts already present in the repo and
summarizes the evidence that the source-model failure is primarily a platform/data
shift problem rather than a lack of model complexity.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
REPORTS = ROOT / "reports"
OUTDIR = Path(__file__).resolve().parent


def load_json(path: Path):
    return json.loads(path.read_text())


def summarize_experiments() -> pd.DataFrame:
    payload = load_json(REPORTS / "gse106817_phase5_results.json")
    rows = []
    for row in payload["experiments"]:
        if row["test_dataset"] == "GSE106817":
            ppos = (row["confusion_matrix"]["tp"] + row["confusion_matrix"]["fp"]) / row["n"]
            row_summary = {
                "train_dataset": row["train_dataset"],
                "test_dataset": row["test_dataset"],
                "normalization": row["normalization"],
                "roc_auc": float(row["roc_auc"]),
                "pr_auc": float(row["pr_auc"]),
                "balanced_accuracy": float(row["balanced_accuracy"]),
                "sensitivity": float(row["sensitivity"]),
                "specificity": float(row["specificity"]),
                "threshold": float(row["threshold"]),
                "positive_rate": float(ppos),
                "accuracy": float(row["accuracy"]),
                "brier_score": float(row["brier_score"]),
            }
            rows.append(row_summary)
    return pd.DataFrame(rows).sort_values(["train_dataset", "normalization"]).reset_index(drop=True)


def save_distribution_shift_plot(distribution: pd.DataFrame) -> None:
    top = distribution.nlargest(8, "ks_statistic").copy()
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(top["feature"][::-1], top["ks_statistic"][::-1], color="#4C72B0")
    ax.set_title("Top feature KS shifts: GSE94533 vs GSE106817")
    ax.set_xlabel("KS statistic")
    ax.set_ylabel("miRNA")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(OUTDIR / "distribution_shift.png", dpi=200)
    plt.close(fig)


def save_threshold_plot(summary: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    for _, row in summary.iterrows():
        ax.scatter(
            row["threshold"],
            row["sensitivity"],
            s=50,
            marker="o",
            label=f"{row['normalization']} sens",
            alpha=0.75,
        )
        ax.scatter(
            row["threshold"],
            row["specificity"],
            s=50,
            marker="x",
            label=f"{row['normalization']} spec",
            alpha=0.75,
        )
    ax.axvline(0.5, color="gray", linestyle="--", linewidth=1, alpha=0.7)
    ax.set_title("Frozen source threshold and operating characteristics on GSE106817")
    ax.set_xlabel("Threshold")
    ax.set_ylabel("Sensitivity / specificity")
    ax.grid(alpha=0.2)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTDIR / "threshold_analysis.png", dpi=200)
    plt.close(fig)


def make_report(summary: pd.DataFrame, distribution: pd.DataFrame, manifest: dict) -> str:
    best_auc = summary.sort_values("roc_auc", ascending=False).iloc[0]
    threshold_failure = summary[summary["sensitivity"] < 0.05].sort_values("sensitivity").iloc[0]
    top_shift = distribution.nlargest(5, "ks_statistic")
    text = f"""# Transportability analysis summary

This research-only analysis uses only the saved OVAI artifacts and does not refit the canonical model or adjust any threshold on external labels. The canonical source artifact remains the frozen production pipeline recorded in [production_artifacts/manifest.json](../../production_artifacts/manifest.json).

## Summary

- Source model: GSE94533 Stage I/II cancer vs benign/control logistic model, encoded in [production_artifacts/manifest.json](../../production_artifacts/manifest.json).
- External stress-test cohort: GSE106817, evaluated using the saved phase-5 results in [reports/gse106817_phase5_results.json](../../reports/gse106817_phase5_results.json).
- Key conclusion: the external failure is not primarily a lack of ranking signal; it is a combination of distribution shift and a probability/threshold mismatch. The ranking signal remains measurable, but the frozen source threshold collapses under cohort shift.

## Evidence from the saved outputs

The best retained GSE94533 -> GSE106817 ranking signal in the saved analysis is:

- normalization: {best_auc['normalization']}
- ROC-AUC: {best_auc['roc_auc']:.4f}
- balanced accuracy: {best_auc['balanced_accuracy']:.4f}
- sensitivity: {best_auc['sensitivity']:.4f}
- specificity: {best_auc['specificity']:.4f}
- threshold: {best_auc['threshold']:.2f}

The strongest evidence of threshold failure is:

- normalization: {threshold_failure['normalization']}
- threshold: {threshold_failure['threshold']:.2f}
- sensitivity: {threshold_failure['sensitivity']:.4f}
- specificity: {threshold_failure['specificity']:.4f}
- predicted positive rate: {threshold_failure['positive_rate']:.4f}

This is exactly the pattern expected when the score distribution shifts but the ranking relationship remains partially preserved.

## Distribution shift diagnostics

The saved shift audit in [reports/gse106817_platform_shift_audit.json](../../reports/gse106817_platform_shift_audit.json) and [reports/gse106817_phase5_distribution.csv](../../reports/gse106817_phase5_distribution.csv) shows severe feature-level shift. The top five KS shifts are:

| feature | KS statistic | median shift (log1p) |
|---|---:|---:|
"""
    for _, row in top_shift.iterrows():
        text += f"| {row['feature']} | {row['ks_statistic']:.3f} | {row['median_shift_log1p_gse106817_minus_gse94533']:.3f} |\n"
    text += """

## Interpretation

- The model still has measurable discrimination on GSE106817, especially after source-fitted scaling; the problem is not pure ranking failure.
- The canonical source threshold cannot be transported unchanged. It causes near-all-negative predictions in the external cohort.
- This is consistent with a platform/assay shift: the miRNA distributions differ substantially, as shown by KS statistics near 1.0 for several features in the audit.
- Therefore, the current computational evidence indicates that no normalization alone reliably fixes transportability without a standardized assay or a truly comparable validation cohort.

## Decision

Did we actually fix transportability, or did we only characterize it?

We only characterized it. The evidence supports a domain-shift bottleneck, not a solved transportability problem. The existing external cohorts are stress tests, not evidence that the canonical model is transportable across assay/platform boundaries.

## Reproducibility metadata

The manifest recorded by this analysis is in [research/transportability_analysis/transportability_manifest.json](transportability_manifest.json).
"""
    return text


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    summary = summarize_experiments()
    distribution = pd.read_csv(REPORTS / "gse106817_phase5_distribution.csv")
    manifest = {
        "analysis_name": "transportability_characterization",
        "root": str(ROOT),
        "source_artifacts": [
            "production_artifacts/manifest.json",
            "reports/gse106817_phase5_results.json",
            "reports/gse106817_phase5_distribution.csv",
            "reports/gse106817_platform_shift_audit.json",
        ],
        "pipeline_policy": "No external labels were used for model selection, feature tuning, or threshold selection.",
        "transformations_evaluated": [
            "raw_values",
            "log1p_values",
            "within_dataset_standard",
            "training_derived_standard",
            "within_dataset_quantile",
            "training_derived_quantile",
        ],
        "seed": 20260919,
        "notes": "This study uses only saved OVAI outputs. It does not refit canonical models or claim improved transportability.",
    }

    summary.to_csv(OUTDIR / "transportability_summary.csv", index=False)
    (OUTDIR / "transportability_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    save_distribution_shift_plot(distribution)
    save_threshold_plot(summary)

    report_text = make_report(summary, distribution, manifest)
    (OUTDIR / "transportability_report.md").write_text(report_text + "\n")

    print(f"Saved summary to {OUTDIR / 'transportability_summary.csv'}")
    print(f"Saved report to {OUTDIR / 'transportability_report.md'}")
    print(f"Saved distribution shift plot to {OUTDIR / 'distribution_shift.png'}")
    print(f"Saved threshold plot to {OUTDIR / 'threshold_analysis.png'}")


if __name__ == "__main__":
    main()
