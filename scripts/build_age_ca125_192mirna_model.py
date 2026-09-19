#!/usr/bin/env python3

from pathlib import Path
import json
import warnings

import numpy as np
import pandas as pd

from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

ROOT = Path.home() / "ovarian-ai"

DATA = (
    ROOT
    / "data"
    / "processed"
    / "gse94533_14mirna_clinical_benchmark.csv.gz"
)

OUTDIR = ROOT / "reports" / "tables"


def main():

    print("=" * 80)
    print("AGE + CA-125 + 192-miRNA MODEL")
    print("=" * 80)

    df = pd.read_csv(DATA)

    mirnas = [
        c for c in df.columns
        if c.startswith("hsa-")
    ]

    if len(mirnas) != 14:
        raise RuntimeError(
            f"Expected 14 miRNA columns in this benchmark, found {len(mirnas)}."
        )

    # The current clinical benchmark contains only the exact 14 PMML
    # miRNAs. Build the full 192-miRNA matrix from Dataset 2.
    xlsx = (
        ROOT
        / "data"
        / "external"
        / "elife_28932"
        / "elife_poa_e28932_Supplementary_File_6.xlsx"
    )

    raw = pd.read_excel(
        xlsx,
        sheet_name="Supplemental Dataset 2",
        header=None,
        engine="openpyxl",
    )

    feature_names = (
        raw.iloc[:, 0]
        .astype(str)
        .str.strip()
    )

    all_features = [
        x for x in feature_names
        if x.startswith("hsa-")
    ]

    # Remove duplicates while preserving order.
    all_features = list(dict.fromkeys(all_features))

    if len(all_features) != 192:
        raise RuntimeError(
            f"Expected 192 miRNA features, found {len(all_features)}."
        )

    # The verified PMML cohort establishes the exact 179 column order.
    pmml = pd.read_csv(
        ROOT
        / "reports"
        / "tables"
        / "original_historical_clean_predictions.csv"
    )

    verified_ids = pmml["sample_id"].astype(str).tolist()

    if len(verified_ids) != 179:
        raise RuntimeError(
            "Verified PMML cohort is not 179 patients."
        )

    matrix = {}

    for feature in all_features:

        idx = np.where(
            feature_names.values == feature
        )[0]

        if len(idx) != 1:
            raise RuntimeError(
                f"Could not uniquely locate {feature}."
            )

        values = pd.to_numeric(
            raw.iloc[idx[0], 1:180],
            errors="coerce",
        ).to_numpy(float)

        if len(values) != 179:
            raise RuntimeError(
                f"{feature}: expected 179 values."
            )

        matrix[feature] = values

    Xmirna = pd.DataFrame(matrix)
    Xmirna["sample_id"] = verified_ids

    # Join age + verified label.
    clinical = df[
        ["sample_id", "age", "ca125", "binary_label"]
    ].copy()

    merged = Xmirna.merge(
        clinical,
        on="sample_id",
        how="inner",
        validate="one_to_one",
    )

    if len(merged) != 179:
        raise RuntimeError(
            f"Expected 179 merged patients, got {len(merged)}."
        )

    features = ["age", "ca125"] + all_features

    X = merged[features]
    y = merged["binary_label"].astype(int)

    print()
    print("Patients:", len(merged))
    print("Features:", len(features))
    print("miRNAs:", len(all_features))
    print("Cancer:", int(y.sum()))
    print("Non-cancer:", int((y == 0).sum()))

    cv = StratifiedKFold(
        n_splits=5,
        shuffle=True,
        random_state=42,
    )

    pipe = Pipeline([
        (
            "imputer",
            SimpleImputer(strategy="median"),
        ),
        (
            "scaler",
            StandardScaler(),
        ),
        (
            "model",
            LogisticRegression(
                C=0.3,
                penalty="l2",
                solver="liblinear",
                max_iter=5000,
            ),
        ),
    ])

    oof_prob = np.full(len(y), np.nan)

    for fold, (train_idx, test_idx) in enumerate(
        cv.split(X, y),
        start=1,
    ):

        pipe.fit(
            X.iloc[train_idx],
            y.iloc[train_idx],
        )

        oof_prob[test_idx] = pipe.predict_proba(
            X.iloc[test_idx]
        )[:, 1]

    oof_pred = (
        oof_prob >= 0.5
    ).astype(int)

    auc = roc_auc_score(y, oof_prob)
    pr_auc = average_precision_score(y, oof_prob)
    accuracy = accuracy_score(y, oof_pred)

    tn, fp, fn, tp = confusion_matrix(
        y,
        oof_pred,
    ).ravel()

    sensitivity = tp / (tp + fn)
    specificity = tn / (tn + fp)

    print()
    print("=" * 80)
    print("RESULTS")
    print("=" * 80)

    print(f"OOF ROC AUC:    {auc:.6f}")
    print(f"OOF PR AUC:     {pr_auc:.6f}")
    print(f"OOF Accuracy:   {accuracy:.6f}")
    print(f"Sensitivity:    {sensitivity:.6f}")
    print(f"Specificity:    {specificity:.6f}")

    print()
    print("Confusion matrix:")
    print(f"  TP={tp} FN={fn}")
    print(f"  FP={fp} TN={tn}")

    predictions = pd.DataFrame({
        "sample_id": merged["sample_id"],
        "true_label": y,
        "oof_probability": oof_prob,
        "oof_prediction": oof_pred,
    })

    pred_path = (
        OUTDIR
        / "age_ca125_plus_192mirna_oof_predictions.csv"
    )

    result_path = (
        OUTDIR
        / "age_ca125_plus_192mirna_results.json"
    )

    predictions.to_csv(
        pred_path,
        index=False,
    )

    result_path.write_text(
        json.dumps(
            {
                "model": "Age + CA-125 + 192 miRNA",
                "patients": len(merged),
                "miRNA_features": 192,
                "roc_auc": float(auc),
                "pr_auc": float(pr_auc),
                "accuracy": float(accuracy),
                "sensitivity": float(sensitivity),
                "specificity": float(specificity),
                "tp": int(tp),
                "fn": int(fn),
                "fp": int(fp),
                "tn": int(tn),
            },
            indent=2,
        )
    )

    print()
    print("Saved:")
    print(pred_path)
    print(result_path)


if __name__ == "__main__":
    main()
