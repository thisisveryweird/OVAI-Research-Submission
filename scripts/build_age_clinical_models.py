#!/usr/bin/env python3

from pathlib import Path
import json
import warnings

import numpy as np
import pandas as pd

from sklearn.compose import ColumnTransformer
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

ROOT = Path(__file__).resolve().parents[1]

DATA = (
    ROOT
    / "data"
    / "processed"
    / "gse94533_14mirna_clinical_benchmark.csv.gz"
)

OUTDIR = ROOT / "reports" / "tables"

MIRNAS = [
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


def run_model(df, features, name):

    X = df[features]
    y = df["binary_label"].astype(int)

    numeric_pipe = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        (
            "model",
            LogisticRegression(
                C=1.0,
                max_iter=5000,
                solver="liblinear",
            ),
        ),
    ])

    cv = StratifiedKFold(
        n_splits=5,
        shuffle=True,
        random_state=42,
    )

    oof_prob = np.full(len(df), np.nan)
    oof_pred = np.full(len(df), -1, dtype=int)

    for train_idx, test_idx in cv.split(X, y):

        model = numeric_pipe

        model.fit(
            X.iloc[train_idx],
            y.iloc[train_idx],
        )

        oof_prob[test_idx] = model.predict_proba(
            X.iloc[test_idx]
        )[:, 1]

        oof_pred[test_idx] = (
            oof_prob[test_idx] >= 0.5
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

    pred = pd.DataFrame({
        "sample_id": df["sample_id"],
        "true_label": y,
        "oof_probability": oof_prob,
        "oof_prediction": oof_pred,
    })

    pred_path = (
        OUTDIR
        / f"{name}_oof_predictions.csv"
    )

    result_path = (
        OUTDIR
        / f"{name}_results.json"
    )

    pred.to_csv(
        pred_path,
        index=False,
    )

    results = {
        "model": name,
        "patients": int(len(df)),
        "features": features,
        "roc_auc": float(auc),
        "pr_auc": float(pr_auc),
        "accuracy": float(accuracy),
        "sensitivity": float(sensitivity),
        "specificity": float(specificity),
        "tp": int(tp),
        "fn": int(fn),
        "fp": int(fp),
        "tn": int(tn),
    }

    result_path.write_text(
        json.dumps(results, indent=2)
    )

    return results


def main():

    print("=" * 80)
    print("AGE / CLINICAL OVARIAN CANCER MODELS")
    print("=" * 80)

    df = pd.read_csv(DATA)

    if len(df) != 179:
        raise RuntimeError(
            f"Expected 179 patients, found {len(df)}"
        )

    required = {
        "sample_id",
        "binary_label",
        "age",
        "ca125",
        *MIRNAS,
    }

    missing = required - set(df.columns)

    if missing:
        raise RuntimeError(
            "Missing columns: "
            + ", ".join(sorted(missing))
        )

    print()
    print("Patients:", len(df))
    print("Cancer:", int(df["binary_label"].sum()))
    print(
        "Non-cancer:",
        int((df["binary_label"] == 0).sum())
    )
    print("Age missing:", int(df["age"].isna().sum()))
    print("CA-125 missing:", int(df["ca125"].isna().sum()))

    # ---------------------------------------------------------------
    # MODEL 1
    # ---------------------------------------------------------------

    r1 = run_model(
        df,
        ["age"] + MIRNAS,
        "age_plus_14mirna",
    )

    # ---------------------------------------------------------------
    # MODEL 2
    # ---------------------------------------------------------------

    r2 = run_model(
        df,
        ["age", "ca125"] + MIRNAS,
        "age_ca125_plus_14mirna",
    )

    print()
    print("=" * 80)
    print("FINAL RESULTS")
    print("=" * 80)

    for r in [r1, r2]:

        print()
        print(r["model"])
        print(
            f"  ROC AUC:    {r['roc_auc']:.6f}"
        )
        print(
            f"  PR AUC:     {r['pr_auc']:.6f}"
        )
        print(
            f"  Accuracy:   {r['accuracy']:.6f}"
        )
        print(
            f"  Sensitivity:{r['sensitivity']:.6f}"
        )
        print(
            f"  Specificity:{r['specificity']:.6f}"
        )

    print()
    print("=" * 80)
    print("FILES SAVED")
    print("=" * 80)

    print(
        OUTDIR
        / "age_plus_14mirna_oof_predictions.csv"
    )

    print(
        OUTDIR
        / "age_ca125_plus_14mirna_oof_predictions.csv"
    )


if __name__ == "__main__":
    main()
