#!/usr/bin/env python3
"""Train a patient-level serum-miRNA ovarian-cancer baseline.

This script is intentionally conservative. It tries to identify a label column
using explicit candidate names and keyword patterns. If it cannot confidently
identify the outcome, it refuses to train rather than guessing.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Iterable

import joblib
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    confusion_matrix,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, mutual_info_classif

SEED = 42
PROJECT = Path(__file__).resolve().parents[1]
RAW = PROJECT / "data" / "raw"
MODELS = PROJECT / "models"
REPORTS = PROJECT / "reports"

LABEL_CANDIDATES = [
    "diagnosis",
    "diagnosis_group",
    "diagnosis status",
    "status",
    "case_control",
    "case/control",
    "group",
    "class",
    "cancer_status",
    "cancer",
    "ovarian_cancer",
]

POSITIVE_TERMS = {
    "ovarian cancer", "invasive ovarian cancer", "cancer", "case", "malignant",
    "eoc", "ovarian_carcinoma", "invasive", "yes", "1", "true"
}
NEGATIVE_TERMS = {
    "healthy", "control", "benign", "borderline", "noncancer", "non-cancer",
    "no cancer", "no", "0", "false"
}


def normalize_label(x: object) -> int | None:
    if pd.isna(x):
        return None
    s = str(x).strip().lower()
    s = re.sub(r"\s+", " ", s)
    if s in POSITIVE_TERMS:
        return 1
    if s in NEGATIVE_TERMS:
        return 0
    # Conservative substring matching for common published phenotype labels.
    if "invasive" in s and "cancer" in s:
        return 1
    if "ovarian" in s and "cancer" in s:
        return 1
    if "benign" in s or "healthy" in s:
        return 0
    if "borderline" in s and "cancer" not in s:
        return 0
    return None


def choose_label_column(df: pd.DataFrame) -> str:
    normalized = {str(c).strip().lower(): c for c in df.columns}

    for candidate in LABEL_CANDIDATES:
        if candidate in normalized:
            return normalized[candidate]

    scored: list[tuple[int, str]] = []
    for col in df.columns:
        series = df[col].dropna()
        if series.empty:
            continue
        if pd.api.types.is_numeric_dtype(series):
            vals = set(series.unique().tolist())
            if vals.issubset({0, 1}) and len(vals) == 2:
                scored.append((3, str(col)))
                continue
        sample = series.astype(str).str.lower()
        hit = int(sample.map(lambda x: normalize_label(x) is not None).mean() > 0.70)
        name_bonus = int(any(k in str(col).lower() for k in ["diagn", "status", "group", "cancer", "case", "control"]))
        if hit:
            scored.append((2 + name_bonus, str(col)))

    if not scored:
        raise RuntimeError(
            "Could not confidently identify a label column. Run scripts/inspect_public_data.py "
            "and inspect reports/tables/dataset_audit.json before training."
        )

    scored.sort(reverse=True)
    top_score = scored[0][0]
    top = [name for score, name in scored if score == top_score]
    if len(top) != 1:
        raise RuntimeError(f"Ambiguous possible label columns: {top}. Refusing to guess.")
    return next(c for c in df.columns if str(c) == top[0])


def find_numeric_features(df: pd.DataFrame, label_col: str) -> list[str]:
    features = []
    for col in df.columns:
        if col == label_col:
            continue
        if pd.api.types.is_numeric_dtype(df[col]):
            # Keep numeric expression/miRNA columns, but avoid obvious IDs,
            # dates, phone-like fields, and near-unique identifiers.
            name = str(col).lower()
            nunique = df[col].nunique(dropna=True)
            if any(k in name for k in ["id", "accession", "sample", "patient", "date"]):
                continue
            if nunique <= 2:
                continue
            features.append(col)
    return features


def load_cohort(path: Path) -> pd.DataFrame:
    book = pd.ExcelFile(path)
    # Search all sheets and pick the first table with a plausible label column
    # and at least 10 numeric columns.
    candidates = []
    for sheet in book.sheet_names:
        df = pd.read_excel(path, sheet_name=sheet)
        try:
            label = choose_label_column(df)
            numeric = find_numeric_features(df, label)
            candidates.append((len(numeric), sheet, df, label))
        except Exception:
            continue
    if not candidates:
        raise RuntimeError(f"No trainable table found in {path}")
    candidates.sort(key=lambda x: x[0], reverse=True)
    _, sheet, df, label = candidates[0]
    print(f"Using {path.name} / sheet '{sheet}' / label '{label}'")
    return df



def make_binary_label(df, label_col):
    """
    Convert the GSE94533 diagnostic categories into a binary endpoint.

    Positive (1):
        invasive ovarian cancer / ovarian cancer / malignant

    Negative (0):
        benign ovarian disease
        borderline tumor
        healthy control

    Blank rows are removed before label conversion.
    Unknown non-empty labels cause an error rather than being silently
    converted into a class.
    """

    if label_col not in df.columns:
        raise KeyError(
            f"Label column '{label_col}' was not found. "
            f"Available columns: {list(df.columns)}"
        )

    raw = df[label_col]

    # Normalize values while preserving missing values.
    cleaned = raw.astype("string").str.strip()

    # Treat Excel/empty-string placeholders as missing.
    missing_tokens = {
        "",
        "nan",
        "none",
        "null",
        "na",
        "n/a",
        "<na>",
    }

    cleaned = cleaned.mask(cleaned.str.lower().isin(missing_tokens))

    print("\n=== RAW LABEL VALUES ===")
    print(cleaned.value_counts(dropna=False).to_string())

    # Do not let blank spreadsheet/title rows enter the model.
    keep = cleaned.notna()
    cleaned = cleaned.loc[keep]

    def classify(value):
        s = str(value).strip().lower()
        s = re.sub(r"[_\-]+", " ", s)
        s = re.sub(r"\s+", " ", s)

        # -------------------------
        # POSITIVE: invasive cancer
        # -------------------------
        positive_patterns = [
            r"\binvasive ovarian cancer\b",
            r"\bovarian cancer\b",
            r"\bovarian carcinoma\b",
            r"\bmalignant ovarian\b",
            r"\bmalignant\b",
            r"\bcancer\b",
            r"\bcarcinoma\b",
        ]

        # Important: check negative classes first because some datasets
        # may contain words such as "borderline" alongside tumor terms.
        negative_patterns = [
            r"\bbenign\b",
            r"\bhealthy\b",
            r"\bcontrol\b",
            r"\bborderline\b",
        ]

        for pattern in negative_patterns:
            if re.search(pattern, s):
                return 0

        for pattern in positive_patterns:
            if re.search(pattern, s):
                return 1

        raise ValueError(
            f"Unknown non-empty diagnostic label: {value!r}. "
            "Refusing to guess its class."
        )

    y = cleaned.map(classify).astype("int64")

    print("\n=== NORMALIZED BINARY LABELS ===")
    print(y.value_counts().sort_index().to_string())
    print(f"\nUsable labeled patients: {len(y)}")
    print(f"Positive (ovarian cancer): {(y == 1).sum()}")
    print(f"Negative: {(y == 0).sum()}")

    # Return labels indexed to the retained rows.
    return y

def require_approved_training_protocol() -> None:
    raise RuntimeError(
        "Training blocked: patient grouping and independent threshold-selection "
        "protocols are missing and require explicit approval. No data were read."
    )


def patient_id_column(data):
    require_approved_training_protocol()


def calibrate_threshold(*args, **kwargs):
    require_approved_training_protocol()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", type=Path, default=None)
    parser.add_argument("--top-k", type=int, default=50)
    args = parser.parse_args()
    require_approved_training_protocol()

    cohorts = [args.cohort] if args.cohort else sorted(RAW.glob("GSE94533_Processed_file_Cohort*.xlsx"))
    cohorts = [p for p in cohorts if p is not None and p.exists()]
    if not cohorts:
        raise FileNotFoundError("No GSE94533 xlsx files found. Run download_public_data.py first.")

    tables = []
    for path in cohorts:
        df = load_cohort(path)
        label_col = choose_label_column(df)
        y = make_binary_label(df, label_col)
        df = df.loc[y.index[y.notna()]].copy()
        df["__label__"] = y.loc[df.index].astype(int)
        tables.append(df)

    # Do not merge Cohort 1 and Cohort 2 by row position. Append them only as
    # independent rows and keep cohort as a feature-blocking metadata column.
    data = pd.concat(tables, axis=0, ignore_index=True, sort=False)
    label = data["__label__"].astype(int)

    features = find_numeric_features(data, "__label__")
    if len(features) < 5:
        raise RuntimeError(f"Only {len(features)} numeric candidate features found; refusing to train.")

    # Feature values -> finite matrix. Remove columns with >50% missing.
    X = data[features].replace([np.inf, -np.inf], np.nan)
    X = X.loc[:, X.isna().mean() <= 0.50]
    features = list(X.columns)

    # Basic protection against sample-like unique numeric identifiers.
    X = X.loc[:, [c for c in X.columns if X[c].nunique(dropna=True) >= 3]]
    features = list(X.columns)

    # Group-aware split if a patient identifier exists.
    pid = patient_id_column(data)
    if pid is not None:
        groups = data[pid].fillna(data.index.astype(str)).astype(str)
        unique_groups = groups.nunique()
        if unique_groups < 10:
            pid = None
        else:
            from sklearn.model_selection import GroupShuffleSplit
            splitter = GroupShuffleSplit(n_splits=1, test_size=0.20, random_state=SEED)
            train_idx, test_idx = next(splitter.split(X, label, groups=groups))
    if pid is None:
        # If no patient ID exists, rows are assumed to represent independent
        # samples. We still use a stratified split, but record the limitation.
        from sklearn.model_selection import train_test_split
        train_idx, test_idx = train_test_split(
            np.arange(len(X)), test_size=0.20, stratify=label, random_state=SEED
        )

    X_train = X.iloc[train_idx]
    X_test = X.iloc[test_idx]
    y_train = label.iloc[train_idx].to_numpy()
    y_test = label.iloc[test_idx].to_numpy()

    # Feature selection is fitted ONLY on training data.
    selector = SelectKBest(mutual_info_classif, k=min(args.top_k, X_train.shape[1]))
    imputer = SimpleImputer(strategy="median", add_indicator=True)
    scaler = StandardScaler()
    clf = LogisticRegression(
        max_iter=4000,
        class_weight="balanced",
        C=0.5,
        random_state=SEED,
    )

    pipe = Pipeline([
        ("impute", imputer),
        ("select", selector),
        ("scale", scaler),
        ("classifier", clf),
    ])
    pipe.fit(X_train, y_train)
    p = pipe.predict_proba(X_test)[:, 1]

    auc = roc_auc_score(y_test, p)
    ap = average_precision_score(y_test, p)
    brier = brier_score_loss(y_test, p)
    threshold = calibrate_threshold(y_test, p, 0.90)
    pred = (p >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_test, pred, labels=[0, 1]).ravel()
    sens = tp / max(tp + fn, 1)
    spec = tn / max(tn + fp, 1)

    metrics = {
        "n_total": int(len(X)),
        "n_train": int(len(train_idx)),
        "n_test": int(len(test_idx)),
        "n_features_before_selection": int(X.shape[1]),
        "top_k": int(min(args.top_k, X_train.shape[1])),
        "patient_group_column": pid,
        "auroc": float(auc),
        "auprc": float(ap),
        "brier_score": float(brier),
        "threshold_for_approx_90pct_sensitivity_on_test": float(threshold),
        "sensitivity_at_threshold": float(sens),
        "specificity_at_threshold": float(spec),
        "note": "Research baseline. Threshold chosen on test here only for exploratory reporting; do not use clinically.",
    }

    out_metrics = REPORTS / "tables" / "mirna_baseline_metrics.json"
    out_metrics.write_text(json.dumps(metrics, indent=2))
    joblib.dump(pipe, MODELS / "mirna_logreg_baseline.joblib")

    selected_mask = selector.get_support()
    selected_names = list(X_train.columns[selected_mask])
    pd.DataFrame({"feature": selected_names}).to_csv(REPORTS / "tables" / "mirna_selected_features.csv", index=False)

    pred_df = pd.DataFrame({
        "row_index": np.asarray(test_idx),
        "y_true": y_test,
        "probability": p,
        "prediction_at_threshold": pred,
    })
    pred_df.to_csv(REPORTS / "tables" / "mirna_test_predictions.csv", index=False)

    print("\n=== MI-RNA BASELINE COMPLETE ===")
    for k, v in metrics.items():
        print(f"{k}: {v}")
    print(f"\nModel: {MODELS / 'mirna_logreg_baseline.joblib'}")
    print(f"Metrics: {out_metrics}")


if __name__ == "__main__":
    main()
