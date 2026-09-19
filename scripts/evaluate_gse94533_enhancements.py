#!/usr/bin/env python3

"""Test leakage-safe algorithmic enhancements on the GSE94533 fusion cohort."""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier, StackingClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import KNNImputer, SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    roc_auc_score,
)
from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC, LinearSVC
from xgboost import XGBClassifier


PROJECT = Path(__file__).resolve().parents[1]
DATA = PROJECT / "data" / "processed" / "gse94533_14mirna_clinical.csv.gz"
OUT_RESULTS = PROJECT / "reports" / "tables" / "gse94533_enhancement_results.json"
OUT_PRED = PROJECT / "reports" / "tables" / "gse94533_enhancement_oof_predictions.csv"
MIRNAS = [
    "hsa-miR-23b-3p", "hsa-miR-29a-3p", "hsa-miR-32-5p", "hsa-miR-92a-3p",
    "hsa-miR-150-5p", "hsa-miR-200a-3p", "hsa-miR-200c-3p", "hsa-miR-203a",
    "hsa-miR-320c", "hsa-miR-320d", "hsa-miR-335-5p", "hsa-miR-450b-5p",
    "hsa-miR-1246", "hsa-miR-1307-5p",
]


def metrics(y, probability, threshold=0.5):
    prediction = (probability >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y, prediction).ravel()
    return {
        "accuracy": float(accuracy_score(y, prediction)),
        "balanced_accuracy": float(balanced_accuracy_score(y, prediction)),
        "roc_auc": float(roc_auc_score(y, probability)),
        "brier_score": float(brier_score_loss(y, probability)),
        "threshold": float(threshold),
        "sensitivity": float(tp / (tp + fn)),
        "specificity": float(tn / (tn + fp)),
        "confusion_matrix": {"tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp)},
    }


def clinical_features(frame):
    output = frame[MIRNAS + ["age", "ca125"]].copy()
    output[MIRNAS] = np.log1p(output[MIRNAS].clip(lower=0))
    output["ratio_200c_29a"] = output["hsa-miR-200c-3p"] - output["hsa-miR-29a-3p"]
    output["ratio_320d_29a"] = output["hsa-miR-320d"] - output["hsa-miR-29a-3p"]
    output["ratio_1246_29a"] = output["hsa-miR-1246"] - output["hsa-miR-29a-3p"]
    output["ratio_200c_320d"] = output["hsa-miR-200c-3p"] - output["hsa-miR-320d"]
    output["ca125_log"] = np.log1p(output["ca125"].clip(lower=0))
    output["age_squared"] = output["age"] ** 2
    output["age_ca125"] = output["age"] * output["ca125_log"]
    return output


def logistic_pipeline(imputer):
    return Pipeline([
        ("impute", imputer),
        ("scale", StandardScaler()),
        ("select", SelectKBest(f_classif)),
        ("model", LogisticRegression(max_iter=5000, solver="liblinear")),
    ])


def fit_search(estimator, grid, X_train, y_train, inner):
    search = GridSearchCV(
        estimator, grid, scoring="balanced_accuracy", cv=inner,
        n_jobs=-1, refit=True, error_score="raise",
    )
    search.fit(X_train, y_train)
    return search


def main():
    warnings.filterwarnings("ignore", category=UserWarning)
    frame = pd.read_csv(DATA)
    X = clinical_features(frame)
    y = frame["binary_label"].astype(int).to_numpy()
    feature_names = list(X.columns)
    ratio_names = [name for name in feature_names if name.startswith("ratio_")]
    X_base = X.drop(columns=ratio_names)
    outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    inner = StratifiedKFold(n_splits=3, shuffle=True, random_state=43)
    model_probabilities = {}
    fold_metrics = {}
    lasso_selection = []
    threshold_predictions = {}
    selected_thresholds = {}
    specificity_predictions = {}
    specificity_thresholds = {}

    for fold, (train_idx, test_idx) in enumerate(outer.split(X, y), start=1):
        X_train, X_test = X_base.iloc[train_idx], X_base.iloc[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        k_values = [k for k in [5, 10, 14, 16] if k <= X_base.shape[1]] + ["all"]

        candidates = {
            "median_logistic": (
                logistic_pipeline(SimpleImputer(strategy="median")),
                {"select__k": k_values, "model__C": [0.03, 0.1, 0.3, 1.0], "model__class_weight": [None, "balanced"]},
            ),
            "knn_logistic": (
                logistic_pipeline(KNNImputer(n_neighbors=5, weights="distance")),
                {"select__k": k_values, "model__C": [0.03, 0.1, 0.3, 1.0], "model__class_weight": [None, "balanced"]},
            ),
            "calibrated_svm": (
                Pipeline([
                    ("impute", SimpleImputer(strategy="median")),
                    ("scale", StandardScaler()),
                    ("model", CalibratedClassifierCV(
                        estimator=LinearSVC(class_weight="balanced", dual="auto", max_iter=10000),
                        method="sigmoid", cv=3,
                    )),
                ]),
                {"model__estimator__C": [0.03, 0.1, 0.3, 1.0]},
            ),
            "random_forest": (
                Pipeline([
                    ("impute", SimpleImputer(strategy="median")),
                    ("model", RandomForestClassifier(n_estimators=300, random_state=42, n_jobs=-1)),
                ]),
                {"model__max_depth": [2, 3, 5], "model__min_samples_leaf": [2, 5, 10], "model__class_weight": [None, "balanced"]},
            ),
            "xgboost": (
                Pipeline([
                    ("impute", SimpleImputer(strategy="median")),
                    ("model", XGBClassifier(
                        n_estimators=200, learning_rate=0.05, subsample=0.8,
                        colsample_bytree=0.8, objective="binary:logistic",
                        eval_metric="logloss", random_state=42, n_jobs=2,
                    )),
                ]),
                {"model__max_depth": [1, 2, 3], "model__min_child_weight": [3, 8]},
            ),
        }
        ratio_estimator = logistic_pipeline(SimpleImputer(strategy="median"))
        ratio_grid = {
            "select__k": [k for k in [5, 10, 14, 16, 20] if k <= len(feature_names)] + ["all"],
            "model__C": [0.03, 0.1, 0.3, 1.0],
            "model__class_weight": [None, "balanced"],
        }
        candidates["ratio_logistic"] = (ratio_estimator, ratio_grid)

        base_estimators = [
            ("xgb", candidates["xgboost"][0]),
            ("rf", candidates["random_forest"][0]),
            ("svm", candidates["calibrated_svm"][0]),
        ]
        candidates["stacked_ensemble"] = (
            StackingClassifier(
                estimators=base_estimators,
                final_estimator=LogisticRegression(C=0.3, class_weight="balanced", max_iter=5000),
                stack_method="predict_proba", cv=3, n_jobs=-1,
            ),
            {},
        )

        for name, (estimator, grid) in candidates.items():
            candidate_train = X.iloc[train_idx] if name == "ratio_logistic" else X_train
            candidate_test = X.iloc[test_idx] if name == "ratio_logistic" else X_test
            search = fit_search(estimator, grid, candidate_train, y_train, inner) if grid else estimator.fit(candidate_train, y_train)
            probability = search.predict_proba(candidate_test)[:, 1]
            model_probabilities.setdefault(name, np.full(len(X), np.nan))[test_idx] = probability
            fold_metrics.setdefault(name, []).append({
                "fold": fold,
                **metrics(y_test, probability),
            })
            if name == "calibrated_svm":
                train_probability = cross_val_predict(
                    search.best_estimator_, candidate_train, y_train,
                    cv=inner, method="predict_proba", n_jobs=-1,
                )[:, 1]
                thresholds = np.linspace(0.10, 0.90, 81)
                threshold = max(
                    thresholds,
                    key=lambda candidate: balanced_accuracy_score(
                        y_train, (train_probability >= candidate).astype(int)
                    ),
                )
                selected_thresholds.setdefault(name, []).append(float(threshold))
                threshold_predictions.setdefault(name, np.full(len(X), np.nan))[test_idx] = (
                    probability >= threshold
                ).astype(int)
                specificity_candidates = []
                for candidate in thresholds:
                    candidate_prediction = (train_probability >= candidate).astype(int)
                    candidate_metrics = metrics(y_train, train_probability, threshold=candidate)
                    specificity_candidates.append((candidate, candidate_metrics["specificity"], candidate_metrics["sensitivity"]))
                eligible = [row for row in specificity_candidates if row[1] >= 0.95]
                selected_specificity = max(
                    eligible or specificity_candidates,
                    key=lambda row: (row[2], row[1]),
                )[0]
                specificity_thresholds.setdefault(name, []).append(float(selected_specificity))
                specificity_predictions.setdefault(name, np.full(len(X), np.nan))[test_idx] = (
                    probability >= selected_specificity
                ).astype(int)

        lasso = Pipeline([
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            ("model", LogisticRegression(penalty="l1", solver="liblinear", C=0.1, class_weight="balanced", max_iter=5000)),
        ])
        lasso.fit(X_train, y_train)
        selected = [name for name, weight in zip(feature_names, lasso.named_steps["model"].coef_[0]) if abs(weight) > 1e-10]
        lasso_selection.append({"fold": fold, "selected_features": selected})

    results = {}
    for name, probability in model_probabilities.items():
        results[name] = {
            **metrics(y, probability),
            "folds": fold_metrics[name],
        }
        if name in threshold_predictions:
            results[name]["inner_selected_thresholds"] = selected_thresholds[name]
            report_threshold = float(np.median(selected_thresholds[name]))
            results[name]["threshold_tuned_metrics"] = metrics(
                y, probability, threshold=report_threshold
            )
            tuned_prediction = threshold_predictions[name].astype(int)
            tn, fp, fn, tp = confusion_matrix(y, tuned_prediction).ravel()
            results[name]["threshold_tuned_metrics"].update({
                "accuracy": float(accuracy_score(y, tuned_prediction)),
                "balanced_accuracy": float(balanced_accuracy_score(y, tuned_prediction)),
                "sensitivity": float(tp / (tp + fn)),
                "specificity": float(tn / (tn + fp)),
            })
            specificity_prediction = specificity_predictions[name].astype(int)
            tn, fp, fn, tp = confusion_matrix(y, specificity_prediction).ravel()
            results[name]["specificity_target_metrics"] = {
                "inner_selected_thresholds": specificity_thresholds[name],
                "accuracy": float(accuracy_score(y, specificity_prediction)),
                "balanced_accuracy": float(balanced_accuracy_score(y, specificity_prediction)),
                "sensitivity": float(tp / (tp + fn)),
                "specificity": float(tn / (tn + fp)),
                "confusion_matrix": {"tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp)},
            }
        print(name, results[name]["accuracy"], results[name]["balanced_accuracy"], results[name]["roc_auc"])

    stability = {}
    for feature in feature_names:
        selected_folds = [row["fold"] for row in lasso_selection if feature in row["selected_features"]]
        stability[feature] = {"folds": selected_folds, "frequency": len(selected_folds) / len(lasso_selection)}

    output = frame[["sample_id_normalized", "binary_label"]].copy()
    for name, probability in model_probabilities.items():
        output[name + "_probability"] = probability
        output[name + "_prediction"] = (probability >= 0.5).astype(int)
    OUT_PRED.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(OUT_PRED, index=False)
    OUT_RESULTS.write_text(json.dumps({
        "endpoint": "GSE94533 invasive ovarian cancer vs non-cancer",
        "patients": len(frame),
        "evaluation": "nested 5-fold out-of-fold predictions; all tuning inside training folds",
        "models": results,
        "lasso_stability": stability,
        "lasso_by_fold": lasso_selection,
        "external_validation": "phs000897 unavailable: dbGaP controlled access",
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()