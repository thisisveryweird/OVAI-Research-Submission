from pathlib import Path
import numpy as np
import pandas as pd

from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    accuracy_score,
    balanced_accuracy_score,
    recall_score,
    confusion_matrix,
)

ROOT = Path("data/external/mendeley_ovarian/raw")
REPORT = Path("reports/tables")
REPORT.mkdir(parents=True, exist_ok=True)

train_file = next(ROOT.rglob("Supplementary data 4.xlsx"))
test_file = next(ROOT.rglob("Supplementary data 5.xlsx"))

train = pd.read_excel(train_file)
test = pd.read_excel(test_file)

# TYPE:
# 0 = ovarian cancer
# 1 = benign ovarian tumor
#
# We define y=1 as ovarian cancer so that:
# sensitivity = cancer sensitivity
# specificity = BOT specificity
y_train = (train["TYPE"] == 0).astype(int)
y_test = (test["TYPE"] == 0).astype(int)

print("=" * 100)
print("MENDELEY CLINICAL BASELINES")
print("=" * 100)

print("TRAIN:", train.shape)
print("TEST :", test.shape)

print("\nTRAIN LABELS:")
print(y_train.value_counts().sort_index().to_string())

print("\nTEST LABELS:")
print(y_test.value_counts().sort_index().to_string())

# Convert all candidate biomarkers safely to numeric.
def numeric_frame(df, cols):
    out = df[cols].copy()

    for c in cols:
        out[c] = (
            out[c]
            .astype(str)
            .str.strip()
            .str.replace("\t", "", regex=False)
            .str.replace(",", "", regex=False)
            .str.replace(r"^[<>]=?", "", regex=True)
        )

        out[c] = pd.to_numeric(out[c], errors="coerce")

    return out


feature_sets = {
    "age_menopause": [
        "Age",
        "Menopause",
    ],

    "ca125_he4": [
        "CA125",
        "HE4",
    ],

    "age_menopause_ca125_he4": [
        "Age",
        "Menopause",
        "CA125",
        "HE4",
    ],

    "clinical_biochemical": [
        "AFP",
        "AG",
        "Age",
        "ALB",
        "ALP",
        "ALT",
        "AST",
        "BASO#",
        "BASO%",
        "BUN",
        "Ca",
        "CA125",
        "CA19-9",
        "CA72-4",
        "CEA",
        "CL",
        "CO2CP",
        "CREA",
        "DBIL",
        "EO#",
        "EO%",
        "GGT",
        "GLO",
        "GLU.",
        "HCT",
        "HE4",
        "HGB",
        "IBIL",
        "K",
        "LYM#",
        "LYM%",
        "MCH",
        "MCV",
        "Menopause",
        "Mg",
        "MONO#",
        "MONO%",
        "MPV",
        "Na",
        "NEU",
        "PCT",
        "PDW",
        "PHOS",
        "PLT",
        "RBC",
        "RDW",
        "TBIL",
        "TP",
        "UA",
    ],
}

results = []

for name, cols in feature_sets.items():

    print("\n" + "=" * 100)
    print("MODEL:", name)
    print("FEATURES:", ", ".join(cols))
    print("=" * 100)

    X_train = numeric_frame(train, cols)
    X_test = numeric_frame(test, cols)

    print("\nTraining missingness:")
    print(
        X_train.isna()
        .sum()
        .sort_values(ascending=False)
        .to_string()
    )

    # Median imputation and scaling are fitted ONLY on training data.
    model = Pipeline([
        (
            "imputer",
            SimpleImputer(strategy="median"),
        ),
        (
            "scaler",
            StandardScaler(),
        ),
        (
            "classifier",
            LogisticRegression(
                max_iter=5000,
                class_weight="balanced",
                solver="liblinear",
            ),
        ),
    ])

    model.fit(X_train, y_train)

    p_train = model.predict_proba(X_train)[:, 1]
    p_test = model.predict_proba(X_test)[:, 1]

    pred_test = (p_test >= 0.5).astype(int)

    tn, fp, fn, tp = confusion_matrix(
        y_test,
        pred_test,
        labels=[0, 1],
    ).ravel()

    specificity = tn / (tn + fp) if (tn + fp) else np.nan
    sensitivity = tp / (tp + fn) if (tp + fn) else np.nan

    result = {
        "model": name,
        "n_train": len(train),
        "n_test": len(test),
        "n_features": len(cols),

        "train_auc": roc_auc_score(y_train, p_train),
        "test_auc": roc_auc_score(y_test, p_test),

        "test_pr_auc": average_precision_score(y_test, p_test),

        "test_accuracy": accuracy_score(
            y_test,
            pred_test,
        ),

        "test_balanced_accuracy": balanced_accuracy_score(
            y_test,
            pred_test,
        ),

        "test_sensitivity_cancer": sensitivity,
        "test_specificity_bot": specificity,

        "tn": tn,
        "fp": fp,
        "fn": fn,
        "tp": tp,
    }

    results.append(result)

    print("\nRESULTS:")
    for k, v in result.items():
        if k != "model":
            print(f"{k}: {v}")

    print("\nConfusion matrix [[TN, FP], [FN, TP]]:")
    print(np.array([[tn, fp], [fn, tp]]))

results_df = pd.DataFrame(results)

output = REPORT / "mendeley_clinical_baselines.csv"
results_df.to_csv(output, index=False)

print("\n" + "=" * 100)
print("SUMMARY")
print("=" * 100)

print(
    results_df[
        [
            "model",
            "test_auc",
            "test_pr_auc",
            "test_accuracy",
            "test_balanced_accuracy",
            "test_sensitivity_cancer",
            "test_specificity_bot",
        ]
    ].to_string(index=False)
)

print("\nSaved:", output)
print("\nDONE")
