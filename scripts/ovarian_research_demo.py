from pathlib import Path
import pandas as pd

ROOT = Path.home() / "ovarian-ai"

PRED = (
    ROOT
    / "reports"
    / "tables"
    / "age_plus_192mirna_oof_predictions.csv"
)

df = pd.read_csv(PRED)

df["actual_text"] = df["true_label"].map({
    0: "Non-cancer",
    1: "Cancer",
})

df["prediction_text"] = df["oof_prediction"].map({
    0: "Non-cancer",
    1: "Cancer",
})

df["correct"] = (
    df["true_label"]
    == df["oof_prediction"]
)

print("=" * 80)
print("OVARIAN CANCER RESEARCH DEMO")
print("=" * 80)
print()
print("Model: Age + 192 miRNAs")
print("Predictions: 5-fold out-of-fold")
print("Research use only — NOT a clinical diagnostic tool.")
print()

print("=" * 80)
print("CORRECTLY IDENTIFIED CANCER CASES")
print("=" * 80)

cases = df[
    (df["true_label"] == 1)
    & (df["oof_prediction"] == 1)
].sort_values(
    "oof_probability",
    ascending=False,
)

for _, r in cases.head(10).iterrows():
    print(
        f"{str(r['sample_id']):>8}  "
        f"cancer probability={r['oof_probability']:.3f}  "
        f"CORRECT"
    )

print()
print("=" * 80)
print("CORRECTLY IDENTIFIED NON-CANCER CASES")
print("=" * 80)

cases = df[
    (df["true_label"] == 0)
    & (df["oof_prediction"] == 0)
].sort_values(
    "oof_probability",
)

for _, r in cases.head(10).iterrows():
    print(
        f"{str(r['sample_id']):>8}  "
        f"cancer probability={r['oof_probability']:.3f}  "
        f"CORRECT"
    )

print()
print("=" * 80)
print("CASE LOOKUP")
print("=" * 80)
print()

lookup = {
    str(r["sample_id"]): r
    for _, r in df.iterrows()
}

while True:

    sid = input(
        "Enter sample ID (or q to quit): "
    ).strip()

    if sid.lower() == "q":
        break

    if sid not in lookup:
        print("Sample not found.")
        print()
        continue

    r = lookup[sid]

    print()
    print("-" * 60)

    print(
        f"Sample:            {sid}"
    )

    print(
        f"Predicted:         {r['prediction_text']}"
    )

    print(
        f"Cancer probability: "
        f"{r['oof_probability']:.1%}"
    )

    print(
        f"Actual outcome:    {r['actual_text']}"
    )

    print(
        "Correct:           "
        + ("YES" if r["correct"] else "NO")
    )

    print("-" * 60)
    print()


if __name__ == "__main__":
    pass
