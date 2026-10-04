#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


ROOT = Path.home() / "ovarian-ai"

OOF = (
    ROOT
    / "reports"
    / "tables"
    / "fair_historical_14mirna_oof_predictions.csv"
)


def build_case_lookup(df) -> dict:
    if "sample_id" not in df.columns:
        raise ValueError("Missing sample_id column")
    lookup = {}
    for _, row in df.iterrows():
        value = row["sample_id"]
        if pd.isna(value) or not str(value).strip():
            raise ValueError("Missing or blank sample_id")
        key = str(value).strip()
        if key in lookup:
            raise ValueError(f"Duplicate sample_id: {key}")
        lookup[key] = row.copy()
    return lookup


def main():

    print("=" * 80)
    print("OVARIAN CANCER miRNA RESEARCH DEMO")
    print("=" * 80)
    print()
    print("Research/validation demonstration only.")
    print("NOT a clinical diagnostic tool.")
    print()

    if not OOF.exists():
        raise SystemExit(
            f"Missing OOF prediction file:\n{OOF}"
        )

    df = pd.read_csv(OOF)
    build_case_lookup(df)

    required = {
        "sample_id",
        "true_label",
        "oof_probability",
        "oof_prediction",
    }

    missing = required - set(df.columns)

    if missing:
        raise SystemExit(
            "Missing columns: "
            + ", ".join(sorted(missing))
        )

    # ---------------------------------------------------------------
    # Basic validation
    # ---------------------------------------------------------------

    if len(df) != 179:
        raise SystemExit(
            f"Expected 179 patients, found {len(df)}."
        )

    df["actual"] = df["true_label"].map(
        {
            0: "Non-cancer",
            1: "Cancer",
        }
    )

    df["predicted"] = df["oof_prediction"].map(
        {
            0: "Non-cancer",
            1: "Cancer",
        }
    )

    df["correct"] = (
        df["true_label"]
        == df["oof_prediction"]
    )

    print(
        f"Patients: {len(df)}"
    )

    print(
        f"Correct OOF predictions: "
        f"{int(df['correct'].sum())}/{len(df)}"
    )

    print()

    # ---------------------------------------------------------------
    # Correctly identified cancers
    # ---------------------------------------------------------------

    cancers = df[
        (df["true_label"] == 1)
        & (df["oof_prediction"] == 1)
    ].copy()

    cancers = cancers.sort_values(
        "oof_probability",
        ascending=False,
    )

    print("=" * 80)
    print("CORRECTLY IDENTIFIED CANCER CASES")
    print("=" * 80)

    print(
        cancers[
            [
                "sample_id",
                "oof_probability",
                "actual",
                "predicted",
            ]
        ]
        .head(10)
        .to_string(index=False)
    )

    # ---------------------------------------------------------------
    # Correctly identified non-cancers
    # ---------------------------------------------------------------

    noncancers = df[
        (df["true_label"] == 0)
        & (df["oof_prediction"] == 0)
    ].copy()

    noncancers = noncancers.sort_values(
        "oof_probability",
        ascending=True,
    )

    print()
    print("=" * 80)
    print("CORRECTLY IDENTIFIED NON-CANCER CASES")
    print("=" * 80)

    print(
        noncancers[
            [
                "sample_id",
                "oof_probability",
                "actual",
                "predicted",
            ]
        ]
        .head(10)
        .to_string(index=False)
    )

    # ---------------------------------------------------------------
    # Interactive case lookup
    # ---------------------------------------------------------------

    by_id = build_case_lookup(df)

    print()
    print("=" * 80)
    print("INTERACTIVE CASE LOOKUP")
    print("=" * 80)

    print(
        "Enter a sample ID such as 50217 or e1056."
    )
    print(
        "Enter 'q' to quit."
    )
    print()

    while True:

        sid = input(
            "Sample ID > "
        ).strip()

        if sid.lower() == "q":
            break

        key = sid

        if key not in by_id:
            print(
                "Sample ID not found."
            )
            print()
            continue

        row = by_id[key]

        probability = float(
            row["oof_probability"]
        )

        prediction = str(
            row["predicted"]
        )

        actual = str(
            row["actual"]
        )

        correct = bool(
            row["correct"]
        )

        print()
        print(
            f"Sample:      {sid}"
        )

        print(
            f"Predicted:   {prediction}"
        )

        print(
            f"Cancer score: {probability:.4f}"
        )

        print(
            f"Actual:      {actual}"
        )

        print(
            "Correct:     "
            + ("YES" if correct else "NO")
        )

        print()


if __name__ == "__main__":
    main()
