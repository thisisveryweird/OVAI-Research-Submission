from pathlib import Path
import pandas as pd

ROOT = Path("data/external/mendeley_ovarian/raw")

files = sorted(ROOT.rglob("*.xlsx"))

print("=" * 100)
print("MENDELEY OVARIAN DATASET — SPREADSHEET AUDIT")
print("=" * 100)

if not files:
    print("ERROR: No XLSX files found.")
    raise SystemExit(1)

for path in files:
    print("\n" + "=" * 100)
    print("FILE:", path.name)
    print("=" * 100)

    try:
        xls = pd.ExcelFile(path)
        print("SHEETS:", xls.sheet_names)

        for sheet in xls.sheet_names:
            print("\n" + "-" * 80)
            print("SHEET:", sheet)
            print("-" * 80)

            df = pd.read_excel(path, sheet_name=sheet)

            print("SHAPE:", df.shape)
            print("COLUMNS:")
            for i, col in enumerate(df.columns):
                print(f"  [{i}] {col!r}")

            print("\nFIRST 5 ROWS:")
            print(df.head().to_string())

            print("\nDTYPES:")
            print(df.dtypes.to_string())

            # Look for likely diagnostic / clinical variables
            keywords = [
                "diagn",
                "cancer",
                "border",
                "benign",
                "stage",
                "figo",
                "he4",
                "ca125",
                "ca-125",
                "ca 125",
                "age",
                "menop",
                "tumor",
                "hist",
                "path",
                "malig",
                "lmp",
            ]

            print("\nLIKELY CLINICAL / LABEL COLUMNS:")
            found = False

            for col in df.columns:
                col_text = str(col).lower()

                if any(k in col_text for k in keywords):
                    print(f"  {col!r}")
                    found = True

            if not found:
                print("  None found by name.")

            # Show unique values for likely categorical columns
            print("\nPOSSIBLE LABEL/CATEGORICAL VALUES:")

            for col in df.columns:
                s = df[col]

                if (
                    s.dtype == "object"
                    or str(s.dtype).startswith("category")
                ):
                    nunique = s.nunique(dropna=True)

                    if 2 <= nunique <= 30:
                        vals = s.dropna().astype(str).value_counts()

                        print(f"\n  {col!r} ({nunique} unique):")
                        print(vals.head(30).to_string())

            print()

    except Exception as e:
        print("ERROR READING FILE:", repr(e))

print("\n" + "=" * 100)
print("AUDIT COMPLETE")
print("=" * 100)
