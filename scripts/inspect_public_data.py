#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
RAW = PROJECT / "data" / "raw"
EXTERNAL = PROJECT / "data" / "external"
REPORT = PROJECT / "reports" / "tables"
REPORT.mkdir(parents=True, exist_ok=True)


def summarize_excel(path: Path) -> dict[str, Any]:
    result: dict[str, Any] = {"path": str(path), "sheets": {}}
    book = pd.ExcelFile(path)
    for sheet in book.sheet_names:
        df = pd.read_excel(path, sheet_name=sheet)
        info: dict[str, Any] = {
            "rows": int(df.shape[0]),
            "columns": int(df.shape[1]),
            "column_names_first_30": [str(x) for x in df.columns[:30]],
            "dtypes": {str(k): str(v) for k, v in df.dtypes.items()},
            "missing_fraction_top_30": {
                str(k): float(v)
                for k, v in df.isna().mean().sort_values(ascending=False).head(30).items()
            },
        }
        # Look for likely phenotype columns without assuming their exact names.
        likely = []
        for col in df.columns:
            text = str(col).lower()
            if any(k in text for k in ["diagn", "case", "status", "group", "cancer", "benign", "control", "class", "stage"]):
                likely.append(str(col))
        info["likely_label_or_metadata_columns"] = likely
        result["sheets"][sheet] = info
    return result


def main() -> None:
    paths = [
        *sorted(RAW.glob("GSE94533_*.xlsx")),
        EXTERNAL / "MethylBERT_EOC_DataS1_mmc2.xlsx",
    ]
    all_results = []
    for path in paths:
        if not path.exists():
            continue
        print(f"\n=== {path.name} ===")
        result = summarize_excel(path)
        print(json.dumps(result, indent=2)[:12000])
        all_results.append(result)

    out = REPORT / "dataset_audit.json"
    out.write_text(json.dumps(all_results, indent=2))
    print(f"\nAudit written to {out}")


if __name__ == "__main__":
    main()
