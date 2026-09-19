from pathlib import Path
import gzip

matrix = Path("data/external/gse9899/GSE9899_series_matrix.txt.gz")

print("=" * 100)
print("GSE9899 GEO SAMPLE METADATA")
print("=" * 100)

lines_found = 0

with gzip.open(matrix, "rt", encoding="utf-8", errors="replace") as f:

    for line in f:

        if line.startswith("!Sample_"):
            print("\n" + "-" * 100)
            print(line.rstrip())
            lines_found += 1

        # Stop once expression matrix begins.
        elif line.startswith("!series_matrix_table_begin"):
            break

print("\n" + "=" * 100)
print("TOTAL !Sample_ LINES:", lines_found)
print("=" * 100)
