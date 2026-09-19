from pathlib import Path
import requests

OUT = Path("data/external/gse9899")
OUT.mkdir(parents=True, exist_ok=True)

headers = {
    "User-Agent": "Mozilla/5.0"
}

candidates = [
    (
        "matrix",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE9nnn/GSE9899/matrix/GSE9899_series_matrix.txt.gz",
    ),
    (
        "clinical",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE9nnn/GSE9899/suppl/GSE9899_clinical_anns.csv",
    ),
]

print("=" * 80)
print("GSE9899 DOWNLOAD")
print("=" * 80)

session = requests.Session()
session.headers.update(headers)

for name, url in candidates:
    print("\n" + "-" * 80)
    print(name.upper())
    print("URL:", url)

    try:
        r = session.get(
            url,
            stream=True,
            timeout=120,
        )

        print("HTTP:", r.status_code)
        print("Content-Type:", r.headers.get("content-type"))
        print("Content-Length:", r.headers.get("content-length"))

        if r.status_code != 200:
            print("NOT DOWNLOADED")
            continue

        filename = url.split("/")[-1]
        path = OUT / filename

        total = 0

        with open(path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
                    total += len(chunk)

        print("Downloaded:", total, "bytes")
        print("Saved:", path)

    except Exception as e:
        print("ERROR:", repr(e))

print("\n" + "=" * 80)
print("FILES PRESENT")
print("=" * 80)

for p in sorted(OUT.iterdir()):
    print(f"{p.stat().st_size:>12,} bytes  {p.name}")

print("\nDONE")
