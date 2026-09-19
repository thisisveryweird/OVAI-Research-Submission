from pathlib import Path
import requests
import zipfile

DATASET_ID = "th7fztbrv9"
VERSION = "11"

url = f"https://data.mendeley.com/public-api/zip/{DATASET_ID}/download/{VERSION}"

out_dir = Path("data/external/mendeley_ovarian")
out_dir.mkdir(parents=True, exist_ok=True)

zip_path = out_dir / "mendeley_ovarian_v11.zip"

print("=" * 80)
print("DOWNLOADING Mendeley ovarian dataset")
print("=" * 80)
print("URL:", url)

try:
    with requests.get(
        url,
        headers={"User-Agent": "Mozilla/5.0"},
        stream=True,
        timeout=120,
    ) as r:
        print("HTTP status:", r.status_code)
        print("Content type:", r.headers.get("content-type"))
        print("Content length:", r.headers.get("content-length"))

        r.raise_for_status()

        total = 0

        with open(zip_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
                    total += len(chunk)

        print("Downloaded bytes:", total)

    print("\nSaved ZIP:")
    print(zip_path)

    print("\n" + "=" * 80)
    print("CHECKING ZIP CONTENTS")
    print("=" * 80)

    with zipfile.ZipFile(zip_path, "r") as z:
        files = z.namelist()

        print("Number of files:", len(files))

        for name in files:
            info = z.getinfo(name)
            print(f"{info.file_size:>12,} bytes  {name}")

    print("\n" + "=" * 80)
    print("EXTRACTING")
    print("=" * 80)

    extract_dir = out_dir / "raw"
    extract_dir.mkdir(exist_ok=True)

    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(extract_dir)

    print("Extracted to:")
    print(extract_dir)

except Exception as e:
    print("\nERROR:", repr(e))

print("\nDONE")
