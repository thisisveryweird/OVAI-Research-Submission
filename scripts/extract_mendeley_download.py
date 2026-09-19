from pathlib import Path
import re
import requests

url = "https://data.mendeley.com/datasets/th7fztbrv9/11"

headers = {
    "User-Agent": "Mozilla/5.0"
}

out = Path("data/external/mendeley_ovarian")
out.mkdir(parents=True, exist_ok=True)

r = requests.get(url, headers=headers, timeout=30)
r.raise_for_status()

html = r.text

print("=" * 80)
print("SEARCHING Mendeley page source for dataset download endpoints")
print("=" * 80)

# Save page
(out / "mendeley_page.html").write_text(html, encoding="utf-8")

patterns = [
    r'https?://[^"\']+',
    r'/api/[^"\']+',
    r'/public-api/[^"\']+',
    r'/datasets/[^"\']+',
    r'/download[^"\']+',
    r'files[^"\']{0,300}',
    r'fileId[^"\']{0,300}',
    r'file_id[^"\']{0,300}',
    r'publicFiles[^"\']{0,300}',
]

seen = set()

for pattern in patterns:
    matches = re.findall(pattern, html, flags=re.I)

    for m in matches:
        m = m.replace("\\/", "/")

        # Ignore enormous/minified junk strings.
        if len(m) > 1000:
            continue

        if m not in seen:
            seen.add(m)

            low = m.lower()

            if any(key in low for key in [
                "download",
                "dataset",
                "file",
                "api",
                "public",
            ]):
                print("\n", m)

print("\n" + "=" * 80)
print("SEARCHING FOR NEXT.JS / REACT DATA OBJECTS")
print("=" * 80)

# Look for useful JSON-ish objects embedded in page.
keywords = [
    "datasetId",
    "dataset_id",
    "versionId",
    "version_id",
    "fileList",
    "filesList",
    "files",
    "currentFile",
    "downloadAll",
    "download-all",
]

for keyword in keywords:
    positions = [m.start() for m in re.finditer(keyword, html, flags=re.I)]

    if positions:
        print(f"\n### {keyword}: {len(positions)} occurrence(s)")

        for pos in positions[:5]:
            start = max(0, pos - 500)
            end = min(len(html), pos + 1500)

            snippet = re.sub(r"\s+", " ", html[start:end])
            print(snippet)
            print("-" * 80)

print("\nSaved page to:")
print(out / "mendeley_page.html")
print("\nDONE")
