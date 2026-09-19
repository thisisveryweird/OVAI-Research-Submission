from pathlib import Path
import requests
import re

BASE = "https://data.mendeley.com"
BUNDLE = BASE + "/datasets/bundle.js?342865b28b070babe69b"

out = Path("data/external/mendeley_ovarian")
out.mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("DOWNLOADING Mendeley dataset JavaScript")
print("=" * 80)

r = requests.get(
    BUNDLE,
    headers={"User-Agent": "Mozilla/5.0"},
    timeout=60,
)

print("HTTP:", r.status_code)
print("Bytes:", len(r.content))

r.raise_for_status()

js = r.text

bundle_path = out / "bundle.js"
bundle_path.write_text(js, encoding="utf-8")

print("Saved:", bundle_path)

print("\n" + "=" * 80)
print("SEARCHING FOR FILE/DOWNLOAD API ROUTES")
print("=" * 80)

# Print strings containing likely API/file terms.
patterns = [
    r'["\']([^"\']{0,250}/api/[^"\']{0,250})["\']',
    r'["\']([^"\']{0,250}/files[^"\']{0,250})["\']',
    r'["\']([^"\']{0,250}/download[^"\']{0,250})["\']',
    r'["\']([^"\']{0,250}datasets-v2[^"\']{0,250})["\']',
    r'["\']([^"\']{0,250}download-all[^"\']{0,250})["\']',
]

seen = set()

for pattern in patterns:
    for match in re.findall(pattern, js, flags=re.I):
        match = match.replace("\\/", "/")

        if len(match) > 500:
            continue

        if match not in seen:
            seen.add(match)
            print("\n", match)

print("\n" + "=" * 80)
print("SEARCHING FOR 'files' REQUEST CODE")
print("=" * 80)

keywords = [
    "files",
    "downloadFile",
    "downloadAll",
    "dataset-files",
    "public-files",
    "file-tree",
    "currentFile",
    "filesLoading",
    "datasets-v2",
]

for keyword in keywords:
    positions = [m.start() for m in re.finditer(re.escape(keyword), js, flags=re.I)]

    print(f"\n### {keyword}: {len(positions)} matches")

    for pos in positions[:8]:
        start = max(0, pos - 600)
        end = min(len(js), pos + 1400)

        snippet = re.sub(r"\s+", " ", js[start:end])

        print(snippet[:2000])
        print("-" * 60)

print("\nDONE")
