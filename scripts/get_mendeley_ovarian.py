from pathlib import Path
import requests
import re
from urllib.parse import urljoin

DATASET_URL = "https://data.mendeley.com/datasets/th7fztbrv9/11"

out_dir = Path("data/external/mendeley_ovarian")
out_dir.mkdir(parents=True, exist_ok=True)

headers = {
    "User-Agent": "Mozilla/5.0"
}

print("=" * 80)
print("Mendeley ovarian dataset audit")
print("=" * 80)
print("Dataset:", DATASET_URL)

try:
    r = requests.get(DATASET_URL, headers=headers, timeout=30)
    r.raise_for_status()

    print("HTTP status:", r.status_code)
    print("Final URL:", r.url)
    print("Downloaded page:", len(r.text), "bytes")

    html = r.text

    # Save raw HTML so we can inspect it if necessary.
    html_path = out_dir / "mendeley_page.html"
    html_path.write_text(html, encoding="utf-8")
    print("Saved:", html_path)

    print("\nPossible downloadable files/links:")
    print("-" * 80)

    links = re.findall(r'href=["\']([^"\']+)["\']', html, flags=re.I)

    seen = set()

    for link in links:
        full = urljoin(r.url, link)

        # Focus on likely file/download links.
        low = full.lower()

        if any(x in low for x in [
            "download",
            ".csv",
            ".xlsx",
            ".xls",
            ".zip",
            ".txt",
            ".json",
            ".sav",
            ".mat",
            "public-files",
            "api.mendeley"
        ]):
            if full not in seen:
                seen.add(full)
                print(full)

    print("\nPossible file names mentioned in page source:")
    print("-" * 80)

    patterns = [
        r'[^"\'<>\s]+\.(?:csv|xlsx|xls|zip|txt|json|sav|mat)',
        r'public-files/[^"\'<>\s]+',
        r'download[^"\'<>\s]+',
    ]

    file_mentions = set()

    for pattern in patterns:
        for match in re.findall(pattern, html, flags=re.I):
            file_mentions.add(match)

    for item in sorted(file_mentions):
        print(item)

    print("\nDataset metadata clues:")
    print("-" * 80)

    for keyword in [
        "Machine Learning for Predicting Ovarian Cancer",
        "version",
        "files",
        "created",
        "modified",
        "license",
        "CC BY",
        "349",
        "borderline",
        "cancer",
    ]:
        positions = [m.start() for m in re.finditer(re.escape(keyword), html, flags=re.I)]

        if positions:
            print(f"\n[{keyword}] found {len(positions)} time(s)")
            for pos in positions[:3]:
                snippet = re.sub(r"\s+", " ", html[max(0, pos-300):pos+700])
                print(snippet[:1000])
                print()

except Exception as e:
    print("ERROR:", repr(e))

print("=" * 80)
print("DONE")
print("=" * 80)
