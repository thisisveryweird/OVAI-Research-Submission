from pathlib import Path
import requests
import re
import json

OUT = Path("data/external/clinical_404")
OUT.mkdir(parents=True, exist_ok=True)

queries = [
    '"Preoperative differentiation of borderline and malignant ovarian tumors using interpretable machine learning"',
    '"199" "205" borderline ovarian malignant 404 HE4 CA-125',
    '"199 BOT" "205 ovarian cancer" HE4 CA-125',
]

print("=" * 80)
print("SEARCHING FOR THE 404-PATIENT CLINICAL DATASET")
print("=" * 80)

for i, q in enumerate(queries, 1):
    print(f"\nQUERY {i}: {q}\n")

    # Crossref API
    url = "https://api.crossref.org/works"
    params = {
        "query.bibliographic": q,
        "rows": 5,
    }

    try:
        r = requests.get(url, params=params, timeout=20)
        r.raise_for_status()
        data = r.json()

        for j, item in enumerate(data["message"]["items"], 1):
            title = item.get("title", [""])[0]
            doi = item.get("DOI", "")
            published = item.get("published", {}).get("date-parts", [[""]])[0]
            print(f"{j}. {title}")
            print(f"   DOI: {doi}")
            print(f"   Date: {published}")
    except Exception as e:
        print("Crossref error:", repr(e))

print("\n" + "=" * 80)
print("SEARCHING ARTICLE PAGE")
print("=" * 80)

article_url = "https://link.springer.com/article/10.1186/s13048-026-02062-5"

try:
    r = requests.get(
        article_url,
        timeout=30,
        headers={"User-Agent": "Mozilla/5.0"}
    )
    print("HTTP status:", r.status_code)
    print("Final URL:", r.url)

    text = r.text

    # Print useful dataset-related phrases if they occur in page HTML
    patterns = [
        r"199.{0,150}205",
        r"205.{0,150}199",
        r"HE4",
        r"CA-125",
        r"Supplementary",
        r"Data availability",
        r"dataset",
        r"GitHub",
        r"Figshare",
        r"Zenodo",
        r"repository",
    ]

    print("\nRelevant snippets:\n")

    found = set()

    for pattern in patterns:
        for m in re.finditer(pattern, text, flags=re.I | re.S):
            start = max(0, m.start() - 250)
            end = min(len(text), m.end() + 500)
            snippet = re.sub(r"\s+", " ", text[start:end])

            key = snippet[:300]
            if key not in found:
                found.add(key)
                print(f"[{pattern}]")
                print(snippet[:1000])
                print("-" * 80)

except Exception as e:
    print("Article request error:", repr(e))

print("\n" + "=" * 80)
print("DONE")
print("=" * 80)
