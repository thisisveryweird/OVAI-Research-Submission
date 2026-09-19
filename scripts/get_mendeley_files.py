import requests
import json

DATASET_ID = "th7fztbrv9"
VERSION = "11"

urls = [
    f"https://data.mendeley.com/public-api/datasets/{DATASET_ID}/versions/{VERSION}/files",
    f"https://data.mendeley.com/public-api/datasets/{DATASET_ID}/{VERSION}/files",
    f"https://api.mendeley.com/datasets/{DATASET_ID}/{VERSION}/files",
]

headers = {
    "User-Agent": "Mozilla/5.0",
    "Accept": "application/json",
}

for url in urls:
    print("=" * 80)
    print("TRYING:", url)
    print("=" * 80)

    try:
        r = requests.get(url, headers=headers, timeout=30)

        print("STATUS:", r.status_code)
        print("CONTENT TYPE:", r.headers.get("content-type"))
        print("FINAL URL:", r.url)

        text = r.text
        print("RESPONSE LENGTH:", len(text))

        if r.status_code == 200:
            try:
                data = r.json()

                print("\nJSON RESPONSE:")
                print(json.dumps(data, indent=2)[:30000])

                # Save successful response
                with open(
                    "data/external/mendeley_ovarian/files_api_response.json",
                    "w",
                    encoding="utf-8",
                ) as f:
                    json.dump(data, f, indent=2)

                print("\nSUCCESS: saved API response.")
                break

            except Exception:
                print("\nResponse was not valid JSON.")
                print(text[:10000])

        else:
            print("\nERROR RESPONSE:")
            print(text[:5000])

    except Exception as e:
        print("REQUEST ERROR:", repr(e))

print("\nDONE")
