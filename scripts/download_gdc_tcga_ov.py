#!/usr/bin/env python3

"""Create a reproducible GDC inventory for TCGA-OV and optionally download data.

The inventory is the source of truth for the raw cohort. The GDC repository
contains multiple releases and file formats for some modalities, so this
script deliberately records file metadata before any downstream parsing.
"""

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import requests
from tqdm import tqdm


PROJECT = "TCGA-OV"
API = "https://api.gdc.cancer.gov"
PAGE_SIZE = 1_000

MODALITY_RULES = {
    "mutations": ("simple nucleotide variation", "masked somatic mutation"),
    "mrna": ("gene expression quantification",),
    "mirna": ("mirna expression quantification",),
    "cnv": ("copy number segment", "copy number variation"),
    "methylation": ("methylation beta value", "methylation"),
    "clinical": ("clinical",),
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/external/gdc/tcga_ov"),
        help="Directory for the inventory, case map, and optional downloads.",
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download open-access files listed in the inventory.",
    )
    parser.add_argument(
        "--modality",
        choices=sorted(MODALITY_RULES),
        action="append",
        help="Limit downloads to one or more modalities; repeat this option.",
    )
    return parser.parse_args()


def gdc_get(session, path, params):
    response = session.get(f"{API}{path}", params=params, timeout=120)
    response.raise_for_status()
    return response.json()


def project_filter():
    return json.dumps({
        "op": "in",
        "content": {
            "field": "cases.project.project_id",
            "value": [PROJECT],
        },
    })


def classify_file(record):
    values = " ".join(
        str(record.get(field) or "")
        for field in ("data_category", "data_type", "experimental_strategy")
    ).lower()
    for modality, keywords in MODALITY_RULES.items():
        if any(keyword in values for keyword in keywords):
            return modality
    return "other"


def fetch_files(session):
    fields = ",".join([
        "file_id",
        "file_name",
        "file_size",
        "md5sum",
        "data_category",
        "data_type",
        "experimental_strategy",
        "access",
        "state",
        "analysis.workflow_type",
        "cases.case_id",
        "cases.submitter_id",
        "associated_entities.entity_submitter_id",
    ])
    records = []
    offset = 0
    while True:
        payload = gdc_get(
            session,
            "/files",
            {
                "filters": project_filter(),
                "format": "JSON",
                "fields": fields,
                "size": PAGE_SIZE,
                "from": offset,
            },
        )
        hits = payload["data"]["hits"]
        records.extend(hits)
        if len(hits) < PAGE_SIZE:
            break
        offset += len(hits)
    return records


def first_case_value(record, field):
    cases = record.get("cases") or []
    if cases:
        return cases[0].get(field, "")
    return ""


def normalize_file(record):
    workflow = record.get("analysis") or {}
    if isinstance(workflow, list):
        workflow = workflow[0] if workflow else {}
    return {
        "file_id": record.get("file_id", ""),
        "file_name": record.get("file_name", ""),
        "file_size": record.get("file_size", ""),
        "md5sum": record.get("md5sum", ""),
        "data_category": record.get("data_category", ""),
        "data_type": record.get("data_type", ""),
        "experimental_strategy": record.get("experimental_strategy", ""),
        "access": record.get("access", ""),
        "state": record.get("state", ""),
        "workflow_type": workflow.get("workflow_type", ""),
        "case_id": first_case_value(record, "case_id"),
        "submitter_id": first_case_value(record, "submitter_id"),
        "modality": classify_file(record),
    }


def fetch_cases(session):
    fields = ",".join([
        "case_id",
        "submitter_id",
        "project.project_id",
        "demographic",
        "diagnoses",
        "sample_ids",
    ])
    payload = gdc_get(
        session,
        "/cases",
        {
            "filters": project_filter(),
            "format": "JSON",
            "fields": fields,
            "size": PAGE_SIZE,
        },
    )
    return payload["data"]["hits"]


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def write_inventory(path, records):
    columns = list(records[0]) if records else [
        "file_id", "file_name", "file_size", "md5sum", "data_category",
        "data_type", "experimental_strategy", "access", "state",
        "workflow_type", "case_id", "submitter_id", "modality",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def download_file(session, record, root):
    destination = root / record["modality"] / record["file_id"] / record["file_name"]
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and record["md5sum"]:
        digest = hashlib.md5(destination.read_bytes()).hexdigest()
        if digest == record["md5sum"]:
            return "verified"
    with session.get(
        f"{API}/data/{record['file_id']}",
        stream=True,
        timeout=300,
    ) as response:
        response.raise_for_status()
        with destination.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
    if record["md5sum"]:
        digest = hashlib.md5(destination.read_bytes()).hexdigest()
        if digest != record["md5sum"]:
            destination.unlink()
            raise RuntimeError(f"Checksum mismatch for {record['file_id']}")
    return "downloaded"


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers["Accept"] = "application/json"

    print(f"Querying GDC file inventory for {PROJECT}...")
    files = [normalize_file(record) for record in fetch_files(session)]
    cases = fetch_cases(session)
    write_inventory(args.output_dir / "files.tsv", files)
    write_json(args.output_dir / "cases.json", cases)
    write_json(args.output_dir / "query.json", {
        "project": PROJECT,
        "api": API,
        "file_count": len(files),
        "case_count": len(cases),
        "modalities": sorted({record["modality"] for record in files}),
    })

    print(f"Cases: {len(cases)}")
    print(f"Files: {len(files)}")
    print(f"Inventory: {args.output_dir / 'files.tsv'}")
    print(f"Case metadata: {args.output_dir / 'cases.json'}")

    if not args.download:
        print("Metadata only. Re-run with --download to fetch open-access files.")
        return

    selected = set(args.modality or MODALITY_RULES)
    downloads = [
        record for record in files
        if record["access"] == "open"
        and record["modality"] in selected
    ]
    print(f"Downloading {len(downloads)} open-access files...")
    for record in tqdm(downloads, unit="file"):
        download_file(session, record, args.output_dir / "raw")


if __name__ == "__main__":
    try:
        main()
    except requests.RequestException as error:
        print(f"GDC request failed: {error}", file=sys.stderr)
        raise SystemExit(1)