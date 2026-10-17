#!/usr/bin/env python3
"""Streamlit research workbench for the frozen ovarian model."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
import streamlit as st

PROJECT = Path(__file__).resolve().parents[1]
ARTIFACT_DIR = PROJECT / "production_artifacts"
FEATURES = (
    "hsa-miR-200b-3p",
    "hsa-miR-296-3p",
    "hsa-miR-6826-5p",
    "hsa-miR-32-3p",
    "hsa-miR-1255a",
)
CLINICAL_FEATURES = ("age", "ca125")
DEFAULTS_PATH = ARTIFACT_DIR / "imputation_defaults.json"


st.set_page_config(
    page_title="Ovarian Research Workbench",
    page_icon=":microscope:",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <style>
    :root { --ink:#102a43; --muted:#627d98; --line:#d9e2ec; --teal:#087f8c; --teal-soft:#e6fffa; --red:#c53030; --red-soft:#fff5f5; --green:#276749; --green-soft:#f0fff4; }
    .stApp { background: linear-gradient(135deg, #f7fbfc 0%, #ffffff 55%, #eef7f7 100%); color: var(--ink); }
    [data-testid="stSidebar"] { background: #102a43; }
    [data-testid="stSidebar"] * { color: #f0f4f8; }
    [data-testid="stSidebar"] input { color: #102a43 !important; background: #ffffff !important; }
    [data-testid="stSidebar"] .stNumberInput label, [data-testid="stSidebar"] .stCheckbox label { color: #d9e2ec !important; }
    .portal-kicker { color: var(--teal); font-size: .78rem; font-weight: 700; letter-spacing: .12em; text-transform: uppercase; margin-bottom: .4rem; }
    .portal-title { color: var(--ink); font-size: 2.35rem; font-weight: 750; line-height: 1.08; margin: 0; }
    .portal-subtitle { color: var(--muted); font-size: 1rem; margin-top: .7rem; max-width: 760px; }
    .section-label { color: var(--ink); font-size: 1.05rem; font-weight: 700; border-bottom: 1px solid var(--line); padding-bottom: .55rem; margin: 1.35rem 0 .8rem; }
    .risk-card { border: 1px solid var(--line); border-radius: 12px; padding: 1.3rem 1.5rem; background: rgba(255,255,255,.9); box-shadow: 0 8px 24px rgba(16,42,67,.07); }
    .risk-label { color: var(--muted); font-size: .78rem; font-weight: 700; letter-spacing: .08em; text-transform: uppercase; }
    .risk-value { color: var(--ink); font-size: 2.65rem; font-weight: 800; margin: .2rem 0 .25rem; }
    .status-badge { border-radius: 8px; padding: .85rem 1rem; font-size: 1rem; font-weight: 750; margin: 1rem 0; }
    .status-low { background: var(--green-soft); border: 1px solid #9ae6b4; color: var(--green); }
    .status-high { background: var(--red-soft); border: 1px solid #feb2b2; color: var(--red); }
    .benchmark { background: #ffffff; border: 1px solid var(--line); border-radius: 10px; padding: 1rem; min-height: 100px; }
    .benchmark-value { color: var(--teal); font-size: 1.55rem; font-weight: 800; }
    .benchmark-label { color: var(--muted); font-size: .82rem; margin-top: .25rem; }
    .disclaimer { color: #486581; font-size: .8rem; border-left: 3px solid #9fb3c8; padding: .65rem .8rem; background: #f5f7fa; margin-top: 1.4rem; }
    div.stButton > button { background: var(--teal); color: white; border: 0; border-radius: 8px; font-weight: 750; min-height: 3rem; }
    div.stButton > button:hover { background: #05616b; color: white; }
    </style>
    """,
    unsafe_allow_html=True,
)


EXPECTED_FEATURES = FEATURES + CLINICAL_FEATURES + ("mirna_200b_ca125", "mirna_296_age")


def validate_artifact_schema(manifest, scaler=None, model=None) -> None:
    expected = list(EXPECTED_FEATURES)
    if not isinstance(manifest, dict) or manifest.get("features") != expected:
        raise ValueError("Artifact manifest must declare exact ordered native features")
    for artifact in (scaler, model):
        if artifact is not None:
            dimension = getattr(artifact, "n_features_in_", None)
            if dimension != len(expected):
                raise ValueError("Artifact feature dimension missing or incompatible")


    if model is not None:
        classes = getattr(model, "classes_", None)
        if classes is None or np.asarray(classes).shape != (2,) or not np.array_equal(classes, [0, 1]):
            raise ValueError("Model classes_ must be exactly [0, 1] in that order")


def validate_probabilities(probabilities):
    values = np.asarray(probabilities, dtype=float)
    if not np.isfinite(values).all() or ((values < 0) | (values > 1)).any():
        raise ValueError("Model probabilities must be finite and within [0, 1]")
    return values


def effective_clinical_inputs(values, defaults):
    effective, flags = {}, {}
    for feature in CLINICAL_FEATURES:
        default = float(defaults[feature])
        measured = float(values[feature])
        if not np.isfinite(default) or default < 0:
            raise ValueError(f"Invalid default for {feature}")
        if np.isinf(measured) or (not np.isnan(measured) and measured < 0):
            raise ValueError(f"Invalid measured {feature}")
        flags[feature] = bool(np.isnan(measured))
        effective[feature] = default if flags[feature] else measured
    return effective, flags


def make_snapshot(values, defaults, artifact_identity):
    effective, flags = effective_clinical_inputs(values, defaults)
    raw = {}
    for feature in FEATURES + CLINICAL_FEATURES:
        value = float(values[feature])
        if feature in FEATURES and (not np.isfinite(value) or value < 0):
            raise ValueError(f"Invalid molecular input: {feature}")
        raw[feature] = None if np.isnan(value) else value
    return json.loads(json.dumps({"raw_inputs": raw, "effective_inputs": {
        **{feature: raw[feature] for feature in FEATURES}, **effective},
        "default_flags": flags, "artifact_identity": artifact_identity}, allow_nan=False))


def snapshot_matches(stored, current) -> bool:
    return isinstance(stored, dict) and stored == current


def artifact_identity():
    return {name: hashlib.sha256((ARTIFACT_DIR / name).read_bytes()).hexdigest()
            for name in ("frozen_model.joblib", "feature_scaler.joblib",
                         "imputation_defaults.json", "manifest.json")}


def clear_result():
    st.session_state.pop("probability", None)
    st.session_state.pop("research_result", None)


def load_artifacts() -> tuple[Any, Any, dict[str, Any], dict[str, Any]]:
    manifest = json.loads((ARTIFACT_DIR / "manifest.json").read_text())
    validate_artifact_schema(manifest)
    model = joblib.load(ARTIFACT_DIR / "frozen_model.joblib")
    scaler = joblib.load(ARTIFACT_DIR / "feature_scaler.joblib")
    defaults = json.loads(DEFAULTS_PATH.read_text())
    validate_artifact_schema(manifest, scaler, model)
    effective_clinical_inputs({"age": np.nan, "ca125": np.nan}, defaults)
    return model, scaler, defaults, manifest


def build_model_input(values: dict[str, float], defaults: dict[str, Any], manifest: dict[str, Any]) -> pd.DataFrame:
    validate_artifact_schema(manifest)
    for feature in FEATURES:
        if not np.isfinite(values[feature]) or values[feature] < 0:
            raise ValueError(f"{feature} must be a finite non-negative value")
    clinical, _ = effective_clinical_inputs(values, defaults)
    age, ca125 = clinical["age"], clinical["ca125"]
    row: dict[str, float] = {
        feature: float(np.log1p(values[feature])) for feature in FEATURES
    }
    row["age"] = age
    row["ca125"] = ca125
    row["mirna_200b_ca125"] = row["hsa-miR-200b-3p"] * ca125
    row["mirna_296_age"] = row["hsa-miR-296-3p"] * age
    feature_order = list(EXPECTED_FEATURES)
    result = pd.DataFrame([[row[feature] for feature in feature_order]], columns=feature_order)
    if not np.isfinite(result.to_numpy(dtype=float)).all():
        raise ValueError("Non-finite transformed model input")
    return result


def normalized_label(value: Any) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).strip().lower())


def read_uploaded_table(uploaded_file: Any) -> pd.DataFrame:
    """Read either patient-per-row or feature-per-row clinical exports."""

    filename = str(uploaded_file.name).lower()
    if filename.endswith(".xlsx"):
        raw = pd.read_excel(uploaded_file, engine="openpyxl")
    else:
        raw = pd.read_csv(uploaded_file)
    if raw.empty:
        raise ValueError("The uploaded file contains no rows.")

    columns = [normalized_label(column) for column in raw.columns]
    if any(not column for column in columns) or len(set(columns)) != len(columns):
        raise ValueError("Blank or normalized duplicate upload columns")
    # Pandas may mangle repeated CSV headers; reject ambiguous mangled names too.
    if any(re.search(r"\.\d+$", str(column)) and str(column).rsplit(".", 1)[0] in raw.columns for column in raw.columns):
        raise ValueError("Duplicate upload columns")
    required = {normalized_label(feature): feature for feature in FEATURES}
    required.update({"age": "age", "ca125": "ca125"})
    column_matches = {
        required[normalized_label(column)]: column
        for column in raw.columns
        if normalized_label(column) in required
    }
    if all(feature in column_matches for feature in FEATURES):
        output = pd.DataFrame(index=raw.index)
        id_columns = {"patientid", "patient", "sampleid", "sample", "subjectid", "id"}
        id_column = next((column for column in raw.columns if normalized_label(column) in id_columns), None)
        output["Patient ID"] = raw[id_column] if id_column else [f"Patient {index + 1}" for index in range(len(raw))]
        for feature in FEATURES:
            output[feature] = pd.to_numeric(raw[column_matches[feature]], errors="raise")
        for clinical in CLINICAL_FEATURES:
            output[clinical] = pd.to_numeric(raw[column_matches[clinical]], errors="raise") if clinical in column_matches else np.nan
        return validate_uploaded_inputs(output).reset_index(drop=True)

    label_column = raw.columns[0]
    labels = [normalized_label(value) for value in raw[label_column]]
    if raw[label_column].isna().any() or any(not label for label in labels) or len(set(labels)) != len(labels):
        raise ValueError("Missing or normalized duplicate feature labels")
    label_lookup = {normalized_label(value): index for index, value in raw[label_column].items()}
    missing = [feature for feature in FEATURES if normalized_label(feature) not in label_lookup]
    if missing:
        raise ValueError("Could not find all Top-5 miRNA labels as columns or rows: " + ", ".join(missing))
    patient_columns = [column for column in raw.columns if column != label_column]
    if not patient_columns:
        raise ValueError("The uploaded feature-oriented file has no patient columns.")
    output = pd.DataFrame({"Patient ID": [str(column) for column in patient_columns]})
    for feature in FEATURES:
        output[feature] = pd.to_numeric(raw.loc[label_lookup[normalized_label(feature)], patient_columns], errors="raise").to_numpy()
    for clinical in CLINICAL_FEATURES:
        index = label_lookup.get(normalized_label(clinical))
        output[clinical] = pd.to_numeric(raw.loc[index, patient_columns], errors="raise").to_numpy() if index is not None else np.nan
    return validate_uploaded_inputs(output).reset_index(drop=True)


def validate_uploaded_inputs(table):
    ids = table["Patient ID"]
    normalized = ids.map(lambda value: str(value).strip())
    keys = normalized.map(normalized_label)
    if ids.isna().any() or (keys == "").any() or keys.duplicated().any():
        raise ValueError("Missing, blank or normalized duplicate patient IDs")
    table["Patient ID"] = normalized
    for feature in FEATURES + CLINICAL_FEATURES:
        values = table[feature].to_numpy(dtype=float)
        if np.isinf(values).any() or (values < 0).any():
            raise ValueError(f"Invalid input: {feature}")
        if feature in FEATURES and np.isnan(values).any():
            raise ValueError(f"Missing molecular input: {feature}")
    return table


def build_batch_model_input(table: pd.DataFrame, defaults: dict[str, Any], manifest: dict[str, Any]) -> pd.DataFrame:
    inputs = []
    for row in table.to_dict(orient="records"):
        values = {feature: float(row[feature]) if pd.notna(row[feature]) else np.nan for feature in FEATURES}
        values["age"] = float(row["age"]) if pd.notna(row["age"]) else np.nan
        values["ca125"] = float(row["ca125"]) if pd.notna(row["ca125"]) else np.nan
        inputs.append(build_model_input(values, defaults, manifest))
    if not inputs:
        raise ValueError("The uploaded file contains no patient records.")
    return pd.concat(inputs, ignore_index=True)


def render_bulk_upload(model: Any, scaler: Any, defaults: dict[str, Any], manifest: dict[str, Any], identity: dict[str, str]) -> None:
    st.markdown('<div class="section-label">Bulk Patient File Upload (.csv / .xlsx)</div>', unsafe_allow_html=True)
    st.write("Upload a patient-per-row or feature-per-row molecular pathology export. Age and CA-125 are optional and use frozen training medians when absent.")
    uploaded_file = st.file_uploader("Drop a CSV or Excel clinical export here", type=["csv", "xlsx"], key="bulk_patient_file")
    if uploaded_file is None:
        st.info("Awaiting a tabular patient export.")
        return
    try:
        table = read_uploaded_table(uploaded_file)
        validate_artifact_schema(manifest, scaler, model)
        if artifact_identity() != identity:
            raise ValueError("Artifacts changed since loading; reload before scoring")
        model_input = build_batch_model_input(table, defaults, manifest)
        probabilities = model.predict_proba(scaler.transform(model_input.to_numpy(dtype=float)))[:, 1]
        validate_probabilities(probabilities)
        if artifact_identity() != identity:
            raise ValueError("Artifacts changed during prediction; results discarded")
        results = pd.DataFrame({
            "Patient ID": table["Patient ID"].astype(str),
            "Research model probability %": np.round(probabilities * 100, 2),
            "Research classification": np.where(probabilities >= 0.5, "Above legacy research cutoff", "Below legacy research cutoff"),
        })
        for clinical in CLINICAL_FEATURES:
            results[f"{clinical} effective"] = table[clinical].fillna(float(defaults[clinical]))
            results[f"{clinical}_imputed"] = table[clinical].isna()
            results[f"{clinical} raw measured"] = table[clinical]
            results[f"{clinical} input source"] = np.where(table[clinical].isna(), "Frozen default", "Measured")
        results["research_only"] = True
        results["legacy_research_cutoff"] = 0.5
        results["threshold_selection_provenance"] = "Unknown"
        results["artifact_identity_sha256"] = json.dumps(identity, sort_keys=True)
        st.write({"artifact_identity_sha256": identity})
        st.caption("Legacy research cutoff: 0.5; selection provenance unknown. Research only.")
        st.write({"source": "production_artifacts/manifest.json",
                  "training_endpoint": manifest.get("endpoint", "Not supplied"),
                  "training_N": manifest.get("n_training_samples", "Not supplied"),
                  "positive_samples": manifest.get("positive_samples", "Not supplied"),
                  "negative_samples": manifest.get("negative_samples", "Not supplied"),
                  "source_reported_exclusions": {key: value for key, value in manifest.items()
                      if "exclu" in key.lower()}})
        if not any("exclu" in key.lower() for key in manifest):
            st.caption("Exclusions: unavailable in manifest.")
        st.caption("Linked evaluation unavailable; not validated. Source-reported endpoint and counts are not clinical validation.")
        st.caption("Prior five-miRNA-only external all-positive failure and historical 14-miRNA all-negative failure concern different pipelines; neither validates this portal. No platform cause established.")
        st.caption("Numbers require training-compatible assays; cross-platform measurements are not validated. No supported-assay certification is established.")
        st.warning("Research only. Not for diagnosis, screening, stage inference, or treatment decisions.")
        st.success(f"Processed {len(results)} patients using the frozen nine-feature model.")
        st.dataframe(results, use_container_width=True, hide_index=True)
        st.download_button("Download Results CSV", results.to_csv(index=False).encode("utf-8"), "ovarian_diagnostic_results.csv", "text/csv", use_container_width=True)
    except Exception as exc:
        clear_result()
        st.error(f"The research file could not be processed: {exc}")


def sidebar_inputs(defaults: dict[str, Any]) -> dict[str, float]:
    st.sidebar.markdown("## Research Input Workbench")
    st.sidebar.caption("Enter the pre-treatment clinical and liquid-biopsy measurements.")
    st.sidebar.markdown("### Clinical Indicators")
    age_pending = st.sidebar.checkbox("Mark Age as Pending/Missing", key="age_pending", value=True)
    age = st.sidebar.number_input("Patient Age (years)", min_value=0.0, max_value=130.0, value=float(defaults["age"]), step=0.1, disabled=age_pending)
    if age_pending:
        age = float("nan")
        st.sidebar.caption(f"Frozen default applied: {defaults['age']}")
    ca125_pending = st.sidebar.checkbox("Mark CA-125 as Pending/Missing", key="ca125_pending", value=True)
    ca125 = st.sidebar.number_input("CA-125 Test Level (U/mL)", min_value=0.0, max_value=100000.0, value=float(defaults["ca125"]), step=0.1, disabled=ca125_pending)
    if ca125_pending:
        ca125 = float("nan")
        st.sidebar.caption(f"Frozen default applied: {defaults['ca125']} U/mL")
    st.sidebar.markdown("### Liquid Biopsy Panel")
    values = {"age": age, "ca125": ca125}
    for feature in FEATURES:
        values[feature] = st.sidebar.number_input(feature, min_value=0.0, value=0.0, step=0.01, format="%.4f", key=feature)
    return values


def render_result(probability: float, manifest: dict[str, Any], identity=None) -> None:
    validate_probabilities(probability)
    st.markdown("### Research-only model output")
    st.metric("Research model probability (not clinical risk)", f"{probability:.2%}")
    st.write("Above legacy research cutoff" if probability >= 0.5 else "Below legacy research cutoff")
    st.write({"artifact_identity_sha256": identity if identity is not None else artifact_identity()})
    st.caption(f"Legacy research cutoff: 0.5 (selection provenance unknown); model class score max(p, 1-p): {max(probability, 1-probability):.2%}; not calibrated confidence.")
    st.write({"source": "production_artifacts/manifest.json",
              "training_endpoint": manifest.get("endpoint", "Not supplied"),
              "training_N": manifest.get("n_training_samples", "Not supplied"),
              "source_reported_counts_and_exclusions": {key: value for key, value in manifest.items()
                  if key.startswith("n_") or key in ("positive_samples", "negative_samples") or "exclu" in key.lower()}})
    if not any("exclu" in key.lower() for key in manifest):
        st.caption("Exclusions: unavailable in manifest.")
    st.caption("Linked evaluation unavailable; not validated. Source-reported endpoint and N are not clinical validation.")
    st.caption("Prior five-miRNA-only external all-positive failure and historical 14-miRNA all-negative failure concern different pipelines; neither validates this portal. No platform cause established.")
    st.caption("Numbers require training-compatible assays; cross-platform measurements are not validated. No supported-assay certification is established.")
    st.warning("Research only. No validated clinical benchmark. Not for diagnosis, screening, stage inference, or treatment decisions.")


def main() -> None:
    st.session_state.pop("probability", None)
    try:
        identity_before = artifact_identity()
        model, scaler, defaults, manifest = load_artifacts()
        identity = artifact_identity()
        if identity_before != identity:
            raise ValueError("Artifacts changed while loading")
    except Exception as exc:
        clear_result()
        st.error(f"Production artifacts could not be loaded: {exc}")
        st.stop()
        return
    try:
        values = sidebar_inputs(defaults)
        current = make_snapshot(values, defaults, identity)
        stored = st.session_state.get("research_result")
        if stored and not snapshot_matches(stored.get("snapshot"), current):
            clear_result()
    except Exception as exc:
        clear_result()
        st.error(f"Invalid inputs: {exc}")
        st.stop()
        return
    st.markdown('<div class="portal-kicker">Ovarian oncology | Frozen Phase 10 model</div><h1 class="portal-title">Ovarian Research Workbench</h1><p class="portal-subtitle">Research-only workspace. Not a clinical diagnostic or screening tool.</p>', unsafe_allow_html=True)
    st.markdown("---")
    mode = st.radio("Research input mode", ["Single patient", "Bulk upload"])
    previous_mode = st.session_state.get("research_mode")
    if previous_mode != mode:
        clear_result()
    st.session_state["research_mode"] = mode
    if mode == "Single patient":
        left, right = st.columns([1.35, 1], gap="large")
        with left:
            st.markdown('<div class="section-label">Ready for Research Execution</div>', unsafe_allow_html=True)
            st.write("Review the patient inputs in the research panel, then execute the frozen model when the specimen and clinical values are available.")
            measured = st.checkbox("I confirm all five miRNA values are entered specimen measurements, not placeholder zeros", key="molecular_measurement_attestation")
            execute = st.button("Execute Research Model", use_container_width=True, disabled=not measured)
            if not measured:
                clear_result()
            if execute and measured:
                try:
                    model_input = build_model_input(values, defaults, manifest)
                    probability = float(model.predict_proba(scaler.transform(model_input.to_numpy(dtype=float)))[:, 1][0])
                    if not np.isfinite(probability) or not 0 <= probability <= 1:
                        raise ValueError("Invalid model probability")
                    if artifact_identity() != identity:
                        raise ValueError("Artifacts changed during prediction; results discarded")
                    st.session_state["research_result"] = {"snapshot": current, "probability": probability}
                except Exception as exc:
                    clear_result()
                    st.error(f"Assessment could not be completed: {exc}")
        with right:
            st.markdown('<div class="section-label">Pipeline Status</div>', unsafe_allow_html=True)
            st.success(f"Frozen artifact set loaded\n\n{manifest.get('n_training_samples', 'Unknown')} curated training patients")
            st.caption("Features: Top-5 serum miRNAs + Age + CA-125 + two prespecified interactions")
    if mode == "Bulk upload":
        clear_result()
        render_bulk_upload(model, scaler, defaults, manifest, identity)
    result = st.session_state.get("research_result")
    if result and snapshot_matches(result.get("snapshot"), current):
        st.write({"effective_clinical_inputs": {key: current["effective_inputs"][key] for key in CLINICAL_FEATURES},
                  "default_flags": current["default_flags"]})
        render_result(float(result["probability"]), manifest, current["artifact_identity"])


if __name__ == "__main__":
    main()
