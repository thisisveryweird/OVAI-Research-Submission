#!/usr/bin/env python3
"""Streamlit doctor-facing portal for the frozen ovarian diagnostic model."""

from __future__ import annotations

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
    page_title="Ovarian Diagnostic Portal",
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


def load_artifacts() -> tuple[Any, Any, dict[str, Any], dict[str, Any]]:
    model = joblib.load(ARTIFACT_DIR / "frozen_model.joblib")
    scaler = joblib.load(ARTIFACT_DIR / "feature_scaler.joblib")
    defaults = json.loads(DEFAULTS_PATH.read_text())
    manifest = json.loads((ARTIFACT_DIR / "manifest.json").read_text())
    return model, scaler, defaults, manifest


def build_model_input(values: dict[str, float], defaults: dict[str, Any], manifest: dict[str, Any]) -> pd.DataFrame:
    for feature in FEATURES:
        if not np.isfinite(values[feature]) or values[feature] < 0:
            raise ValueError(f"{feature} must be a finite non-negative value")
    age = values["age"] if np.isfinite(values["age"]) else float(defaults["age"])
    ca125 = values["ca125"] if np.isfinite(values["ca125"]) else float(defaults["ca125"])
    row: dict[str, float] = {
        feature: float(np.log1p(values[feature])) for feature in FEATURES
    }
    row["age"] = age
    row["ca125"] = ca125
    row["mirna_200b_ca125"] = row["hsa-miR-200b-3p"] * ca125
    row["mirna_296_age"] = row["hsa-miR-296-3p"] * age
    feature_order = manifest.get("features", list(FEATURES) + ["age", "ca125", "mirna_200b_ca125", "mirna_296_age"])
    missing = [feature for feature in feature_order if feature not in row]
    if missing:
        raise ValueError(f"Artifact manifest requires unsupported features: {missing}")
    return pd.DataFrame([[row[feature] for feature in feature_order]], columns=feature_order)


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
        output["Patient ID"] = raw[id_column].astype(str) if id_column else [f"Patient {index + 1}" for index in range(len(raw))]
        for feature in FEATURES:
            output[feature] = pd.to_numeric(raw[column_matches[feature]], errors="coerce")
        for clinical in CLINICAL_FEATURES:
            output[clinical] = pd.to_numeric(raw[column_matches[clinical]], errors="coerce") if clinical in column_matches else np.nan
        return output.reset_index(drop=True)

    label_column = raw.columns[0]
    label_lookup = {normalized_label(value): index for index, value in raw[label_column].items()}
    missing = [feature for feature in FEATURES if normalized_label(feature) not in label_lookup]
    if missing:
        raise ValueError("Could not find all Top-5 miRNA labels as columns or rows: " + ", ".join(missing))
    patient_columns = [column for column in raw.columns if column != label_column]
    if not patient_columns:
        raise ValueError("The uploaded feature-oriented file has no patient columns.")
    output = pd.DataFrame({"Patient ID": [str(column) for column in patient_columns]})
    for feature in FEATURES:
        output[feature] = pd.to_numeric(raw.loc[label_lookup[normalized_label(feature)], patient_columns], errors="coerce").to_numpy()
    for clinical in CLINICAL_FEATURES:
        index = label_lookup.get(normalized_label(clinical))
        output[clinical] = pd.to_numeric(raw.loc[index, patient_columns], errors="coerce").to_numpy() if index is not None else np.nan
    return output.reset_index(drop=True)


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


def render_bulk_upload(model: Any, scaler: Any, defaults: dict[str, Any], manifest: dict[str, Any]) -> None:
    st.markdown('<div class="section-label">Bulk Patient File Upload (.csv / .xlsx)</div>', unsafe_allow_html=True)
    st.write("Upload a patient-per-row or feature-per-row molecular pathology export. Age and CA-125 are optional and use frozen training medians when absent.")
    uploaded_file = st.file_uploader("Drop a CSV or Excel clinical export here", type=["csv", "xlsx"], key="bulk_patient_file")
    if uploaded_file is None:
        st.info("Awaiting a tabular patient export.")
        return
    try:
        table = read_uploaded_table(uploaded_file)
        model_input = build_batch_model_input(table, defaults, manifest)
        probabilities = model.predict_proba(scaler.transform(model_input.to_numpy(dtype=float)))[:, 1]
        results = pd.DataFrame({
            "Patient ID": table["Patient ID"].astype(str),
            "Malignancy Probability %": np.round(probabilities * 100, 2),
            "Risk Status": np.where(probabilities >= 0.5, "CRITICAL: Elevated Risk", "LOW RISK: Healthy / Benign"),
        })
        st.success(f"Processed {len(results)} patients using the frozen nine-feature model.")
        st.dataframe(results, use_container_width=True, hide_index=True)
        st.download_button("Download Results CSV", results.to_csv(index=False).encode("utf-8"), "ovarian_diagnostic_results.csv", "text/csv", use_container_width=True)
    except Exception as exc:
        st.error(f"The clinical file could not be processed: {exc}")


def sidebar_inputs(defaults: dict[str, Any]) -> dict[str, float]:
    st.sidebar.markdown("## Patient Diagnostic Workup")
    st.sidebar.caption("Enter the pre-treatment clinical and liquid-biopsy measurements.")
    st.sidebar.markdown("### Clinical Indicators")
    age_pending = st.sidebar.checkbox("Mark Age as Pending/Missing", key="age_pending")
    age = st.sidebar.number_input("Patient Age (years)", min_value=0.0, max_value=130.0, value=float(defaults["age"]), step=0.1, disabled=age_pending)
    if age_pending:
        age = float("nan")
        st.sidebar.caption(f"Frozen default applied: {defaults['age']}")
    ca125_pending = st.sidebar.checkbox("Mark CA-125 as Pending/Missing", key="ca125_pending")
    ca125 = st.sidebar.number_input("CA-125 Test Level (U/mL)", min_value=0.0, max_value=100000.0, value=float(defaults["ca125"]), step=0.1, disabled=ca125_pending)
    if ca125_pending:
        ca125 = float("nan")
        st.sidebar.caption(f"Frozen default applied: {defaults['ca125']} U/mL")
    st.sidebar.markdown("### Liquid Biopsy Panel")
    values = {"age": age, "ca125": ca125}
    for feature in FEATURES:
        values[feature] = st.sidebar.number_input(feature, min_value=0.0, value=0.0, step=0.01, format="%.4f", key=feature)
    return values


def render_result(probability: float, manifest: dict[str, Any]) -> None:
    st.markdown('<div class="section-label">Risk Assessment</div>', unsafe_allow_html=True)
    if probability >= 0.5:
        badge_class = "status-high"
        badge = "CRITICAL: Elevated Risk / Actionable Stage I/II Malignancy Warning"
    else:
        badge_class = "status-low"
        badge = "LOW RISK: Healthy / Benign Clinical Profile"
    st.markdown(
        f'<div class="risk-card"><div class="risk-label">Malignancy Probability</div><div class="risk-value">{probability:.2%}</div><div class="status-badge {badge_class}">{badge}</div><div class="risk-label">Model threshold: 50.0% | Training endpoint: {manifest.get("endpoint", "GSE94533 Stage I/II")}</div></div>',
        unsafe_allow_html=True,
    )
    st.markdown('<div class="section-label">Frozen Model Benchmarks</div>', unsafe_allow_html=True)
    columns = st.columns(3)
    benchmarks = [("86.43%", "Balanced Accuracy"), ("82.70%", "Early-Stage Sensitivity"), ("90.20%", "Surgical Specificity")]
    for column, (value, label) in zip(columns, benchmarks):
        with column:
            st.markdown(f'<div class="benchmark"><div class="benchmark-value">{value}</div><div class="benchmark-label">{label}</div></div>', unsafe_allow_html=True)
    st.markdown('<div class="disclaimer"><strong>Clinical use notice:</strong> This AI screening assistant evaluates molecular expression and clinical telemetry for Stage I/II invasive epithelial tumors. It is a decision-support aid and should be integrated with broader radiological imaging, pathology, examination findings, and specialist clinical judgment.</div>', unsafe_allow_html=True)


def main() -> None:
    try:
        model, scaler, defaults, manifest = load_artifacts()
    except Exception as exc:
        st.error(f"Production artifacts could not be loaded: {exc}")
        st.stop()
    values = sidebar_inputs(defaults)
    st.markdown('<div class="portal-kicker">Ovarian oncology | Frozen Phase 10 model</div><h1 class="portal-title">Ovarian Diagnostic Portal</h1><p class="portal-subtitle">A focused clinical decision-support workspace for integrating serum microRNA expression with Age and CA-125 telemetry.</p>', unsafe_allow_html=True)
    st.markdown("---")
    manual_tab, bulk_tab = st.tabs(["Single Patient Assessment", "Bulk Patient File Upload (.csv / .xlsx)"])
    with manual_tab:
        left, right = st.columns([1.35, 1], gap="large")
        with left:
            st.markdown('<div class="section-label">Ready for Assessment</div>', unsafe_allow_html=True)
            st.write("Review the patient inputs in the workup panel, then execute the frozen model when the specimen and clinical values are available.")
            execute = st.button("Execute Malignancy Risk Assessment", use_container_width=True)
            if execute:
                try:
                    model_input = build_model_input(values, defaults, manifest)
                    probability = float(model.predict_proba(scaler.transform(model_input.to_numpy(dtype=float)))[:, 1][0])
                    st.session_state["probability"] = probability
                except Exception as exc:
                    st.error(f"Assessment could not be completed: {exc}")
        with right:
            st.markdown('<div class="section-label">Pipeline Status</div>', unsafe_allow_html=True)
            st.success(f"Frozen artifact set loaded\n\n{manifest.get('n_training_samples', 'Unknown')} curated training patients")
            st.caption("Features: Top-5 serum miRNAs + Age + CA-125 + two prespecified interactions")
    with bulk_tab:
        render_bulk_upload(model, scaler, defaults, manifest)
    if "probability" in st.session_state:
        render_result(float(st.session_state["probability"]), manifest)


if __name__ == "__main__":
    main()
