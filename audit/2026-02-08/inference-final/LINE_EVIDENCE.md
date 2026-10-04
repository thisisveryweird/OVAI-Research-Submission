# Numbered static source evidence
Original source line numbers; bounded spans reproduced from previously extracted evidence.

## scripts/ovarian_diagnostic_portal.py

```text
67: def load_artifacts() -> tuple[Any, Any, dict[str, Any], dict[str, Any]]:
68:     model = joblib.load(ARTIFACT_DIR / "frozen_model.joblib")
69:     scaler = joblib.load(ARTIFACT_DIR / "feature_scaler.joblib")
70:     defaults = json.loads(DEFAULTS_PATH.read_text())
71:     manifest = json.loads((ARTIFACT_DIR / "manifest.json").read_text())
72:     return model, scaler, defaults, manifest
73: 
74: 
75: def build_model_input(values: dict[str, float], defaults: dict[str, Any], manifest: dict[str, Any]) -> pd.DataFrame:
76:     for feature in FEATURES:
77:         if not np.isfinite(values[feature]) or values[feature] < 0:
78:             raise ValueError(f"{feature} must be a finite non-negative value")
79:     age = values["age"] if np.isfinite(values["age"]) else float(defaults["age"])
80:     ca125 = values["ca125"] if np.isfinite(values["ca125"]) else float(defaults["ca125"])
81:     row: dict[str, float] = {
82:         feature: float(np.log1p(values[feature])) for feature in FEATURES
83:     }
84:     row["age"] = age
85:     row["ca125"] = ca125
86:     row["mirna_200b_ca125"] = row["hsa-miR-200b-3p"] * ca125
87:     row["mirna_296_age"] = row["hsa-miR-296-3p"] * age
88:     feature_order = manifest.get("features", list(FEATURES) + ["age", "ca125", "mirna_200b_ca125", "mirna_296_age"])
89:     missing = [feature for feature in feature_order if feature not in row]
90:     if missing:
91:         raise ValueError(f"Artifact manifest requires unsupported features: {missing}")
92:     return pd.DataFrame([[row[feature] for feature in feature_order]], columns=feature_order)
93: 
```
```text
95: def normalized_label(value: Any) -> str:
96:     return re.sub(r"[^a-z0-9]", "", str(value).strip().lower())
97: 
98: 
99: def read_uploaded_table(uploaded_file: Any) -> pd.DataFrame:
100:     """Read either patient-per-row or feature-per-row clinical exports."""
101: 
102:     filename = str(uploaded_file.name).lower()
103:     if filename.endswith(".xlsx"):
104:         raw = pd.read_excel(uploaded_file, engine="openpyxl")
105:     else:
106:         raw = pd.read_csv(uploaded_file)
107:     if raw.empty:
108:         raise ValueError("The uploaded file contains no rows.")
109: 
110:     required = {normalized_label(feature): feature for feature in FEATURES}
111:     required.update({"age": "age", "ca125": "ca125"})
112:     column_matches = {
113:         required[normalized_label(column)]: column
114:         for column in raw.columns
115:         if normalized_label(column) in required
116:     }
117:     if all(feature in column_matches for feature in FEATURES):
118:         output = pd.DataFrame(index=raw.index)
119:         id_columns = {"patientid", "patient", "sampleid", "sample", "subjectid", "id"}
120:         id_column = next((column for column in raw.columns if normalized_label(column) in id_columns), None)
121:         output["Patient ID"] = raw[id_column].astype(str) if id_column else [f"Patient {index + 1}" for index in range(len(raw))]
122:         for feature in FEATURES:
123:             output[feature] = pd.to_numeric(raw[column_matches[feature]], errors="coerce")
124:         for clinical in CLINICAL_FEATURES:
125:             output[clinical] = pd.to_numeric(raw[column_matches[clinical]], errors="coerce") if clinical in column_matches else np.nan
126:         return output.reset_index(drop=True)
127: 
128:     label_column = raw.columns[0]
129:     label_lookup = {normalized_label(value): index for index, value in raw[label_column].items()}
130:     missing = [feature for feature in FEATURES if normalized_label(feature) not in label_lookup]
131:     if missing:
132:         raise ValueError("Could not find all Top-5 miRNA labels as columns or rows: " + ", ".join(missing))
133:     patient_columns = [column for column in raw.columns if column != label_column]
134:     if not patient_columns:
135:         raise ValueError("The uploaded feature-oriented file has no patient columns.")
136:     output = pd.DataFrame({"Patient ID": [str(column) for column in patient_columns]})
137:     for feature in FEATURES:
138:         output[feature] = pd.to_numeric(raw.loc[label_lookup[normalized_label(feature)], patient_columns], errors="coerce").to_numpy()
139:     for clinical in CLINICAL_FEATURES:
140:         index = label_lookup.get(normalized_label(clinical))
141:         output[clinical] = pd.to_numeric(raw.loc[index, patient_columns], errors="coerce").to_numpy() if index is not None else np.nan
142:     return output.reset_index(drop=True)
143: 
144: 
145: def build_batch_model_input(table: pd.DataFrame, defaults: dict[str, Any], manifest: dict[str, Any]) -> pd.DataFrame:
146:     inputs = []
147:     for row in table.to_dict(orient="records"):
148:         values = {feature: float(row[feature]) if pd.notna(row[feature]) else np.nan for feature in FEATURES}
149:         values["age"] = float(row["age"]) if pd.notna(row["age"]) else np.nan
150:         values["ca125"] = float(row["ca125"]) if pd.notna(row["ca125"]) else np.nan
151:         inputs.append(build_model_input(values, defaults, manifest))
152:     if not inputs:
153:         raise ValueError("The uploaded file contains no patient records.")
154:     return pd.concat(inputs, ignore_index=True)
155: 
156: 
157: def render_bulk_upload(model: Any, scaler: Any, defaults: dict[str, Any], manifest: dict[str, Any]) -> None:
158:     st.markdown('<div class="section-label">Bulk Patient File Upload (.csv / .xlsx)</div>', unsafe_allow_html=True)
159:     st.write("Upload a patient-per-row or feature-per-row molecular pathology export. Age and CA-125 are optional and use frozen training medians when absent.")
160:     uploaded_file = st.file_uploader("Drop a CSV or Excel clinical export here", type=["csv", "xlsx"], key="bulk_patient_file")
161:     if uploaded_file is None:
162:         st.info("Awaiting a tabular patient export.")
163:         return
164:     try:
165:         table = read_uploaded_table(uploaded_file)
166:         model_input = build_batch_model_input(table, defaults, manifest)
167:         probabilities = model.predict_proba(scaler.transform(model_input.to_numpy(dtype=float)))[:, 1]
168:         results = pd.DataFrame({
169:             "Patient ID": table["Patient ID"].astype(str),
170:             "Malignancy Probability %": np.round(probabilities * 100, 2),
171:             "Risk Status": np.where(probabilities >= 0.5, "CRITICAL: Elevated Risk", "LOW RISK: Healthy / Benign"),
172:         })
173:         st.success(f"Processed {len(results)} patients using the frozen nine-feature model.")
174:         st.dataframe(results, use_container_width=True, hide_index=True)
175:         st.download_button("Download Results CSV", results.to_csv(index=False).encode("utf-8"), "ovarian_diagnostic_results.csv", "text/csv", use_container_width=True)
176:     except Exception as exc:
177:         st.error(f"The clinical file could not be processed: {exc}")
178: 
```
```text
180: def sidebar_inputs(defaults: dict[str, Any]) -> dict[str, float]:
181:     st.sidebar.markdown("## Patient Diagnostic Workup")
182:     st.sidebar.caption("Enter the pre-treatment clinical and liquid-biopsy measurements.")
183:     st.sidebar.markdown("### Clinical Indicators")
184:     age_pending = st.sidebar.checkbox("Mark Age as Pending/Missing", key="age_pending")
185:     age = st.sidebar.number_input("Patient Age (years)", min_value=0.0, max_value=130.0, value=float(defaults["age"]), step=0.1, disabled=age_pending)
186:     if age_pending:
187:         age = float("nan")
188:         st.sidebar.caption(f"Frozen default applied: {defaults['age']}")
189:     ca125_pending = st.sidebar.checkbox("Mark CA-125 as Pending/Missing", key="ca125_pending")
190:     ca125 = st.sidebar.number_input("CA-125 Test Level (U/mL)", min_value=0.0, max_value=100000.0, value=float(defaults["ca125"]), step=0.1, disabled=ca125_pending)
191:     if ca125_pending:
192:         ca125 = float("nan")
193:         st.sidebar.caption(f"Frozen default applied: {defaults['ca125']} U/mL")
194:     st.sidebar.markdown("### Liquid Biopsy Panel")
195:     values = {"age": age, "ca125": ca125}
196:     for feature in FEATURES:
197:         values[feature] = st.sidebar.number_input(feature, min_value=0.0, value=0.0, step=0.01, format="%.4f", key=feature)
198:     return values
199: 
200: 
201: def render_result(probability: float, manifest: dict[str, Any]) -> None:
202:     st.markdown('<div class="section-label">Risk Assessment</div>', unsafe_allow_html=True)
203:     if probability >= 0.5:
204:         badge_class = "status-high"
205:         badge = "CRITICAL: Elevated Risk / Actionable Stage I/II Malignancy Warning"
206:     else:
207:         badge_class = "status-low"
208:         badge = "LOW RISK: Healthy / Benign Clinical Profile"
209:     st.markdown(
210:         f'<div class="risk-card"><div class="risk-label">Malignancy Probability</div><div class="risk-value">{probability:.2%}</div><div class="status-badge {badge_class}">{badge}</div><div class="risk-label">Model threshold: 50.0% | Training endpoint: {manifest.get("endpoint", "GSE94533 Stage I/II")}</div></div>',
211:         unsafe_allow_html=True,
212:     )
213:     st.markdown('<div class="section-label">Frozen Model Benchmarks</div>', unsafe_allow_html=True)
214:     columns = st.columns(3)
215:     benchmarks = [("86.43%", "Balanced Accuracy"), ("82.70%", "Early-Stage Sensitivity"), ("90.20%", "Surgical Specificity")]
216:     for column, (value, label) in zip(columns, benchmarks):
217:         with column:
218:             st.markdown(f'<div class="benchmark"><div class="benchmark-value">{value}</div><div class="benchmark-label">{label}</div></div>', unsafe_allow_html=True)
219:     st.markdown('<div class="disclaimer"><strong>Clinical use notice:</strong> This AI screening assistant evaluates molecular expression and clinical telemetry for Stage I/II invasive epithelial tumors. It is a decision-support aid and should be integrated with broader radiological imaging, pathology, examination findings, and specialist clinical judgment.</div>', unsafe_allow_html=True)
```
```text
222: def main() -> None:
223:     try:
224:         model, scaler, defaults, manifest = load_artifacts()
225:     except Exception as exc:
226:         st.error(f"Production artifacts could not be loaded: {exc}")
227:         st.stop()
228:     values = sidebar_inputs(defaults)
229:     st.markdown('<div class="portal-kicker">Ovarian oncology | Frozen Phase 10 model</div><h1 class="portal-title">Ovarian Diagnostic Portal</h1><p class="portal-subtitle">A focused clinical decision-support workspace for integrating serum microRNA expression with Age and CA-125 telemetry.</p>', unsafe_allow_html=True)
230:     st.markdown("---")
231:     manual_tab, bulk_tab = st.tabs(["Single Patient Assessment", "Bulk Patient File Upload (.csv / .xlsx)"])
232:     with manual_tab:
233:         left, right = st.columns([1.35, 1], gap="large")
234:         with left:
235:             st.markdown('<div class="section-label">Ready for Assessment</div>', unsafe_allow_html=True)
236:             st.write("Review the patient inputs in the workup panel, then execute the frozen model when the specimen and clinical values are available.")
237:             execute = st.button("Execute Malignancy Risk Assessment", use_container_width=True)
238:             if execute:
239:                 try:
240:                     model_input = build_model_input(values, defaults, manifest)
241:                     probability = float(model.predict_proba(scaler.transform(model_input.to_numpy(dtype=float)))[:, 1][0])
242:                     st.session_state["probability"] = probability
243:                 except Exception as exc:
244:                     st.error(f"Assessment could not be completed: {exc}")
245:         with right:
246:             st.markdown('<div class="section-label">Pipeline Status</div>', unsafe_allow_html=True)
247:             st.success(f"Frozen artifact set loaded\n\n{manifest.get('n_training_samples', 'Unknown')} curated training patients")
248:             st.caption("Features: Top-5 serum miRNAs + Age + CA-125 + two prespecified interactions")
249:     with bulk_tab:
250:         render_bulk_upload(model, scaler, defaults, manifest)
251:     if "probability" in st.session_state:
252:         render_result(float(st.session_state["probability"]), manifest)
253: 
254: 
255: if __name__ == "__main__":
256:     main()
```

## scripts/export_production_pipeline.py

```text
80: def transformed_values(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
81:     mirna = np.log1p(frame[list(FEATURES)].clip(lower=0).to_numpy(dtype=float))
82:     clinical = frame[list(CLINICAL_FEATURES)].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
83:     imputer = SimpleImputer(strategy="median")
84:     imputed = imputer.fit_transform(clinical)
85:     engineered = np.column_stack([
86:         mirna[:, 0] * imputed[:, 1],
87:         mirna[:, 1] * imputed[:, 0],
88:     ])
89:     return np.column_stack([mirna, imputed, engineered]), imputer.statistics_, clinical
90: 
```
```text
110:     scaled = scaler.transform(values)
111:     model = LogisticRegression(C=0.3, class_weight="balanced", solver="liblinear", max_iter=5000, random_state=SEED)
112:     model.fit(scaled, labels)
113:     args.output_dir.mkdir(parents=True, exist_ok=True)
114:     joblib.dump(model, args.output_dir / "frozen_model.joblib")
115:     joblib.dump(scaler, args.output_dir / "feature_scaler.joblib")
116:     (args.output_dir / "imputation_defaults.json").write_text(json.dumps({
117:         "age": float(medians[0]),
118:         "ca125": float(medians[1]),
119:         "clinical_features": list(CLINICAL_FEATURES),
120:         "policy": "median values fitted on the curated GSE94533 Stage I/II cancer versus benign/control training cohort",
121:     }, indent=2) + "\n")
122:     manifest = {
123:         "source": source,
124:         "endpoint": "GSE94533 Stage I/II invasive cancer versus benign/control",
125:         "n_training_samples": int(len(frame)),
126:         "positive_samples": int(labels.sum()),
127:         "negative_samples": int((labels == 0).sum()),
128:         "features": list(MODEL_FEATURES),
129:         "miRNA_transform": "log1p(max(value, 0))",
130:         "interaction_features": {
131:             "mirna_200b_ca125": "hsa-miR-200b-3p_log1p * ca125_imputed",
132:             "mirna_296_age": "hsa-miR-296-3p_log1p * age_imputed",
133:         },
134:         "clinical_transform": "median imputation using imputation_defaults.json",
135:         "scaler": "StandardScaler fitted after imputation",
136:         "model": "LogisticRegression(C=0.3, class_weight=balanced, solver=liblinear)",
137:         "random_state": SEED,
138:     }
139:     (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
140:     LOGGER.info("Exported production artifacts to %s", args.output_dir)
141:     LOGGER.info("Training samples=%d; positives=%d; negatives=%d", len(frame), int(labels.sum()), int((labels == 0).sum()))
142:     return 0
143: 
144: 
145: if __name__ == "__main__":
146:     raise SystemExit(main())
```

## scripts/predict_ovarian_risk.py

```text
26: def parse_number(value: str) -> float:
27:     if value.strip().lower() in {"", "none", "nan", "na", "null"}:
28:         return np.nan
29:     number = float(value)
30:     if not np.isfinite(number):
31:         raise ValueError(f"Non-finite input: {value}")
32:     return number
33: 
34: 
35: def main() -> int:
36:     parser = argparse.ArgumentParser(description=__doc__)
37:     for feature in FEATURES:
38:         parser.add_argument(f"--{feature}", required=True, type=parse_number)
39:     parser.add_argument("--age", required=True, type=parse_number)
40:     parser.add_argument("--ca125", required=True, type=parse_number)
41:     parser.add_argument("--artifact-dir", type=Path, default=DEFAULT_ARTIFACT_DIR)
42:     parser.add_argument("--high-confidence", type=float, default=0.80, help="Confidence threshold for the high-confidence flag.")
43:     args = parser.parse_args()
44:     if not 0.5 <= args.high_confidence < 1.0:
45:         raise ValueError("--high-confidence must be at least 0.5 and below 1.0")
46:     artifact_dir = args.artifact_dir
47:     model = joblib.load(artifact_dir / "frozen_model.joblib")
48:     scaler = joblib.load(artifact_dir / "feature_scaler.joblib")
49:     defaults = json.loads((artifact_dir / "imputation_defaults.json").read_text())
50:     row = {feature: getattr(args, feature.replace("-", "_")) for feature in FEATURES}
51:     row.update({"age": args.age, "ca125": args.ca125})
52:     values = pd.DataFrame([row], columns=list(FEATURES) + list(CLINICAL_FEATURES))
53:     values.loc[:, FEATURES] = np.log1p(values.loc[:, FEATURES].clip(lower=0))
54:     for feature in CLINICAL_FEATURES:
55:         values[feature] = values[feature].fillna(float(defaults[feature]))
56:     probability = float(model.predict_proba(scaler.transform(values.to_numpy(dtype=float)))[:, 1][0])
57:     classification = "Cancer" if probability >= 0.5 else "Healthy/Benign"
58:     confidence = max(probability, 1.0 - probability)
59:     result = {
60:         "cancer_probability": probability,
61:         "classification": classification,
62:         "confidence": confidence,
63:         "high_confidence": bool(confidence >= args.high_confidence),
64:         "threshold": 0.5,
65:         "high_confidence_threshold": args.high_confidence,
66:     }
67:     print(json.dumps(result, indent=2))
68:     return 0
69: 
```

## scripts/ovarian_research_demo.py

```text
4: ROOT = Path.home() / "ovarian-ai"
5: 
6: PRED = (
7:     ROOT
8:     / "reports"
9:     / "tables"
10:     / "age_plus_192mirna_oof_predictions.csv"
11: )
12: 
13: df = pd.read_csv(PRED)
14: 
15: df["actual_text"] = df["true_label"].map({
16:     0: "Non-cancer",
17:     1: "Cancer",
18: })
19: 
20: df["prediction_text"] = df["oof_prediction"].map({
21:     0: "Non-cancer",
22:     1: "Cancer",
23: })
24: 
25: df["correct"] = (
26:     df["true_label"]
27:     == df["oof_prediction"]
28: )
29: 
30: print("=" * 80)
31: print("OVARIAN CANCER RESEARCH DEMO")
32: print("=" * 80)
33: print()
34: print("Model: Age + 192 miRNAs")
35: print("Predictions: 5-fold out-of-fold")
36: print("Research use only — NOT a clinical diagnostic tool.")
```
```text
43: cases = df[
44:     (df["true_label"] == 1)
45:     & (df["oof_prediction"] == 1)
46: ].sort_values(
47:     "oof_probability",
48:     ascending=False,
49: )
50: 
51: for _, r in cases.head(10).iterrows():
52:     print(
53:         f"{str(r['sample_id']):>8}  "
54:         f"cancer probability={r['oof_probability']:.3f}  "
55:         f"CORRECT"
```
```text
82: 
83: lookup = {
84:     str(r["sample_id"]): r
85:     for _, r in df.iterrows()
```

## scripts/demo_ovarian_cases.py

```text
7: ROOT = Path.home() / "ovarian-ai"
8: 
9: OOF = (
10:     ROOT
11:     / "reports"
12:     / "tables"
13:     / "fair_historical_14mirna_oof_predictions.csv"
14: )
15: 
16: 
17: def main():
18: 
19:     print("=" * 80)
20:     print("OVARIAN CANCER miRNA RESEARCH DEMO")
21:     print("=" * 80)
22:     print()
23:     print("Research/validation demonstration only.")
24:     print("NOT a clinical diagnostic tool.")
25:     print()
26: 
27:     if not OOF.exists():
28:         raise SystemExit(
29:             f"Missing OOF prediction file:\n{OOF}"
30:         )
31: 
32:     df = pd.read_csv(OOF)
33: 
34:     required = {
35:         "sample_id",
36:         "true_label",
37:         "oof_probability",
38:         "oof_prediction",
39:     }
40: 
41:     missing = required - set(df.columns)
42: 
43:     if missing:
44:         raise SystemExit(
45:             "Missing columns: "
46:             + ", ".join(sorted(missing))
```
```text
52: 
53:     if len(df) != 179:
54:         raise SystemExit(
55:             f"Expected 179 patients, found {len(df)}."
56:         )
```
```text
154: 
155:     by_id = {
156:         str(x): row
157:         for _, row in df.iterrows()
158:     }
```

## scripts/test_high_capacity_models.py

```text
56: def load_training_data() -> tuple[pd.DataFrame, np.ndarray]:
57:     if CLINICAL_DATA.exists():
58:         frame = pd.read_csv(CLINICAL_DATA)
59:         source = CLINICAL_DATA
60:     else:
61:         if not FALLBACK_REPORT.exists() or not BASE_MATRIX.exists():
62:             raise FileNotFoundError(f"Missing clinical fusion and fallback inputs: {CLINICAL_DATA}, {FALLBACK_REPORT}, {BASE_MATRIX}")
63:         clinical = pd.read_csv(FALLBACK_REPORT).drop_duplicates("sample_id_normalized")
64:         matrix = pd.read_csv(BASE_MATRIX, usecols=["sample_id_normalized", *SIGNATURE])
65:         frame = clinical.merge(matrix, on="sample_id_normalized", how="inner", validate="one_to_one")
66:         source = FALLBACK_REPORT
67:         LOGGER.warning("Canonical clinical fusion absent; using explicit derived source %s", source)
68:     frame = validate_endpoint(frame)
69:     values = frame[list(BASE_FEATURES)].copy()
70:     values[list(SIGNATURE)] = np.log1p(values[list(SIGNATURE)].clip(lower=0))
71:     return values, frame["binary_label"].astype(int).to_numpy()
```
```text
85: def base_pipeline(estimator: Any) -> Pipeline:
86:     return Pipeline([
87:         ("impute", SimpleImputer(strategy="median")),
88:         ("scale", StandardScaler()),
89:         ("model", estimator),
90:     ])
91: 
92: 
93: def engineered_frame(values: pd.DataFrame) -> pd.DataFrame:
94:     engineered = values.copy()
95:     engineered["mirna_200b_ca125"] = engineered["hsa-miR-200b-3p"] * engineered["ca125"]
96:     engineered["mirna_296_age"] = engineered["hsa-miR-296-3p"] * engineered["age"]
97:     return engineered
98: 
```
```text
117: def run(values: pd.DataFrame, labels: np.ndarray) -> pd.DataFrame:
118:     factories = model_factories()
119:     outer = StratifiedKFold(5, shuffle=True, random_state=SEED)
120:     predictions = {name: np.full(len(labels), np.nan) for name in factories}
121:     for fold, (train_index, validation_index) in enumerate(outer.split(values, labels), start=1):
122:         LOGGER.info("Outer fold %d/5", fold)
123:         train = values.iloc[train_index]
124:         validation = values.iloc[validation_index]
125:         engineered_train = engineered_frame(train)
126:         engineered_validation = engineered_frame(validation)
127:         for name, factory in factories.items():
128:             if name == "engineered_logistic":
129:                 estimator = factory()
130:                 estimator.fit(engineered_train, labels[train_index])
131:                 probability = estimator.predict_proba(engineered_validation)[:, 1]
132:             else:
133:                 estimator = factory()
134:                 estimator.fit(train, labels[train_index])
135:                 probability = estimator.predict_proba(validation)[:, 1]
136:             predictions[name][validation_index] = probability
137:     rows = []
138:     for name, probability in predictions.items():
139:         result = metrics(labels, probability)
140:         result["model"] = name
141:         result["exceeds_90_auc"] = bool(result["roc_auc"] > 0.90)
142:         rows.append(result)
143:     return pd.DataFrame(rows)[["model", "roc_auc", "balanced_accuracy", "sensitivity", "specificity", "exceeds_90_auc"]]
144: 
```

## Launch_Diagnostic_Portal.bat

```text
1: @echo off
2: setlocal
3: cd /d "%~dp0"
4: python -m pip install -r requirements.txt >nul 2>&1
5: if errorlevel 1 exit /b 1
6: python -m streamlit run scripts\ovarian_diagnostic_portal.py --browser.gatherUsageStats=false
```

## Launch_Diagnostic_Portal.command

```text
1: #!/bin/bash
2: set -e
3: cd "$(dirname "$0")"
4: python3 -m pip install -r requirements.txt >/dev/null 2>&1
5: exec python3 -m streamlit run scripts/ovarian_diagnostic_portal.py --browser.gatherUsageStats=false
```

## production_artifacts/imputation_defaults.json

```text
1: {
2:   "age": 55.5,
3:   "ca125": 87.2,
4:   "clinical_features": [
5:     "age",
6:     "ca125"
7:   ],
8:   "policy": "median values fitted on the curated GSE94533 Stage I/II cancer versus benign/control training cohort"
9: }
```

## production_artifacts/manifest.json

```text
1: {
2:   "source": "reports/error_analysis_predictions.csv",
3:   "endpoint": "GSE94533 Stage I/II invasive cancer versus benign/control",
4:   "n_training_samples": 113,
5:   "positive_samples": 52,
6:   "negative_samples": 61,
7:   "features": [
8:     "hsa-miR-200b-3p",
9:     "hsa-miR-296-3p",
10:     "hsa-miR-6826-5p",
11:     "hsa-miR-32-3p",
12:     "hsa-miR-1255a",
13:     "age",
14:     "ca125",
15:     "mirna_200b_ca125",
16:     "mirna_296_age"
17:   ],
18:   "miRNA_transform": "log1p(max(value, 0))",
19:   "interaction_features": {
20:     "mirna_200b_ca125": "hsa-miR-200b-3p_log1p * ca125_imputed",
21:     "mirna_296_age": "hsa-miR-296-3p_log1p * age_imputed"
22:   },
23:   "clinical_transform": "median imputation using imputation_defaults.json",
24:   "scaler": "StandardScaler fitted after imputation",
25:   "model": "LogisticRegression(C=0.3, class_weight=balanced, solver=liblinear)",
26:   "random_state": 20260919
27: }
```