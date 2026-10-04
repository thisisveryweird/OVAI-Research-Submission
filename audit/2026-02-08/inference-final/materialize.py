import pathlib,json,hashlib,csv
root=pathlib.Path(__file__).parent
es=json.loads((root.parent/'inference-review/evidence.json').read_text())
spans={'scripts/ovarian_diagnostic_portal.py':[(67,93),(95,178),(180,219),(222,256)],'scripts/export_production_pipeline.py':[(80,90),(110,155)],'scripts/predict_ovarian_risk.py':[(26,69)],'scripts/ovarian_research_demo.py':[(4,36),(43,55),(82,85)],'scripts/demo_ovarian_cases.py':[(7,46),(52,56),(154,158)],'scripts/test_high_capacity_models.py':[(56,71),(85,98),(117,144)],'Launch_Diagnostic_Portal.bat':[(1,6)],'Launch_Diagnostic_Portal.command':[(1,5)],'production_artifacts/manifest.json':[(1,27)],'production_artifacts/imputation_defaults.json':[(1,9)]}
out=['# Numbered static source evidence','Original source line numbers; bounded spans reproduced from previously extracted evidence.']
rows=[]
art=['# Current artifact consistency','Static byte hashes only. No model loading or pickle/joblib deserialization.']
for e in es:
    p=e['path'].replace('\\','/')
    if p.startswith('audit/'): continue
    data=pathlib.Path(p).read_bytes()
    h=hashlib.sha256(data).hexdigest()
    same=h==e['sha256']
    rows.append([p,e['mode'],len(e.get('lines',[])) or '', 'all extracted source lines read' if 'lines' in e else 'prior binary string hits read only','none; current bytes hashed only',h,same,'; '.join(f'{a}-{b}' for a,b in spans.get(p,[]))])
    if p in spans:
        out.append('\n## '+p+'\n')
        for a,b in spans[p]:
            out.append('```text\n'+'\n'.join(f'{i}: {e["lines"][i-1]}' for i in range(a,min(b,len(e['lines']))+1))+'\n```')
    if p.startswith('production_artifacts/'):
        art.append(f'\n- `{p}`: {len(data)} bytes; SHA-256 `{h}`; matches prior extraction: {same}')
(root/'LINE_EVIDENCE.md').write_text('\n'.join(out),encoding='utf-8')
with (root/'COVERAGE.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.writer(f)
    w.writerow(['path','prior_extraction_mode','extracted_source_line_count','actual_review','direct_source_read','current_sha256','matches_prior_hash','materialized_deep_spans'])
    w.writerows(rows)
art.append('''
## Order, class and dimensions
Manifest lines 7–17 specify nine features, in order: hsa-miR-200b-3p, hsa-miR-296-3p, hsa-miR-6826-5p, hsa-miR-32-3p, hsa-miR-1255a, age, ca125, mirna_200b_ca125, mirna_296_age. Defaults lines 2–6 specify age=55.5, ca125=87.2 and clinical order age, ca125.

Prior safe binary string evidence identifies `sklearn.preprocessing._data`, `StandardScaler` and `n_features_in_` in the scaler; `sklearn.linear_model._logistic` and `n_features_in_` in the model. These strings suggest serialized estimator types, not authenticated object schemas. Actual serialized n_features_in_ values, classes_ order, coefficient/intercept dimensions and scaler array dimensions are **unknown in this final audit**: the provided evidence contains string offsets but no decoded dimensions/classes. Other workers' notes are not substituted for inspected static evidence. No dynamic compatibility claim is made.

Export source `main` fits StandardScaler and balanced liblinear LogisticRegression(C=0.3) on nine constructed values; expected binary labels are 0/1 and endpoint composition is 52 positive, 61 negative. This is source/manifest intent, not proof of serialized contents. Portal and CLI both assume predict_proba column 1 is cancer without checking classes_. Portal follows manifest order; CLI supplies seven features and does not read it. Hash agreement authenticates equality to prior extraction, not scientific validity, environment reproducibility or training provenance. Manifest source is reports/error_analysis_predictions.csv, not the export script's default curated path. Manifest lacks dependency versions and validation provenance; binary contents were not decoded to recover these.
''')
(root/'ARTIFACT_CONSISTENCY.md').write_text('\n'.join(art),encoding='utf-8')
