# Current artifact consistency
Static byte hashes only. No model loading or pickle/joblib deserialization.

- `production_artifacts/imputation_defaults.json`: 213 bytes; SHA-256 `ee2d67c92ebe224b1829d230234d16ff6decf117c76b11d85d16ffa083e287c6`; matches prior extraction: True

- `production_artifacts/manifest.json`: 879 bytes; SHA-256 `fa411c200093dca0306faac4e54df60c049e7f07d7dca13f508aaf60bfa4c1e6`; matches prior extraction: True

- `production_artifacts/feature_scaler.joblib`: 783 bytes; SHA-256 `5a63f905e2744238cc376ac2aaa776de8d5bcac8574146f2322354d4267c623f`; matches prior extraction: True

- `production_artifacts/frozen_model.joblib`: 971 bytes; SHA-256 `0a3698e97e083d723566127c44abddf671f6cd8722e9b2509ed926e76ddd1989`; matches prior extraction: True

## Order, class and dimensions
Manifest lines 7–17 specify nine features, in order: hsa-miR-200b-3p, hsa-miR-296-3p, hsa-miR-6826-5p, hsa-miR-32-3p, hsa-miR-1255a, age, ca125, mirna_200b_ca125, mirna_296_age. Defaults lines 2–6 specify age=55.5, ca125=87.2 and clinical order age, ca125.

Prior safe binary string evidence identifies `sklearn.preprocessing._data`, `StandardScaler` and `n_features_in_` in the scaler; `sklearn.linear_model._logistic` and `n_features_in_` in the model. These strings suggest serialized estimator types, not authenticated object schemas. Actual serialized n_features_in_ values, classes_ order, coefficient/intercept dimensions and scaler array dimensions are **unknown in this final audit**: the provided evidence contains string offsets but no decoded dimensions/classes. Other workers' notes are not substituted for inspected static evidence. No dynamic compatibility claim is made.

Export source `main` fits StandardScaler and balanced liblinear LogisticRegression(C=0.3) on nine constructed values; expected binary labels are 0/1 and endpoint composition is 52 positive, 61 negative. This is source/manifest intent, not proof of serialized contents. Portal and CLI both assume predict_proba column 1 is cancer without checking classes_. Portal follows manifest order; CLI supplies seven features and does not read it. Hash agreement authenticates equality to prior extraction, not scientific validity, environment reproducibility or training provenance. Manifest source is reports/error_analysis_predictions.csv, not the export script's default curated path. Manifest lacks dependency versions and validation provenance; binary contents were not decoded to recover these.
