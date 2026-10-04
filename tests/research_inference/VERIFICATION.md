# Offline containment verification

## Final commands and results
Run from repository root, with bytecode disabled:
- `python -B -m unittest discover -s tests/research_inference -v`: 24 tests, 15 passed, 9 dependency skips, no failures/errors, exit 0. Default Python 3.14 / NumPy 2.4.6; pandas unavailable. Full output: `test_results.txt`.
- `py -3.11 -B -m unittest discover -s tests/research_inference -v`: 24 tests, 23 passed, 1 genuine XLSX dependency skip, no failures/errors, exit 0. Installed Python 3.11.9 / NumPy 1.26.4 / pandas 2.3.3; openpyxl absent. Full output: `test_results_py311.txt`.
- Interpreter discovery: `py -0p`; `py -3.11 -B -c "import importlib.util,sys; print(sys.executable); print({m:bool(importlib.util.find_spec(m)) for m in ['numpy','pandas','openpyxl']})"` discovered the usable existing Python 3.11 installation. No install/environment creation/network.
- `git diff --name-only`: exactly the two Launch_Diagnostic_Portal launchers and scripts/{demo_ovarian_cases,ovarian_diagnostic_portal,predict_ovarian_risk,train_mirna_baseline}. No tracked edits beyond those six containment paths. This worker writes only this tests directory. Inventory includes all six test/document/result files.

## Coverage
`test_containment.py` exercises actual AST-extracted helpers and main functions, not whole product imports:
- Exact native ordered CLI7/portal9 manifests; missing/reversed/wrong/duplicated features; fake scaler/model dimension rejection; model missing/reversed/wrong classes rejection; scalar/vector probabilities finite and [0,1].
- CLI nine-feature manifest fails before mocked joblib; valid seven-feature hyphenated arguments reach a DataFrame sentinel with the exact row. Separate real-pandas synthetic CLI transform/default/prediction test, neutral research output, nonfinite/out-of-range probability rejection. Fake paths return synthetic JSON; fake joblib returns synthetic objects only.
- Baseline main blocked before data IO and help exits 0; no top-level mkdir execution. In-memory AST compilation of four Python files; existence/static noninstall/telemetry checks for both launchers (launchers are not Python and are not AST-compiled).
- Clinical NaN defaulting; invalid defaults and measured negative/infinite/non-numeric rejection; detached JSON-safe snapshots; NaN stability; changed raw values/defaults/artifact invalidation; clear both result keys. Actual portal main with fake Streamlit exercises artifact-loading failure, changed input snapshot, artifact replacement/load race, mode change, missing measurement attestation, and failed prediction clearing. No real UI is started.
- Real pandas CSV patient-row/feature-row parsing, duplicated headers/labels/IDs, malformed/negative/infinite/missing molecular rejection, single/batch input parity, expected log1p/interactions/default values; demo invalid IDs.
- SHA-256-only reads of every production_artifacts entry in `audit/resume-verification/integrity/initial_snapshot.json`: all four match. No artifact deserialization and no audit writes.

## Skips and boundaries
The XLSX both-orientation test is materialized but genuinely skipped because installed openpyxl is absent. No XLSX integration claim. Default-Python additional skips are real pandas dependency skips; the installed Python 3.11 run verifies those CSV/numeric tests. Fake DataFrame sentinel is only argument-regression coverage, not pandas emulation.

`ast_helpers.py` excludes imports, decorators, top-level calls and UI side effects; only selected functions and safe literal/simple tuple-concatenation constants execute with synthetic globals. No patient files, real model loads, training, server/UI, installation or bytecode writes. Hash reading immutable artifacts is the sole intentional product-artifact read. Source fixes remain owned by other workers; the initial old test fixture lacked classes_ after source hardening and was corrected to supply [0,1], with explicit rejected-class tests added rather than weakening assertions.

Containment verification is not scientific validation. Full UI rendering, real artifact inference, XLSX parsing, calibration, clinical performance, broader scientific pipeline and deployment remain unverified/unfinished. Only today's minimal six-path containment work is covered.

## Final scope checklist and blocked follow-up
- Added real-pandas BytesIO CSV/fake-UI bulk success, pre-score identity mismatch, and post-predict identity race regression. Captured frame/download verify raw/effective/default source, both imputation flags, research-only, cutoff 0.5/provenance Unknown, and exactly the supplied loaded synthetic identity. Both races clear manual/stale score and suppress success/frame/download; pre-score race never transforms/predicts.
- Final git diff --check clean; name-only exactly six authorized paths listed above; staged name-only empty. All 356 initial snapshot hashes compared: exactly those six changed, other 350 unchanged (including source PDF, reports/audit and all four production artifacts). Prior resume-verification inventory exactly 21 files; all 20 closure-manifest recorded hashes match. Untracked files outside snapshot/prior audit/approved test tree: none. Existing untracked audit and source PDF are expected.
- Implemented fail-closed baseline, NOT operational training: original baseline gate executes before IO, preventing new training and outputs. CLI intentionally rejects current nine-feature artifact. Portal/evaluation transformation mismatch unchanged. No calibration choice or new threshold.
- Global filters, nested selection, threshold provenance, calibration, external independence, full-provenance logging, governance and broader roadmap remain unresolved; original comprehensive roadmap is NOT complete.
- XLSX, integration, real-artifact inference and real UI unverified. No launcher execution, deserialization, server, training, installation or network; no clinical claim.
- Initial added test had an import-time newline-quoting SyntaxError, reported immediately and corrected in the test literal only. Final exact run outputs replace both results files; no assertions weakened.

### Immutable production artifact SHA-256 evidence
- `production_artifacts/feature_scaler.joblib`: `5a63f905e2744238cc376ac2aaa776de8d5bcac8574146f2322354d4267c623f` (current matches initial).
- `production_artifacts/frozen_model.joblib`: `0a3698e97e083d723566127c44abddf671f6cd8722e9b2509ed926e76ddd1989` (current matches initial).
- `production_artifacts/imputation_defaults.json`: `ee2d67c92ebe224b1829d230234d16ff6decf117c76b11d85d16ffa083e287c6` (current matches initial).
- `production_artifacts/manifest.json`: `fa411c200093dca0306faac4e54df60c049e7f07d7dca13f508aaf60bfa4c1e6` (current matches initial).

## Final POSIX launcher portability recheck
Raw-byte regression requires .command LF/no CR (does not impose LF on .bat). It caught the prior CRLF file: Python 3.11 initial 24-test run failed exactly this assertion. Authorized launcher worker then corrected only .command line endings. Assertion retained unchanged; both final 24-test runs pass with dependency skips listed above. Final .command raw CR count is zero.
Existing /usr/bin/bash executed `bash -n Launch_Diagnostic_Portal.command` against the actual corrected file: exit 0, no output. This is static parsing only, not launcher/server execution. Earlier normalized-text syntax checking was not raw-file portability evidence. Final git diff --check clean, staged empty, six changed/350 unchanged snapshot and prior 21-file audit integrity rechecked after LF fix.
