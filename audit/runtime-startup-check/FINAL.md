# Runtime startup results

**Portal blocked by missing dependencies; no running UI verified.** Actual runtime imports and CLI help were executed, not AST-only tests.

| Check | Result | Evidence |
|---|---|---|
| Default `python` | Python 3.14.2 | default-version.json |
| Default readiness | Exit 1: `ModuleNotFoundError` for streamlit, joblib, pandas, openpyxl, sklearn; numpy 2.4.6 imports | default-readiness.json |
| Existing `py -3.11` | Python 3.11.9 | py311-version.json |
| Python 3.11 readiness | Exit 1: `ModuleNotFoundError: No module named 'streamlit'` and `No module named 'openpyxl'`; joblib 1.5.3, numpy 1.26.4, pandas 2.3.3, sklearn 1.7.2 import successfully | py311-readiness.json |
| Prediction CLI `--help` (3.11) | Exit 0; empty stderr | predict-help.json |
| Baseline CLI `--help` (3.11) | Exit 0; empty stderr | baseline-help.json |

Each JSON records command, stdout, stderr, and exit code. `readiness_check.py` is a reusable stdlib probe, compiled in memory successfully and run with bytecode disabled. Environment flags disabled bytecode and Streamlit usage telemetry.

## Startup limitation
The Windows launcher uses default `python` and its readiness imports would fail in this environment. The macOS launcher uses `python3`; it was inspected but not executed on Windows. Neither existing tested interpreter is portal-ready. Following the dependency-failure stop condition, no Streamlit process, AppTest, HTTP health check, or app initialization was attempted. There is no PID or localhost URL. Successful CLI help establishes imports/argument parsing only, not model inference or training viability.

No models were loaded; artifact before/after hashing is therefore not applicable. No datasets, clinical predictions, training run, demo patient-file reads, installs, environment creation, product edits, or network requests were performed. Existing repository modifications were not changed. An initial shell command using Unix separators failed on the Windows command shell before executing version probes; it was corrected to Windows syntax.

## Suggested next steps — not executed
After separate approval to prepare dependencies in the existing Python 3.11 environment:

```bat
py -3.11 -m pip install streamlit openpyxl
py -3.11 -B audit\runtime-startup-check\readiness_check.py
```

Only after readiness passes, rerun startup verification with bytecode disabled and telemetry off. A candidate local-only startup command (not executed):

```bat
set PYTHONDONTWRITEBYTECODE=1
set STREAMLIT_BROWSER_GATHER_USAGE_STATS=false
py -3.11 -B -m streamlit run scripts\ovarian_diagnostic_portal.py --server.address=127.0.0.1 --server.port=8501 --server.headless=true --server.fileWatcherType=none --browser.gatherUsageStats=false
```

A future check must distinguish server health from actual app initialization and must not trigger prediction controls. These proposed commands are not evidence that installation or app initialization will succeed.
