# Completed: local research portal

- Dependency installation succeeded using only `py -3.11 -m pip install streamlit openpyxl`. Streamlit 1.65.0/openpyxl 3.1.5 import successfully. Full newly installed dependency/version inventory is in 01-install.json; no broad upgrades/project requirements installation. Script-directory-not-on-PATH warnings occurred.
- Synthetic containment suite: 24/24 passed, including XLSX both orientations. This does not establish real inference validity.
- Actual Streamlit AppTest: passed initial rendering; no exceptions or st.error; research-only title present, schema/artifacts loaded, Execute Research Model disabled absent attestation. No buttons clicked, uploads supplied, or real predictions executed. Initial helper-relative path error was recorded in 04-apptest.json; retry used the requested relative app path from repository context and succeeded (05-apptest-retry.json).
- Compatibility warning: saved LogisticRegression and StandardScaler report sklearn 1.9.1, whereas runtime is 1.7.2. Loading succeeded but inference compatibility/results remain unverified. No attempt to fix or upgrade sklearn.
- Detached server left running: PID 3976, http://127.0.0.1:8501 . Health and root each returned HTTP 200; HTTP alone is not an initialization test. Interpreter/full detached command in launch.json. Existing processes were not stopped/overwritten.
- Reviewed complete portal source: initial main only loads the existing local production artifact set, validates schema/dimensions/classes/defaults and hashes, renders default inputs. Prediction requires attested execution or supplied bulk upload, neither triggered. Existing joblib artifacts are treated as trusted project inputs, not independently authenticated.
- All four production artifact SHA256 values unchanged before/after install/tests/AppTest/startup (artifact JSON records). No app/data/model edits, training, real patient data, or clinical inference performed. Research-only; not validated for clinical use.

The unchanged .bat launcher still defaults to Python 3.14. From the repository use explicit Python 3.11 instead (do not run a second server while this one is available):

```bat
set PYTHONDONTWRITEBYTECODE=1
set STREAMLIT_BROWSER_GATHER_USAGE_STATS=false
py -3.11 -B -m streamlit run scripts/ovarian_diagnostic_portal.py --server.address=127.0.0.1 --server.port=8501 --server.headless=true --server.fileWatcherType=none --browser.gatherUsageStats=false
```

Reusable stdlib orchestration/probe helpers: probe.py and retry.py; they compile in memory without bytecode. Exclusive-new evidence files only in this directory. Server stdout/stderr are live logs. FILES_CREATED.csv inventories this directory.
