@echo off
setlocal
cd /d "%~dp0"
if errorlevel 1 exit /b 1
where python >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python is unavailable. Prepare an approved Python environment separately and add python to PATH. 1>&2
    exit /b 1
)
echo Checking research portal environment readiness...
python -B -c "import streamlit, joblib, numpy, pandas, openpyxl"
if errorlevel 1 (
    echo ERROR: Research portal environment is not ready. Prepare an approved environment separately with working streamlit, joblib, numpy, pandas, and openpyxl imports. No dependencies were installed. 1>&2
    exit /b 1
)
python -m streamlit run scripts\ovarian_diagnostic_portal.py --browser.gatherUsageStats=false