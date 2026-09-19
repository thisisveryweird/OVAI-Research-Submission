@echo off
setlocal
cd /d "%~dp0"
python -m pip install -r requirements.txt >nul 2>&1
if errorlevel 1 exit /b 1
python -m streamlit run scripts\ovarian_diagnostic_portal.py --browser.gatherUsageStats=false