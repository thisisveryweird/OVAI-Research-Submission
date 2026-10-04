#!/bin/bash
set -e
cd "$(dirname "$0")"
if ! command -v python3 >/dev/null 2>&1; then
    echo "ERROR: Python is unavailable. Prepare an approved Python environment separately and add python3 to PATH." >&2
    exit 1
fi
echo "Checking research portal environment readiness..."
if ! python3 -B -c "import streamlit, joblib, numpy, pandas, openpyxl"; then
    echo "ERROR: Research portal environment is not ready. Prepare an approved environment separately with working streamlit, joblib, numpy, pandas, and openpyxl imports. No dependencies were installed." >&2
    exit 1
fi
exec python3 -m streamlit run scripts/ovarian_diagnostic_portal.py --browser.gatherUsageStats=false