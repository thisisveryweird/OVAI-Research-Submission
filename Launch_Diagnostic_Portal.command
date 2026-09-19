#!/bin/bash
set -e
cd "$(dirname "$0")"
python3 -m pip install -r requirements.txt >/dev/null 2>&1
exec python3 -m streamlit run scripts/ovarian_diagnostic_portal.py --browser.gatherUsageStats=false