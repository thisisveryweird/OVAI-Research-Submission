#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_DIR"
source .venv/bin/activate

mkdir -p reports/logs
LOG="reports/logs/public_mirna_baseline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "=== OVARIAN AI PUBLIC MI-RNA PIPELINE ==="
echo "Project: $PROJECT_DIR"
python --version

python scripts/download_public_data.py --dataset gse94533_cohort1
python scripts/download_public_data.py --dataset gse94533_cohort2
python scripts/inspect_public_data.py
python scripts/train_mirna_baseline.py --top-k 50
python scripts/smoke_test.py

echo
 echo "=== DONE ==="
echo "Audit: reports/tables/dataset_audit.json"
echo "Metrics: reports/tables/mirna_baseline_metrics.json"
echo "Selected features: reports/tables/mirna_selected_features.csv"
echo "Model: models/mirna_logreg_baseline.joblib"
echo "Log: $LOG"
