# Ovarian Multimodal AI

Research prototype for studying whether combining symptomatic/clinical information with multiple genomic modalities can improve detection of epithelial ovarian cancer, with a particular focus on early-stage disease.

## Current implementation

1. GSE94533 serum miRNA-seq ingestion and dataset audit.
2. Patient-level train/test split and miRNA baseline model.
3. Multimodal patient table schema with modality presence masks.
4. Gated multimodal neural network scaffold in PyTorch.
5. Evaluation including AUROC, AUPRC, sensitivity/specificity, calibration and confusion matrix.
6. MethylBERT-related methylation data scaffold.
7. TCGA-OV reference-data scaffold.

## The model should ultimately compare

- Clinical/symptoms only
- Clinical + conventional biomarkers/imaging
- Each genomic modality separately
- Clinical + all available genomic modalities

## Important scientific rule

Never concatenate data from unrelated patients across different studies and pretend that the resulting row represents one patient. Multimodal fusion requires patient-linked modalities. Separate cohorts are useful for representation learning, discovery, reproduction, or external validation, but not for inventing patient-level multimodal rows.

## Run

```bash
cd ~/ovarian-ai
source .venv/bin/activate
bash scripts/run_public_mirna_baseline.sh
```

For the future multimodal model, place patient-aligned modality tables in `data/processed/` using the schema documented in `configs/multimodal_schema.yaml`.
