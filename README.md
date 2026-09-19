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
8. GDC-native TCGA-OV file inventory and case metadata acquisition.

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

## Add TCGA-OV from GDC

Create a reproducible inventory of the TCGA-OV cases and available files:

```bash
python scripts/download_gdc_tcga_ov.py
```

This writes `data/external/gdc/tcga_ov/files.tsv`, `cases.json`, and
`query.json`. The inventory classifies GDC files as clinical, mutations,
mRNA, miRNA, CNV, methylation, or other. It preserves GDC file IDs, checksums,
case IDs, submitter IDs, access level, and workflow metadata so downstream
parsers can select one defined release and assay representation.

To download the open-access files for selected modalities:

```bash
python scripts/download_gdc_tcga_ov.py --download \
	--modality mrna --modality mirna --modality cnv
```

Downloaded payloads are stored below `data/external/gdc/tcga_ov/raw/`. The
default command only retrieves metadata; it does not download the large cohort.
Controlled-access files are recorded in the inventory but are not downloaded.
Raw GDC files must be normalized to patient-level `patient_id` tables before
being joined into a multimodal model, as required by the schema.

## Evaluate TCGA-OV stage model

Run the leakage-safe stage evaluation with:

```bash
python scripts/train_tcga_ov_stage_model.py
```

The script evaluates Stage I/II versus Stage III/IV using nested five-fold
out-of-fold predictions and saves the feature table and metrics under
`data/processed/` and `reports/tables/`. It reports raw accuracy, balanced
accuracy, ROC-AUC, sensitivity, specificity, and a majority-class baseline.
This endpoint is tumor-stage research only; it is not a cancer-detection or
clinical diagnostic model.

## Evaluate GSE94533 diagnosis

GSE94533 is kept as an independent diagnostic cohort rather than merged with
TCGA patients. Build it and evaluate the historical signature with:

```bash
python scripts/download_public_data.py --dataset gse94533_cohort1
python scripts/download_public_data.py --dataset gse94533_cohort2
python scripts/build_gse94533_dataset.py
python scripts/build_gse94533_clinical_fusion.py
python scripts/evaluate_gse94533_diagnostic_model.py
```

The evaluation uses current GEO sample annotations and records the resulting
96 cancer / 83 non-cancer composition. The older hard-coded publication count
check was removed because it does not match the current GEO annotations; the
excluded E-1056 sample and all source labels remain auditable in the manifest.
The eLife supplement provides matched Age for 178 patients and CA-125 for 119;
missing values are imputed inside each training fold. In the current run,
adding Age and CA-125 to the 14-miRNA signature improved out-of-fold accuracy
from 70.39% to 73.74% and ROC-AUC from 0.779 to 0.813.

The requested `phs000897` dbGaP study is controlled-access. It requires the
user's dbGaP authorization and cannot be downloaded or added automatically
from this environment. It should be added later as an independent early-stage
validation/training cohort only after its consent and patient-level modality
mapping are verified.

Proteomic repositories such as CPTAC, iProX, and GSA are not merged into this
diagnostic table: their tissue/serum assays, patient identifiers, and labels do
not establish a verified one-to-one overlap with GSE94533. Likewise, TCGA-OV
is not used to pretrain diagnostic weights, because its available endpoint is
tumor stage and its tissue expression distribution is not the serum cancer vs
non-cancer endpoint. Such transfer would risk improving apparent performance
through domain leakage rather than improving diagnosis.

## Enhancement benchmark

Run the algorithm comparison with:

```bash
python scripts/evaluate_gse94533_enhancements.py
```

This evaluates median versus KNN imputation, calibrated SVM, regularized
random forest, shallow XGBoost, and an XGBoost/RF/SVM stack. Age/CA-125
interaction terms and fold-contained LASSO stability are included. The
current results do not support replacing the probability-ranking model with
the stack: the stack ROC-AUC is 0.800 versus 0.813 for the fused logistic
baseline. A calibrated SVM with an inner-fold-selected threshold is the best
classification option observed so far, reaching 76.54% accuracy and 77.06%
balanced accuracy in outer-fold predictions. These are research metrics, not
clinical validation.

No `phs000897` rows are present, so independent dbGaP testing and ComBat
correction cannot be run. ComBat should only be added when an external serum
cohort with batch metadata and compatible measurements is actually available.

The additional validation experiments are recorded by
`scripts/evaluate_gse94533_enhancements.py` and
`scripts/evaluate_gse106817_external.py`. The four stable-miRNA log-ratios did
not improve balanced accuracy, and a 95%-specificity threshold reduced
sensitivity too far for a useful overall classifier. GSE106817 was tested as a
frozen external serum cohort; its ranking ROC-AUC was 0.823 with standard
scaling, but its training-derived threshold produced 0% sensitivity. This is
platform shift evidence, not a validated diagnostic result.

## Methylation benchmark

The public MethylBERT-EOC supplement can be evaluated independently with:

```bash
python scripts/evaluate_methylbert_methylation.py
```

This downloads the Europe PMC supplementary archive and runs nested five-fold
Ridge classification with fold-contained imputation and top-feature selection.
The processed methylation matrix reached 88.14% accuracy, 86.90% balanced
accuracy, and 0.941 ROC-AUC for ovarian cancer versus healthy controls. The
early-stage-only versus healthy comparison reached 84.13% accuracy, 77.26%
balanced accuracy, and 0.861 ROC-AUC. These are independent cohort results;
the methylation rows were not merged with GSE94533 patients.

GSE94533 provides only processed miRNA Excel files in its public supplement,
not raw FASTQ/BAM files, so an isomiR sequence transformer cannot be trained
from that accession. The MethylBERT supplement is also processed data, not raw
BAM/FASTQ; the benchmark therefore does not claim raw-fragment modeling.

## Frozen methylation validation

The MethylBERT Ridge pipeline can be frozen with:

```bash
python scripts/freeze_methylbert_pipeline.py
python scripts/audit_methylbert_external_targets.py
```

The frozen artifact is written to `models/` and is never refit on a validation
cohort. The proposed GSE144365 target is invalid: it is a *Pseudomonas
aeruginosa* transcriptome study. TCGA-OV HumanMethylation450 is also not a
valid exact-feature test: it has `cg...` Illumina probe IDs, while the
MethylBERT matrix has genomic-coordinate features, producing zero exact
overlap, and it is tissue rather than blood. No external accuracy is reported
for either target. HRA011451 remains unavailable as a separate, feature-aligned
validation matrix.
