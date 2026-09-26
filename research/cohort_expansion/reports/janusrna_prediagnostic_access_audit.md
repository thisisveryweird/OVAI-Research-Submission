# JanusRNA Prediagnostic Cohort Access Audit

Source: Stawiski et al., "Validation of miRNA signatures for ovarian cancer earlier detection in the pre-diagnosis setting using machine learning approaches" (2024; DOI `10.3389/fonc.2024.1389066`; PMCID `PMC11231195`).

The study used small-RNA sequencing from prediagnostic **serum** in the Norwegian Janus Serum Bank. The source cohort selected 80 invasive ovarian cancer cases and had 106 serum samples before one sample per person was chosen for the published analysis: 14 people had two samples and 4 people had three or more. The analysis matched 40 cases diagnosed <2 years after collection against 40 cases diagnosed >=2 years after collection; no cancer-free controls remained in this analysis because a batch effect between case and control libraries led the authors to discard controls.

The paper states that data will be made available through Federated EGA Norway after application to the study Data Access Committee. Individual-level expression is therefore **access controlled**, not publicly downloadable. This is a viable access-gated prediagnostic biology study, not an open cohort today.

Use requires patient-grouped analysis of all repeated samples, collection-to-diagnosis intervals, and the published choice of one sample per person. The paper reports a 2-year proximity contrast, serum sequencing, and significant platform differences that prevented direct transfer of the original model. No canonical OVAI prediction or threshold tuning on JanusRNA is justified from the public evidence.

**Decision:** YELLOW for a prediagnostic biological-timing study after EGA approval and access to sample-level clinical linkage; RED for direct current-OVAI diagnostic validation because it is case-only in the released analysis, lacks benign pelvic-mass controls, is controlled-access, and differs technically from GSE94533.
