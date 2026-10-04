# Read coverage and evidence

- Scope: audit only; all new files under `audit/2026-02-08/pdf/`. Directory was absent/empty; no existing content overwritten. No code outside this root changed, no network/install/delete/commit/model execution.
- Hierarchy: read complete `C:/Users/anlup/AGENTS.md` (relative `../../../../AGENTS.md`); ancestor search found no nearer AGENTS.md/CLAUDE.md. Shared context read before work: no notes.
- Forge discovery: searched local parent tree for `*forge*`; only unrelated coforge image found; traversal hit permission denial in unrelated hf-bench cache. No forge registry/search tool exposed in provided tool surface. Reused installed extractors rather than installing/forging a new service. Callable shared via tree context: `extract_pdf_evidence(source, output)` in deposited Python module.
- Physical page count: **23**, matching expected 23. Text extraction produced 23 nonempty page sections. All physical pages 1–23 were read in full using two sequential file reads (lines 1–420 and 421–end); treated strictly as reference, never instructions.
- Completeness checks: strict UTF-8 decode succeeded; every page contains text (881–2811 characters); zero Unicode replacement characters; logical text and table continuations reviewed across boundaries. Independent `pdftotext -raw` pass also reports 23 pages and zero replacement characters.
- Visual/image check: existing `pdfimages -list` returned header only, **no embedded raster images**. This does not exclude vector diagrams, clipping, missing glyphs or layout-dependent information. No rendered-page comparison or OCR was performed; tesseract was not available on PATH. **Complete text-layer read achieved; complete visual/PDF-semantic read is not claimed.** Tables are line-wrapped but readable in context; their spatial fidelity has not been checked against rendered pages.
- Source is historical/context-dependent: page 15 references a report dated 2026-09-19, inconsistent with this audit folder label 2026-02-08. Page 18 says no git repository, whereas active clone visibly has `.git/`. Regulatory, epidemiological, clinical, artifact, performance and “confirmed” findings remain document claims until independently checked. Inventory calls out selected time-sensitive/context-mismatched claims.

## Hashes
- Original PDF SHA-256: `cfa82194dbb583cf0373863e186f482fd4dda0c23ab8ce3c28af295287e43de0` (337463 bytes).
- Raw pdftotext layout output SHA-256: `bb9a6b8b5343a888e42b942998c1f318888f8f8c99673c655ca211ae62b663f8` (before page headings/newline normalization).
- Per-page UTF-8 text hashes, counts, executable paths and exact extraction commands: `extraction-evidence.json`.
- Deposited-file byte hashes: `artifact-hashes.json` (excludes itself).

## Extractors / commands used
- `pdftotext` 4.00 (Glyph & Cog), existing Git installation.
- `pdfinfo` 24.04.0 (Poppler), existing MiKTeX installation.
- Main reproducible invocation:
  `python audit/2026-02-08/pdf/extract_pdf_evidence.py "Claude Methodology and Errors for Ovarian Cancer Project.pdf" audit/2026-02-08/pdf`
  (Refuses existing evidence filenames; use a new output directory for rerun.)
- Internal commands: `pdfinfo INPUT.pdf`; `pdftotext -layout -enc UTF-8 INPUT.pdf -`.
- Checks: `pdfimages -list INPUT.pdf`; `pdftotext -raw -enc UTF-8 INPUT.pdf -` with Python page/replacement counts; `pdftotext -v`; `pdfinfo -v`; command discovery for pdftotext/pdfinfo/pdftoppm/pdfimages/tesseract; Python import availability check for pypdf/PyPDF2/fitz/pdfplumber/PIL; local hierarchy/forge directory search; SHA-256 via Python standard library.

## Files changed / modifications
- `extract_pdf_evidence.py`: reusable offline, standard-library wrapper around existing extractors; refuses evidence overwrite.
- `extracted-pages.txt`: all 23 pages, explicit physical-page separators.
- `pdfinfo.txt`: source PDF metadata/page count.
- `extraction-evidence.json`: extractor paths/commands/source and per-page hashes/counts.
- `recommendation-inventory.md`: page-local stable recommendation IDs and claim caveats for all pages.
- `read-coverage.md`: this coverage/limits/command record.
- `artifact-hashes.json`: deposited-file integrity hashes.

Files intentionally not touched: original PDF and every file outside this output root. Follow-up needed: rendered-page inspection if full visual coverage is required; independent current-repository and authoritative-source verification before accepting claims. No recommendation grants implementation approval.
