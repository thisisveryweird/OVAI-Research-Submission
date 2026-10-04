# Manifest / preservation verification

Repair authorization remains **false**. No project imports, generators, model deserialization, network, installs, training, git mutation, deletion or existing-file edits.

- Current review manifest: **11/11 SHA-256 matches**, zero missing or mismatched entries. Cached `review_audit.cpython-314.pyc` exists and matches. Review contains 12 files: exactly the 11 entries plus `HASH_MANIFEST.json`. Its intentional self-exclusion avoids a circular content hash; it does not cover later session files or the whole repository.
- `FILES_CHANGED.md`: 141 unique listed paths, all present. Review coverage is complete including itself and manifest. Across the current dated audit tree, three paths are unlisted: `integration/FILES_CHANGED.md`, `integration/HASH_MANIFEST.json`, `session/RESUME.md`. Thus it is not a complete current audit-tree inventory.
- Session directory contains only `RESUME.md`. Independent dated-audit recursive search found neither `FILES_CREATED.csv` nor `handoff_inventory.py`; these handoff references are absent artifacts, **not failures of the review hash manifest**.
- Historical preservation mapping independently read from `integration/PRESERVATION_CHECKS.json`: **212/212 current byte hashes match recorded initial hashes**, including owner PDF. This corroborates those bounded recorded paths, not historical capture authenticity or whole-tree coverage.
- HEAD: `5fd00c949aa89ed44446f5084a9fe03800ae4def`. All five read-only git commands exit 0. Tracked unstaged stat/name diff and staged name diff empty. Status has 145 lines, recorded verbatim in evidence; untracked output is not an original-product modification.
- Initial snapshot: **356 preexisting regular files**, excluding `.git` and all `audit/resume-verification/` sibling/output trees. Final recheck: **356/356 match**. Snapshot was taken after reading handoff/manifest and creating helper, before substantive checks; it demonstrates byte preservation within its stated boundary, not metadata preservation or absence of concurrent additions.

## Reusable helper and tests
`python -B audit/resume-verification/integrity/integrity_helper.py ROOT MANIFEST` reads a JSON path-to-SHA256 mapping, streams bytes, emits JSON to stdout and returns 0 only for all matches. No output-file mutation. Four synthetic assertions passed: match, mismatch, absent file, escaping root. Test fixture remains under this output root. `check_current.py` is a new offline driver (exclusive output creation); do not rerun in this populated directory.

## Commands and evidence
Working directory: repository root. Exact substantive commands:
1. `python -B -c "import runpy,json,pathlib; r=pathlib.Path('.'); o=r/'audit/resume-verification/integrity'; h=runpy.run_path(str(o/'integrity_helper.py')); files=[p for p in r.rglob('*') if p.is_file() and '.git' not in p.parts and not p.as_posix().startswith('audit/resume-verification/')]; d={p.as_posix():h['digest'](p) for p in files}; (o/'initial_snapshot.json').open('x').write(json.dumps(d,indent=2)); print('snapshot',len(d))"` — exit 0, 356.
2. `python -B audit/resume-verification/integrity/check_current.py` — exit 0; exact git argv, stdout/stderr and return codes in `evidence.json`.
3. Summary inspection via stdlib `json`/`collections` — exit 0; counts above. Preliminary stdlib inspection read manifest/list, session/root directory listings and bounded first 4500 characters of historical mapping. An attempted multiline inline helper creation failed parsing (SyntaxError, exit 1), before executing any filesystem change. New helper and driver then created through new-file write tool on previously absent paths; all driver outputs use exclusive `x` creation.

Machine evidence: `evidence.json`, `initial_snapshot.json`, `preservation_recheck.json`, `synthetic_results.json`. Exact created path inventory: `FILES_CREATED.csv`, including itself. New `HASH_MANIFEST.json` hashes every created file except itself; its self-exclusion is explicit and inventory still includes it. No semantic quality or whole-tree historical hash claim is made.
