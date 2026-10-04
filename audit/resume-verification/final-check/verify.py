"""Bounded, offline structural verification; stdlib only. No source writes."""
import csv
import json
import re
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent

def load(relative):
    with (ROOT / relative).open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        return reader.fieldnames, rows

def duplicates(rows, key):
    groups = {}
    for row in rows:
        groups.setdefault(row[key], []).append(row)
    return {k: {'count': len(v), 'conflict': any(x != v[0] for x in v[1:])}
            for k, v in groups.items() if len(v) > 1}

def target_check(references, targets):
    return sorted(set(references) - set(targets))

assert duplicates([{'id':'a','v':'1'},{'id':'a','v':'2'}], 'id') == {'a': {'count':2,'conflict':True}}
assert duplicates([{'id':'a'},{'id':'a'}], 'id')['a']['conflict'] is False
assert target_check(['a','unknown'], ['a']) == ['unknown']

base = 'audit/2026-02-08/review/'
tables = {}
for name, key in [('PDF_CHECKLIST.csv','recommendation_id'),('ENDPOINT_REGISTRY.csv','family_id'),('ENDPOINT_ALIASES.csv','alias_id')]:
    headers, rows = load(base + name)
    tables[name] = (headers, rows)
summary = {'authorization_for_repairs':False, 'scope':'offline structural checks plus attributed semantic spot-checks', 'fixture_tests':'passed', 'execution_note':'Corrected copy run; prior helper failures preserved; malformed citation fields isolated.', 'tables':{}}
for name, (headers, rows) in tables.items():
    key = {'PDF_CHECKLIST.csv':'recommendation_id','ENDPOINT_REGISTRY.csv':'family_id','ENDPOINT_ALIASES.csv':'alias_id'}[name]
    summary['tables'][name] = {'headers':headers, 'rows':len(rows), 'unique_ids':len({r[key] for r in rows}), 'duplicates':duplicates(rows,key)}
checklist = tables['PDF_CHECKLIST.csv'][1]
registry = tables['ENDPOINT_REGISTRY.csv'][1]
aliases = tables['ENDPOINT_ALIASES.csv'][1]
requested = ROOT / 'audit/2026-02-08/pdf/recommendation-inventory.md'
inventory = requested if requested.is_file() else None
candidates = [requested] if inventory else []
ids = set(re.findall(r'PDF-P\d+-R\d+', inventory.read_text(encoding='utf-8-sig'))) if inventory else set()
check_ids = {r['recommendation_id'] for r in checklist}
summary['inventory'] = {'requested_path_exists':requested.exists(), 'resolved_path':str(inventory.relative_to(ROOT)).replace('\\','/') if inventory else None, 'candidate_count':len(candidates), 'unique_ids':len(ids), 'checklist_id_equality':ids == check_ids if inventory else None, 'only_checklist':sorted(check_ids-ids), 'only_inventory':sorted(ids-check_ids)}
refs = [x for r in registry for x in re.findall(r'EP-\d+',r['evidence_aliases'])]
summary['alias_relationship'] = {'direction':'registry.evidence_aliases -> aliases.alias_id (actual schema has no target-family column)', 'reference_count':len(refs), 'unique_references':len(set(refs)), 'unknown_targets':target_check(refs,[r['alias_id'] for r in aliases]), 'unreferenced_alias_ids':target_check([r['alias_id'] for r in aliases],refs)}
issues = []
counts = Counter()
def citation(table, row_id, path, lines=None, physical_row=False):
    try:
        return checked_citation(table, row_id, path, lines, physical_row)
    except (OSError, ValueError, TypeError) as exc:
        issues.append({'table':table,'id':row_id,'issue':'malformed citation','detail':str(exc)[:120]})

def checked_citation(table, row_id, path, lines=None, physical_row=False):
    path = str(path).strip().strip('`')
    if not path:
        return
    counts['path_citations'] += 1
    p = ROOT / path
    if not p.is_file():
        issues.append({'table':table,'id':row_id,'path':path,'issue':'path missing or not a single file'})
        counts['missing_paths'] += 1
        return
    counts['existing_paths'] += 1
    if lines:
        nums = re.fullmatch(r'\s*(\d+)\s*(?:[-â€“]\s*(\d+))?\s*', str(lines))
        if not nums:
            counts['non_numeric_line_specs'] += 1
            issues.append({'table':table,'id':row_id,'path':path,'issue':'line specification not numeric','lines':lines})
            return
        a, b = int(nums[1]), int(nums[2] or nums[1])
        n = len(p.read_text(encoding='utf-8-sig',errors='replace').splitlines())
        counts['physical_row_checks' if physical_row else 'numeric_line_checks'] += 1
        if not 1 <= a <= b <= n:
            counts['out_of_bounds'] += 1
            issues.append({'table':table,'id':row_id,'path':path,'issue':'physical CSV row bounds invalid' if physical_row else 'line bounds invalid','lines':lines,'file_lines':n})
for r in checklist:
    citation('PDF_CHECKLIST',r['recommendation_id'],r['source_path'],r['inclusive_lines'])
for r in aliases:
    citation('ENDPOINT_ALIASES',r['alias_id'],r['lineage_file'],r['lineage_row'],physical_row=True)
def decode_fields(value, family):
    text = (value or '').strip()
    if not text or text.lower() in {'unknown','unresolved','not available','unavailable','n/a','na','none','null'} or text.lower().startswith(('unknown:', 'unresolved:')):
        counts['source_fields_unavailable'] += 1
        issues.append({'table':'ENDPOINT_REGISTRY','id':family,'issue':'source fields unavailable','raw':text[:120]})
        return {}
    try:
        obj = json.loads(text)
        if not isinstance(obj, dict):
            raise ValueError('JSON value is not an object')
        return obj
    except (ValueError, TypeError) as exc:
        counts['source_fields_parse_issues'] += 1
        issues.append({'table':'ENDPOINT_REGISTRY','id':family,'issue':'source_fields_json parse issue','raw':text[:120],'detail':str(exc)[:120]})
        return {}
assert decode_fields('{}','fixture') == {}
assert decode_fields('unknown','fixture') == {}
assert decode_fields('literal placeholder','fixture') == {}
assert decode_fields('[]','fixture') == {}
for key in ['source_fields_unavailable','source_fields_parse_issues']: counts.pop(key,None)
issues[:] = [x for x in issues if x['id'] != 'fixture']
for r in registry:
    raw = decode_fields(r['source_fields_json'], r['family_id'])
    try:
        if raw.get('evidence_file'):
            citation('ENDPOINT_REGISTRY',r['family_id'],raw['evidence_file'],raw.get('evidence_row'),physical_row=True)
        if isinstance(raw.get('source'),str) and (ROOT / raw['source']).is_file():
            citation('ENDPOINT_REGISTRY',r['family_id'],raw['source'],raw.get('lines'))
    except (OSError, ValueError, TypeError) as exc:
        issues.append({'table':'ENDPOINT_REGISTRY','id':r['family_id'],'issue':'malformed citation','detail':str(exc)[:120]})
summary['citation_scope'] = 'checklist source_path/inclusive_lines; alias lineage_file/lineage_row (physical CSV row existence only, not source-code relevance); registry source_fields_json evidence_file/evidence_row and existing literal source/lines. Free-text evidence/lineage and nested arbitrary citations not exhaustively parsed.'
summary['citations'] = dict(counts)
summary['citation_issues'] = issues
implemented = [r['recommendation_id'] for r in checklist if r['status']=='implemented_verified']
static = {
 'PDF-P05-R01':'Root-read static evidence: implemented_verified too broad for multimodal_table.py: duplicates rejected and ID reindex done; no donor/visit provenance; modality label not excluded; values converted to zero rather than NaN despite docstring. Suggest partial/source-only.',
 'PDF-P07-R02':'Root-read static evidence: methylation is a standalone processed-matrix endpoint with global eligibility filters, not linked miRNA improvement.',
 'PDF-P07-R03':'Root-read static evidence: TCGA source docstring explicitly tumor-stage, not detection; code labels III/IV versus I/II, not Stage I alone.',
 'PDF-P18-R02':'Root-read static evidence: TCGA source docstring explicitly tumor-stage, not detection; code labels III/IV versus I/II, not Stage I alone.',
 'PDF-P14-R05':'Root-read static evidence: companion inventory contains page15 introduction; inventory consulted; PDF visuals remain unverified.'}
summary['implemented_verified'] = {'count':len(implemented),'ids':implemented,'not_substantively_checked':[i for i in implemented if i not in static]}
summary['root_read_static_evidence'] = static
summary['integrity_attributed_not_rerun'] = {'report':'audit/resume-verification/integrity/REPORT.md','review_hashes':'11/11','historical_hashes':'212/212','preexisting_hashes':'356/356 unchanged','tracked_changes':'no staged/unstaged tracked changes'}
placeholder = ROOT / 'audit/resume-verification/semantics/REPORT.md'
summary['prior_placeholder'] = {'report_exists':placeholder.exists(),'preserved':True,'referenced_findings_json_exists':(placeholder.parent/'findings.json').exists()}
spotchecks = {
 'PDF-P08-R06': 'Direct static read scripts/evaluate_gse94533_diagnostic_model.py:62-96: full-frame notna/nunique eligibility precedes outer CV; unresolved remains supported, not repaired.',
 'PDF-P08-R05': 'ATTRIBUTED parent-read discovery-final/REPORT.md section 2: Phase8 inner contamination supported; no independent recomputation.',
 'PDF-P06-R03': 'ATTRIBUTED parent-read external-final/REPORT.md Critical reconciliations: bulk all-positive differs from selected historical14 all-negative arms.',
 'PDF-P18-R01': 'ATTRIBUTED parent-read historical-reconstruction/REPORT.md Counts/splits: original 135/44 split not recovered.',
 'inference feature identity': 'ATTRIBUTED parent-read inference-final/REPORT.md finding 2: seven versus nine feature discrepancy supported.'}
summary['semantic_spotchecks'] = spotchecks
sample_ids = ['F-010','F-023','G-04-06','V-001','V-006','G-01-01','F-034','V-018','V-007','V-008','V-009','V-010']
alias_map = {r['alias_id']:r for r in aliases}
summary['identity_samples'] = [{'family_id':r['family_id'],'pipeline':r['pipeline'],'endpoint':r['endpoint'],'limits':r['limits'],'aliases':[{'alias_id':a,'exists':a in alias_map,'lineage_file':alias_map.get(a,{}).get('lineage_file'),'physical_row':alias_map.get(a,{}).get('lineage_row')} for a in re.findall(r'EP-\d+',r['evidence_aliases'])]} for r in registry if r['family_id'] in sample_ids]
summary['identity_sample_limit'] = 'Purposeful identity/alias consistency samples across five topics, not clinical counts, fitted-artifact compatibility, or whole-registry semantic validation. Historical14 arm is one sampled arm, not all arms.'
summary['pending'] = ['Whole-registry substantive semantics beyond reviewed spot-checks','Fresh raw donor/data validation','PDF visual verification','Other implemented_verified row substantive checking','Free-text citation interpretation']
summary['completed'] = ['CSV row/unique-ID and duplicate-conflict counts','Inventory/checklist ID comparison','Registry alias target membership','Scoped cited-file existence and numeric line bounds','Tiny in-memory duplicate/unknown-target fixture tests','Materialization of attributed root-read and integrity evidence; no rerun']
for name in ['SUMMARY.json','REPORT.md','FILES_CREATED.csv']:
    if (OUT/name).exists():
        raise RuntimeError('Refusing to overwrite '+name)
with (OUT/'SUMMARY.json').open('x',encoding='utf-8') as f: f.write(json.dumps(summary,indent=2)+'\n')
lines = ['# Bounded offline verification','', '**authorization_for_repairs=false**. No existing files modified. Structural checks do not establish scientific correctness.', '', '## Computed results']
for name, data in summary['tables'].items():
    lines.append(f"- {name}: {data['rows']} rows, {data['unique_ids']} unique IDs; {len(data['duplicates'])} duplicate-ID groups.")
lines += [f"- Inventory: {summary['inventory']}",f"- Alias relationship: {summary['alias_relationship']}",f"- Citation checks: {dict(counts)}; {len(issues)} issues (details in SUMMARY.json).", '- Citation scope: '+summary['citation_scope'], f'- implemented_verified actual count: {len(implemented)}; IDs: '+', '.join(implemented), '', '## Root-read static evidence (attributed; not independently rerun)']
lines += ['- '+v for v in static.values()]
lines += ['', '## Semantic spot-checks (attributed except direct diagnostic source)'] + ['- '+k+': '+v for k,v in spotchecks.items()]
lines += ['', '## Purposeful identity samples', '- '+summary['identity_sample_limit']]
for sample in summary['identity_samples']:
    lines.append('- '+sample['family_id']+': '+sample['pipeline']+'; aliases '+', '.join(a['alias_id']+' (exists='+str(a['exists'])+', physical row='+str(a['physical_row'])+')' for a in sample['aliases'])+'; limits: '+sample['limits'])
lines += ['- Other implemented_verified rows not substantively checked: '+(', '.join(summary['implemented_verified']['not_substantively_checked']) or 'none'), '', '## Completed'] + ['- '+s for s in summary['completed']]
lines += ['', '## Pending â€” not newly verified'] + ['- '+s for s in summary['pending']]
lines += ['', '## Prior evidence preserved', '- Existing semantics/REPORT.md remains untouched; its referenced findings.json is absent (existence result: '+str(summary['prior_placeholder']['referenced_findings_json_exists'])+').', '- Integrity evidence attributed to audit/resume-verification/integrity/REPORT.md, not rerun: 11/11 review hashes; 212/212 historical hashes; 356/356 preexisting hashes unchanged; no staged/unstaged tracked changes.', '- No network, models, project execution, installs, generators, or git changes. Helper fixture assertions passed.']
with (OUT/'REPORT.md').open('x',encoding='utf-8') as f: f.write('\n'.join(lines)+'\n')
with (OUT/'FILES_CREATED.csv').open('x',encoding='utf-8',newline='') as f:
    w = csv.writer(f)
    w.writerow(['path','purpose'])
    for name, purpose in [('verify.py','Standalone stdlib bounded verifier with tiny fixtures'),('SUMMARY.json','Compact computed counts issues and attributed evidence'),('REPORT.md','Completed and pending scope with repair authorization false'),('FILES_CREATED.csv','Exact inventory of newly created artifacts')]:
        w.writerow([(OUT/name).relative_to(ROOT).as_posix(),purpose])
print(json.dumps({'rows':{n:d['rows'] for n,d in summary['tables'].items()},'inventory':summary['inventory']['resolved_path'],'id_equality':summary['inventory']['checklist_id_equality'],'implemented_verified':len(implemented),'citation_issues':len(issues),'artifacts':4}))
