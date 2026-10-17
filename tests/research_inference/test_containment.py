import argparse
import ast
import io
import json
import math
import re
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from ast_helpers import ROOT, source, extract
try:
    import numpy as np
except ImportError:
    np = None
try:
    import pandas as pd
except ImportError:
    pd = None

PORTAL = 'ovarian_diagnostic_portal.py'
CLI = 'predict_ovarian_risk.py'
MIRNA = ('hsa-miR-200b-3p', 'hsa-miR-296-3p', 'hsa-miR-6826-5p', 'hsa-miR-32-3p', 'hsa-miR-1255a')
SEVEN = MIRNA + ('age', 'ca125')
NINE = SEVEN + ('mirna_200b_ca125', 'mirna_296_age')

def helpers(name, funcs):
    return extract(name, funcs, dict(np=np, pd=pd, re=re, json=json, math=math, Path=Path,
                                    EXPECTED_FEATURES=SEVEN if name == CLI else NINE))

class SchemaTests(unittest.TestCase):
    @unittest.skipIf(np is None, 'numpy unavailable: model class validation runtime')
    def test_model_class_order_and_probability_range(self):
        for name, expected in ((CLI, SEVEN), (PORTAL, NINE)):
            ns = helpers(name, ['validate_artifact_schema', 'validate_probabilities'])
            for classes in (None, [1, 0], [0], [0, 1, 2], ['0', '1']):
                with self.subTest(name=name, classes=classes), self.assertRaises(ValueError):
                    ns['validate_artifact_schema']({'features': list(expected)}, model=SimpleNamespace(n_features_in_=len(expected), classes_=classes))
            for value in (float('nan'), float('inf'), -0.1, 1.1):
                with self.assertRaises(ValueError): ns['validate_probabilities'](value)
                with self.assertRaises(ValueError): ns['validate_probabilities']([.2, value])
            for value in (0., .5, 1.): ns['validate_probabilities'](value)

    @unittest.skipIf(np is None, 'numpy unavailable: CLI runtime')
    def test_cli_hyphenated_arguments_reach_dataframe(self):
        ns = helpers(CLI, ['main', 'parse_number', 'validate_artifact_schema'])
        class FakePath:
            def __init__(self, name=''): self.name = name
            def __truediv__(self, name): return FakePath(name)
            def read_text(self):
                return json.dumps({'features': list(SEVEN)} if self.name == 'manifest.json' else {'age': 50., 'ca125': 20.})
        loads = []
        def load(path):
            loads.append(path.name)
            return SimpleNamespace(n_features_in_=7, classes_=[0, 1])
        class Reached(Exception): pass
        def dataframe(rows, columns):
            self.assertEqual(rows, [dict(zip(SEVEN, range(1, 8)))])
            self.assertEqual(columns, list(SEVEN))
            raise Reached()
        ns.update(argparse=argparse, DEFAULT_ARTIFACT_DIR=FakePath(), joblib=SimpleNamespace(load=load), pd=SimpleNamespace(DataFrame=dataframe))
        args = ['synthetic']
        for value, feature in enumerate(SEVEN, 1): args += ['--' + feature, str(value)]
        with patch('sys.argv', args), self.assertRaises(Reached): ns['main']()
        self.assertEqual(loads, ['frozen_model.joblib', 'feature_scaler.joblib'])

    @unittest.skipIf(pd is None or np is None, 'pandas/numpy unavailable: CLI numeric runtime')
    def test_cli_synthetic_numeric_prediction_and_nonfinite_rejection(self):
        ns = helpers(CLI, ['main', 'parse_number', 'validate_artifact_schema', 'validate_probabilities'])
        class FakePath:
            def __init__(self, name=''): self.name = name
            def __truediv__(self, name): return FakePath(name)
            def read_text(self): return json.dumps({'features':list(SEVEN)} if self.name == 'manifest.json' else {'age':50., 'ca125':20.})
        captured = []
        class Scaler:
            n_features_in_ = 7
            def transform(self, values): captured.append(values.copy()); return values
        class Model:
            n_features_in_, classes_ = 7, [0, 1]
            score = .7
            def predict_proba(self, values): return np.array([[1-self.score,self.score]])
        model = Model()
        ns.update(argparse=argparse, DEFAULT_ARTIFACT_DIR=FakePath(), joblib=SimpleNamespace(load=lambda p: model if p.name == 'frozen_model.joblib' else Scaler()))
        args = ['synthetic']
        for feature, value in zip(SEVEN, [1,2,3,4,5,'nan',12]): args += ['--'+feature, str(value)]
        with patch('sys.argv', args), patch('sys.stdout', new=io.StringIO()) as out:
            self.assertEqual(ns['main'](), 0)
        result = json.loads(out.getvalue())
        self.assertTrue(result['research_only'])
        self.assertEqual(result['research_classification'], 'Above legacy research cutoff')
        np.testing.assert_allclose(captured[0], [[*np.log1p([1,2,3,4,5]),50.,12.]])
        for score in (float('nan'), float('inf'), -1, 2):
            model.score = score
            with patch('sys.argv', args), self.assertRaises(ValueError): ns['main']()

    def test_six_paths_compile_without_import_or_bytecode(self):
        for name in (CLI, PORTAL, 'demo_ovarian_cases.py', 'train_mirna_baseline.py'):
            compile(ast.parse(source(name)), name, 'exec')
        for name in ('Launch_Diagnostic_Portal.bat', 'Launch_Diagnostic_Portal.command'):
            self.assertTrue((ROOT / name).is_file())

    def test_immutable_artifact_hashes(self):
        import hashlib
        snapshot = json.loads((ROOT / 'audit/resume-verification/integrity/initial_snapshot.json').read_text())
        artifacts = {name: digest for name, digest in snapshot.items() if name.startswith('production_artifacts/')}
        self.assertGreaterEqual(len(artifacts), 4)
        for name, digest in artifacts.items():
            self.assertEqual(hashlib.sha256((ROOT / name).read_bytes()).hexdigest(), digest, name)

    @unittest.skipIf(np is None, 'numpy unavailable: model class validation runtime')
    def test_native_schema_and_fake_dimensions(self):
        for name, expected in ((CLI, SEVEN), (PORTAL, NINE)):
            validate = helpers(name, ['validate_artifact_schema'])['validate_artifact_schema']
            validate({'features': list(expected)}, SimpleNamespace(n_features_in_=len(expected)), SimpleNamespace(n_features_in_=len(expected), classes_=[0, 1]))
            for manifest in (None, {}, {'features': list(reversed(expected))}, {'features': list(expected[:-1])}, {'features': [expected[0]] * len(expected)}, {'features': list(expected[:-1]) + ['wrong']}):
                with self.subTest(name=name, manifest=manifest), self.assertRaises(ValueError):
                    validate(manifest)
            for dimension in (None, len(expected)-1, len(expected)+1):
                with self.assertRaises(ValueError):
                    validate({'features': list(expected)}, SimpleNamespace(n_features_in_=dimension))
                with self.assertRaises(ValueError):
                    validate({'features': list(expected)}, model=SimpleNamespace(n_features_in_=dimension))

    def test_cli_schema_precedes_joblib_without_artifact_reads(self):
        ns = helpers(CLI, ['main', 'parse_number', 'validate_artifact_schema'])
        reads = []
        class FakePath:
            def __truediv__(self, other): return self
            def read_text(self):
                reads.append('manifest')
                return json.dumps({'features': list(NINE)})
        ns.update(argparse=argparse, DEFAULT_ARTIFACT_DIR=FakePath(), joblib=SimpleNamespace(load=lambda *a: self.fail('joblib load forbidden')))
        args = ['synthetic']
        for feature in SEVEN: args += ['--' + feature, '1']
        with patch('sys.argv', args), self.assertRaises(ValueError): ns['main']()
        self.assertEqual(reads, ['manifest'])

    def test_baseline_blocked_before_io_and_help(self):
        ns = extract('train_mirna_baseline.py', ['main', 'require_approved_training_protocol'], dict(argparse=argparse, Path=Path))
        class Forbidden:
            def glob(self, *args): raise AssertionError('data IO forbidden')
        ns['RAW'] = Forbidden()
        with patch('sys.argv', ['synthetic']), self.assertRaisesRegex(RuntimeError, 'Training blocked'):
            ns['main']()
        with patch('sys.argv', ['synthetic', '--help']), patch('sys.stdout', new=io.StringIO()) as out, self.assertRaises(SystemExit) as exit:
            ns['main']()
        self.assertEqual(exit.exception.code, 0)
        self.assertIn('--cohort', out.getvalue())

    def test_no_top_level_mkdir(self):
        for name in ('train_mirna_baseline.py', CLI, PORTAL):
            tree = ast.parse(source(name))
            for node in tree.body:
                if isinstance(node, (ast.FunctionDef, ast.ClassDef, ast.If)): continue
                self.assertFalse(any(isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr == 'mkdir' for n in ast.walk(node)))

    def test_posix_launcher_raw_lf(self):
        raw = (ROOT / 'Launch_Diagnostic_Portal.command').read_bytes()
        self.assertNotIn(bytes([13]), raw, 'POSIX launcher must use LF, not CRLF')
        self.assertIn(bytes([10]), raw)
        self.assertTrue(raw.startswith(b'#!/'))

    def test_launchers_noninstall(self):
        for filename in ('Launch_Diagnostic_Portal.bat', 'Launch_Diagnostic_Portal.command'):
            text = (ROOT / filename).read_text().lower()
            self.assertNotRegex(text, r'pip\s+install')
            self.assertIn('gatherusagestats', text)
            self.assertIn('false', text)
            self.assertIn('import', text)

@unittest.skipIf(np is None, 'numpy unavailable')
class ClinicalTests(unittest.TestCase):
    def setUp(self):
        self.ns = helpers(PORTAL, ['effective_clinical_inputs', 'make_snapshot', 'snapshot_matches', 'build_model_input', 'validate_artifact_schema'])
        self.values = dict(zip(MIRNA, [0., 1., 2., 3., 4.]), age=42., ca125=12.)
        self.defaults = {'age': 50., 'ca125': 20.}

    def test_effective_defaults_and_rejections(self):
        fn = self.ns['effective_clinical_inputs']
        effective, flags = fn(dict(self.values, age=float('nan')), self.defaults)
        self.assertEqual(effective, {'age': 50., 'ca125': 12.})
        self.assertEqual(flags, {'age': True, 'ca125': False})
        for value in (float('inf'), -1, 'bad'):
            for key in ('age', 'ca125'):
                with self.assertRaises((ValueError, TypeError)): fn(dict(self.values, **{key: value}), self.defaults)
        for value in (float('nan'), float('inf'), -1, 'bad'):
            with self.assertRaises((ValueError, TypeError)): fn(self.values, dict(self.defaults, age=value))

    def test_session_cleared_on_failure_and_stale_snapshot(self):
        ns = helpers(PORTAL, ['main', 'clear_result', 'make_snapshot', 'snapshot_matches', 'effective_clinical_inputs'])
        class Stop(Exception): pass
        class FakeST:
            def __init__(self): self.session_state = {'probability': .9, 'research_result': {'snapshot': {'old': True}, 'probability': .9}}
            def error(self, *a): pass
            def stop(self): raise Stop()
            def markdown(self, *a, **k): pass
            def radio(self, *a, **k): raise Stop()
        for failure in (True, False):
            st = FakeST(); ns['st'] = st
            ns['artifact_identity'] = lambda: {'synthetic': 'identity'}
            def load():
                if failure: raise ValueError('synthetic artifact failure')
                return None, None, self.defaults, {'features': list(NINE)}
            ns['load_artifacts'] = load
            ns['sidebar_inputs'] = lambda defaults: self.values
            with self.assertRaises(Stop): ns['main']()
            self.assertNotIn('probability', st.session_state)
            self.assertNotIn('research_result', st.session_state)

    def test_main_mode_attestation_and_artifact_change(self):
        ns = helpers(PORTAL, ['main', 'clear_result', 'make_snapshot', 'snapshot_matches', 'effective_clinical_inputs'])
        class Box:
            def __enter__(self): return self
            def __exit__(self, *args): return False
        class FakeST:
            def __init__(self, mode, measured):
                self.mode, self.measured = mode, measured
                self.session_state = {'probability': .9, 'research_mode': 'Single patient', 'research_result': {'snapshot': ns['make_snapshot'](self_values, defaults, {'id': 'a'}), 'probability': .9}}
            def markdown(self, *args, **kwargs): pass
            def write(self, *args, **kwargs): pass
            def error(self, *args): pass
            def success(self, *args): pass
            def caption(self, *args): pass
            def radio(self, *args): return self.mode
            def columns(self, *args, **kwargs): return Box(), Box()
            def checkbox(self, *args, **kwargs): return self.measured
            def button(self, *args, **kwargs):
                self.disabled = kwargs['disabled']
                return False
            def stop(self): raise RuntimeError('synthetic stop')
        self_values, defaults = self.values, self.defaults
        ns['load_artifacts'] = lambda: (None, None, defaults, {'features': list(NINE)})
        ns['sidebar_inputs'] = lambda d: self_values
        ns['render_bulk_upload'] = lambda *args: None
        ns['render_result'] = lambda *args: self.fail('stale result rendered')
        for mode, measured, identity in (('Single patient', False, {'id':'a'}), ('Bulk upload', True, {'id':'a'}), ('Single patient', True, {'id':'b'})):
            st = FakeST(mode, measured); ns['st'] = st
            ns['artifact_identity'] = lambda: identity
            ns['main']()
            self.assertNotIn('research_result', st.session_state)
            self.assertNotIn('probability', st.session_state)
            if mode == 'Single patient': self.assertEqual(st.disabled, not measured)
        st = FakeST('Single patient', True); ns['st'] = st
        st.button = lambda *args, **kwargs: True
        ns['artifact_identity'] = lambda: {'id':'a'}
        def reject_prediction(*args): raise ValueError('synthetic prediction failure')
        ns['build_model_input'] = reject_prediction
        ns['main']()
        self.assertNotIn('research_result', st.session_state)
        st = FakeST('Single patient', True); ns['st'] = st
        identities = iter([{'id':'a'}, {'id':'b'}])
        ns['artifact_identity'] = lambda: next(identities)
        with self.assertRaisesRegex(RuntimeError, 'synthetic stop'): ns['main']()
        self.assertNotIn('research_result', st.session_state)

    def test_clear_result_removes_both_keys(self):
        st = SimpleNamespace(session_state={'probability': .2, 'research_result': {}, 'keep': 1})
        fn = extract(PORTAL, ['clear_result'], {'st': st})['clear_result']
        fn(); fn()
        self.assertEqual(st.session_state, {'keep': 1})

    def test_snapshot_detached_and_invalidated(self):
        make, matches = self.ns['make_snapshot'], self.ns['snapshot_matches']
        values = dict(self.values, age=float('nan'))
        identity = {'model': 'synthetic-a'}
        stored = make(values, self.defaults, identity)
        json.dumps(stored, allow_nan=False)
        self.assertTrue(matches(stored, make(values, self.defaults, identity)))
        for v,d,a in ((dict(values, ca125=13), self.defaults, identity), (values, dict(self.defaults, age=51), identity), (values, self.defaults, {'model': 'synthetic-b'})):
            self.assertFalse(matches(stored, make(v,d,a)))
        frozen = json.dumps(stored, sort_keys=True)
        values['ca125'] = 99; self.defaults['age'] = 99; identity['model'] = 'mutated'
        self.assertEqual(json.dumps(stored, sort_keys=True), frozen)

    @unittest.skipIf(pd is None, 'pandas unavailable: transform runtime not verified')
    def test_transform_parity_and_molecular_rejection(self):
        fn = self.ns['build_model_input']; manifest = {'features': list(NINE)}
        for age in (42., float('nan')):
            values = dict(self.values, age=age)
            result = fn(values, self.defaults, manifest)
            expected = [np.log1p(values[f]) for f in MIRNA] + [50. if math.isnan(age) else age, 12.]
            expected += [expected[0]*12., expected[1]*expected[5]]
            self.assertEqual(list(result.columns), list(NINE))
            np.testing.assert_allclose(result.iloc[0].to_numpy(), expected)
        for invalid in (-1, float('inf'), float('nan'), 'bad'):
            with self.assertRaises((ValueError, TypeError)): fn(dict(self.values, **{MIRNA[0]:invalid}), self.defaults, manifest)
        missing = dict(self.values); del missing[MIRNA[0]]
        with self.assertRaises((ValueError, KeyError)): fn(missing, self.defaults, manifest)

@unittest.skipIf(pd is None or np is None, 'pandas/numpy unavailable: table runtime not verified')
class TableTests(unittest.TestCase):
    def setUp(self):
        self.ns = helpers(PORTAL, ['read_uploaded_table','normalized_label','build_batch_model_input','build_model_input','effective_clinical_inputs','validate_artifact_schema','validate_uploaded_inputs'])
    def upload(self, text):
        stream = io.BytesIO(text.encode()); stream.name = 'synthetic.csv'
        return self.ns['read_uploaded_table'](stream)
    def test_rows_transpose_and_duplicates(self):
        header = 'Patient ID,' + ','.join(SEVEN)
        row = 'SYNTHETIC,1,2,3,4,5,42,12'
        table = self.upload(header+'\n'+row+'\n')
        self.assertEqual(len(table), 1)
        transposed = 'Feature,SYNTHETIC\n'+'\n'.join(f'{f},{v}' for f,v in zip(SEVEN,[1,2,3,4,5,42,12]))
        self.assertEqual(len(self.upload(transposed)),1)
        for text in (header+'\n'+row+'\n'+row, header+','+MIRNA[0]+'\n'+row+',1', transposed+'\n'+MIRNA[0]+',1'):
            with self.assertRaises(ValueError): self.upload(text)
    def test_batch_single_parity_and_missing_molecular(self):
        table = self.upload('Patient ID,' + ','.join(SEVEN) + '\nSYNTHETIC,1,2,3,4,5,,12\n')
        defaults, manifest = {'age':50., 'ca125':20.}, {'features':list(NINE)}
        batch = self.ns['build_batch_model_input'](table, defaults, manifest)
        single = self.ns['build_model_input'](table.iloc[0].to_dict(), defaults, manifest)
        pd.testing.assert_frame_equal(batch, single)
        for text in ('Patient ID,' + ','.join(SEVEN) + '\nSYNTHETIC,,2,3,4,5,42,12', 'Patient ID,' + ','.join(SEVEN[1:]) + '\nSYNTHETIC,2,3,4,5,42,12'):
            with self.assertRaises(ValueError): self.upload(text)

    def test_bulk_fake_ui_success_and_artifact_races(self):
        ns = helpers(PORTAL, ['render_bulk_upload', 'clear_result', 'read_uploaded_table', 'normalized_label', 'validate_uploaded_inputs', 'validate_artifact_schema', 'build_batch_model_input', 'build_model_input', 'effective_clinical_inputs', 'validate_probabilities'])
        identity = {name: 'synthetic-loaded-' + name for name in ('frozen_model.joblib', 'feature_scaler.joblib', 'manifest.json', 'imputation_defaults.json')}
        defaults, manifest = {'age':50., 'ca125':20.}, {'features':list(NINE)}
        for scenario in ('success', 'before', 'after'):
            with self.subTest(scenario=scenario):
                calls, frames, downloads, successes, errors = [], [], [], [], []
                stream = io.BytesIO(('Patient ID,' + ','.join(SEVEN) + chr(10) + 'SYNTHETIC,1,2,3,4,5,,12' + chr(10)).encode())
                stream.name = 'synthetic.csv'
                class UI:
                    session_state = {'probability':.99, 'research_result':{'probability':.99}}
                    def file_uploader(self, *a, **k): return stream
                    def dataframe(self, frame, **k): frames.append(frame.copy())
                    def download_button(self, label, data, *a, **k): downloads.append(data)
                    def success(self, text): successes.append(text)
                    def error(self, text): errors.append(text)
                    def __getattr__(self, name): return lambda *a, **k: None
                class Scaler:
                    n_features_in_ = 9
                    def transform(self, values): calls.append('scale'); return values
                class Model:
                    n_features_in_, classes_ = 9, [0,1]
                    def predict_proba(self, values): calls.append('predict'); return np.array([[.3,.7]])
                changed = dict(identity, manifest='synthetic-replacement')
                sequence = iter([changed] if scenario == 'before' else [identity, changed if scenario == 'after' else identity])
                ui = UI()
                ns.update(st=ui, artifact_identity=lambda: next(sequence))
                ns['render_bulk_upload'](Model(), Scaler(), defaults, manifest, identity)
                if scenario != 'success':
                    self.assertEqual(calls, [] if scenario == 'before' else ['scale','predict'])
                    self.assertTrue(errors)
                    self.assertEqual((frames,downloads,successes), ([],[],[]))
                    self.assertNotIn('probability', ui.session_state)
                    self.assertNotIn('research_result', ui.session_state)
                    continue
                self.assertEqual(errors, [])
                self.assertEqual(calls, ['scale','predict'])
                self.assertEqual((len(frames),len(downloads),len(successes)), (1,1,1))
                exported = pd.read_csv(io.BytesIO(downloads[0]))
                for frame in (frames[0], exported):
                    row = frame.iloc[0]
                    self.assertTrue(row['age_imputed']); self.assertFalse(row['ca125_imputed'])
                    self.assertEqual(row['age effective'],50.)
                    self.assertTrue(pd.isna(row['age raw measured']))
                    self.assertEqual(row['age input source'],'Frozen default')
                    self.assertEqual(row['ca125 effective'],12.)
                    self.assertEqual(row['ca125 raw measured'],12.)
                    self.assertEqual(row['ca125 input source'],'Measured')
                    self.assertEqual(row['legacy_research_cutoff'],.5)
                    self.assertEqual(row['threshold_selection_provenance'],'Unknown')
                    self.assertEqual(json.loads(row['artifact_identity_sha256']), identity)
                    self.assertTrue(row['research_only'])

    def test_xlsx_both_orientations(self):
        try:
            import openpyxl
        except ImportError:
            self.skipTest('openpyxl unavailable: XLSX runtime not verified')
        values = [1,2,3,4,5,42,12]
        rows = pd.DataFrame([['SYNTHETIC'] + values], columns=['Patient ID'] + list(SEVEN))
        transposed = pd.DataFrame({'Feature':list(SEVEN), 'SYNTHETIC':values})
        outputs = []
        for table in (rows, transposed):
            stream = io.BytesIO(); stream.name = 'synthetic.xlsx'
            table.to_excel(stream, index=False, engine='openpyxl'); stream.seek(0)
            outputs.append(self.ns['read_uploaded_table'](stream))
        pd.testing.assert_frame_equal(outputs[0], outputs[1])

    def test_invalid_upload_ids(self):
        for ids in ([''], [None], ['x', ' x '], ['a-b', 'ab']):
            table = pd.DataFrame([dict(zip(SEVEN, [1]*7), **{'Patient ID': value}) for value in ids])
            with self.assertRaises(ValueError): self.ns['validate_uploaded_inputs'](table)

    def test_malformed_batch_cannot_pass(self):
        for key in SEVEN:
            for bad in ('bad', '-1', 'inf'):
                values = dict(zip(SEVEN, ['1']*7)); values[key] = bad
                text = 'Patient ID,'+','.join(SEVEN)+'\nSYNTHETIC,'+','.join(values[k] for k in SEVEN)
                with self.assertRaises((ValueError, TypeError)):
                    table = self.upload(text)
                    self.ns['build_batch_model_input'](table, {'age':50,'ca125':20}, {'features':list(NINE)})

    def test_demo_lookup(self):
        fn = helpers('demo_ovarian_cases.py',['build_case_lookup'])['build_case_lookup']
        self.assertEqual(set(fn(pd.DataFrame([{'sample_id':' SYNTHETIC '}]))), {'SYNTHETIC'})
        for ids in ([''],[None],['x',' x '],[1,'1']):
            with self.assertRaises(ValueError): fn(pd.DataFrame({'sample_id':ids}))

if __name__ == '__main__': unittest.main()
