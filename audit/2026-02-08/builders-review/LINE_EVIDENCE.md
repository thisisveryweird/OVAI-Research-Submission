# Static numbered source evidence

Complete numbered sources, read-only snapshots for precise evidence. Hashes in COVERAGE.csv bind originals.

## scripts/build_14mirna_benchmark.py

~~~python
1: #!/usr/bin/env python3
2: 
3: """
4: Create the exact 14-miRNA benchmark matrix used by the historical
5: eLife ovarian-cancer neural-network model.
6: 
7: This script DOES NOT train a model.
8: 
9: It:
10: 
11: 1. Loads the frozen 179-patient GSE94533 dataset.
12: 2. Uses the exact 14 miRNAs recovered from the original PMML.
13: 3. Verifies every feature exists exactly once.
14: 4. Joins diagnosis labels.
15: 5. Writes a compact patient-level benchmark matrix.
16: 6. Produces detailed provenance/audit information.
17: 
18: The goal is to create a reproducible starting point for:
19:     historical reproduction
20:         ->
21:     modern baseline
22:         ->
23:     multimodal model
24: """
25: 
26: from __future__ import annotations
27: 
28: import json
29: from pathlib import Path
30: 
31: import numpy as np
32: import pandas as pd
33: 
34: 
35: # ================================================================
36: # PATHS
37: # ================================================================
38: 
39: ROOT = Path.home() / "ovarian-ai"
40: 
41: MATRIX_PATH = (
42:     ROOT
43:     / "data"
44:     / "processed"
45:     / "gse94533_mirna_tpm.csv.gz"
46: )
47: 
48: METADATA_PATH = (
49:     ROOT
50:     / "data"
51:     / "processed"
52:     / "gse94533_metadata.csv"
53: )
54: 
55: FEATURE_JSON = (
56:     ROOT
57:     / "reports"
58:     / "tables"
59:     / "original_14_mirna_features.json"
60: )
61: 
62: OUTPUT_DIR = (
63:     ROOT
64:     / "data"
65:     / "processed"
66: )
67: 
68: REPORT_DIR = (
69:     ROOT
70:     / "reports"
71:     / "tables"
72: )
73: 
74: OUTPUT_MATRIX = (
75:     OUTPUT_DIR
76:     / "gse94533_14mirna_benchmark.csv.gz"
77: )
78: 
79: OUTPUT_AUDIT = (
80:     REPORT_DIR
81:     / "gse94533_14mirna_benchmark_audit.json"
82: )
83: 
84: 
85: # ================================================================
86: # FALLBACK EXACT FEATURE LIST
87: # ================================================================
88: #
89: # This is intentionally written explicitly as a second safety
90: # mechanism. The PMML-derived JSON remains the provenance source,
91: # but these are the exact 14 inputs we already verified.
92: #
93: # The script will compare the JSON against this list.
94: #
95: 
96: EXPECTED_14 = [
97:     "hsa-miR-23b-3p",
98:     "hsa-miR-29a-3p",
99:     "hsa-miR-32-5p",
100:     "hsa-miR-92a-3p",
101:     "hsa-miR-150-5p",
102:     "hsa-miR-200a-3p",
103:     "hsa-miR-200c-3p",
104:     "hsa-miR-203a",
105:     "hsa-miR-320c",
106:     "hsa-miR-320d",
107:     "hsa-miR-335-5p",
108:     "hsa-miR-450b-5p",
109:     "hsa-miR-1246",
110:     "hsa-miR-1307-5p",
111: ]
112: 
113: 
114: # ================================================================
115: # LOAD FEATURE PROVENANCE
116: # ================================================================
117: 
118: def load_pmml_features() -> list[str]:
119: 
120:     if not FEATURE_JSON.exists():
121: 
122:         raise RuntimeError(
123:             "Could not find:\n"
124:             f"{FEATURE_JSON}\n\n"
125:             "Run extract_original_14mirna_model.py first."
126:         )
127: 
128:     data = json.loads(
129:         FEATURE_JSON.read_text()
130:     )
131: 
132:     fields = data.get(
133:         "mirna_named_fields",
134:         [],
135:     )
136: 
137:     fields = list(
138:         dict.fromkeys(
139:             str(x).strip()
140:             for x in fields
141:         )
142:     )
143: 
144:     return fields
145: 
146: 
147: # ================================================================
148: # LOAD DATA
149: # ================================================================
150: 
151: def load_matrix() -> pd.DataFrame:
152: 
153:     if not MATRIX_PATH.exists():
154: 
155:         raise RuntimeError(
156:             f"Missing expression matrix:\n{MATRIX_PATH}"
157:         )
158: 
159:     matrix = pd.read_csv(
160:         MATRIX_PATH,
161:         index_col=0,
162:     )
163: 
164:     return matrix
165: 
166: 
167: def load_metadata() -> pd.DataFrame:
168: 
169:     if not METADATA_PATH.exists():
170: 
171:         raise RuntimeError(
172:             f"Missing metadata:\n{METADATA_PATH}"
173:         )
174: 
175:     metadata = pd.read_csv(
176:         METADATA_PATH,
177:     )
178: 
179:     required = {
180:         "sample_id_normalized",
181:         "diagnosis",
182:         "binary_label",
183:     }
184: 
185:     missing = required - set(
186:         metadata.columns
187:     )
188: 
189:     if missing:
190: 
191:         raise RuntimeError(
192:             f"Metadata missing columns: {missing}"
193:         )
194: 
195:     return metadata
196: 
197: 
198: # ================================================================
199: # MAIN BUILD
200: # ================================================================
201: 
202: def main() -> None:
203: 
204:     print("=" * 90)
205:     print("BUILDING EXACT 14-miRNA BENCHMARK DATASET")
206:     print("=" * 90)
207: 
208:     OUTPUT_DIR.mkdir(
209:         parents=True,
210:         exist_ok=True,
211:     )
212: 
213:     REPORT_DIR.mkdir(
214:         parents=True,
215:         exist_ok=True,
216:     )
217: 
218:     # ------------------------------------------------------------
219:     # Load PMML-derived features.
220:     # ------------------------------------------------------------
221: 
222:     pmml_features = load_pmml_features()
223: 
224:     print()
225:     print("Features recovered from original PMML:")
226: 
227:     for feature in pmml_features:
228:         print(
229:             f"  {feature}"
230:         )
231: 
232:     # ------------------------------------------------------------
233:     # Require exactly the known 14 inputs.
234:     # ------------------------------------------------------------
235: 
236:     if set(pmml_features) != set(EXPECTED_14):
237: 
238:         print()
239:         print("EXPECTED:")
240:         for x in EXPECTED_14:
241:             print(" ", x)
242: 
243:         print()
244:         print("PMML:")
245:         for x in pmml_features:
246:             print(" ", x)
247: 
248:         raise RuntimeError(
249:             "PMML-derived feature set does not match the expected "
250:             "14-feature historical model."
251:         )
252: 
253:     if len(pmml_features) != 14:
254: 
255:         raise RuntimeError(
256:             f"Expected 14 features, found {len(pmml_features)}."
257:         )
258: 
259:     # ------------------------------------------------------------
260:     # Load patient expression data.
261:     # ------------------------------------------------------------
262: 
263:     matrix = load_matrix()
264: 
265:     print()
266:     print(
267:         f"Loaded expression matrix: "
268:         f"{matrix.shape[0]} patients x "
269:         f"{matrix.shape[1]} features"
270:     )
271: 
272:     # ------------------------------------------------------------
273:     # Confirm all 14 features exist.
274:     # ------------------------------------------------------------
275: 
276:     missing = [
277:         feature
278:         for feature in pmml_features
279:         if feature not in matrix.columns
280:     ]
281: 
282:     if missing:
283: 
284:         raise RuntimeError(
285:             "Historical model features missing from cleaned dataset:\n"
286:             + "\n".join(
287:                 f"  {x}"
288:                 for x in missing
289:             )
290:         )
291: 
292:     duplicate_features = [
293:         feature
294:         for feature in pmml_features
295:         if list(matrix.columns).count(feature) != 1
296:     ]
297: 
298:     if duplicate_features:
299: 
300:         raise RuntimeError(
301:             "Historical features are duplicated in expression matrix:\n"
302:             + "\n".join(
303:                 f"  {x}"
304:                 for x in duplicate_features
305:             )
306:         )
307: 
308:     # ------------------------------------------------------------
309:     # Extract exact 14 columns IN PMML ORDER.
310:     # ------------------------------------------------------------
311: 
312:     benchmark = matrix[
313:         pmml_features
314:     ].copy()
315: 
316:     # ------------------------------------------------------------
317:     # Load metadata.
318:     # ------------------------------------------------------------
319: 
320:     metadata = load_metadata()
321: 
322:     if len(metadata) != len(matrix):
323: 
324:         raise RuntimeError(
325:             "Metadata patient count does not match expression matrix."
326:         )
327: 
328:     # ------------------------------------------------------------
329:     # Verify patient alignment.
330:     # ------------------------------------------------------------
331: 
332:     metadata_ids = (
333:         metadata[
334:             "sample_id_normalized"
335:         ]
336:         .astype(str)
337:         .tolist()
338:     )
339: 
340:     matrix_ids = [
341:         str(x)
342:         for x in matrix.index
343:     ]
344: 
345:     if metadata_ids != matrix_ids:
346: 
347:         raise RuntimeError(
348:             "Patient ordering differs between expression matrix "
349:             "and metadata."
350:         )
351: 
352:     # ------------------------------------------------------------
353:     # Add metadata to benchmark.
354:     # ------------------------------------------------------------
355: 
356:     benchmark.insert(
357:         0,
358:         "binary_label",
359:         metadata["binary_label"].astype(int).values,
360:     )
361: 
362:     benchmark.insert(
363:         0,
364:         "diagnosis",
365:         metadata["diagnosis"].astype(str).values,
366:     )
367: 
368:     benchmark.insert(
369:         0,
370:         "sample_id",
371:         metadata[
372:             "sample_id_normalized"
373:         ].astype(str).values,
374:     )
375: 
376:     # Add stage/histology where available.
377:     if "stage" in metadata.columns:
378: 
379:         benchmark["stage"] = (
380:             metadata["stage"]
381:             .astype("string")
382:             .values
383:         )
384: 
385:     if "histology" in metadata.columns:
386: 
387:         benchmark["histology"] = (
388:             metadata["histology"]
389:             .astype("string")
390:             .values
391:         )
392: 
393:     # ------------------------------------------------------------
394:     # Numeric checks.
395:     # ------------------------------------------------------------
396: 
397:     feature_values = benchmark[
398:         pmml_features
399:     ].to_numpy(
400:         dtype=float
401:     )
402: 
403:     if not np.isfinite(
404:         feature_values
405:     ).all():
406: 
407:         raise RuntimeError(
408:             "14-miRNA benchmark contains NaN or infinite values."
409:         )
410: 
411:     if (
412:         feature_values < 0
413:     ).any():
414: 
415:         raise RuntimeError(
416:             "14-miRNA TPM matrix contains negative values."
417:         )
418: 
419:     # ------------------------------------------------------------
420:     # Confirm patient uniqueness.
421:     # ------------------------------------------------------------
422: 
423:     if benchmark[
424:         "sample_id"
425:     ].duplicated().any():
426: 
427:         raise RuntimeError(
428:             "Duplicate patient IDs found."
429:         )
430: 
431:     if len(benchmark) != 179:
432: 
433:         raise RuntimeError(
434:             f"Expected 179 patients, found {len(benchmark)}."
435:         )
436: 
437:     # ------------------------------------------------------------
438:     # Summary.
439:     # ------------------------------------------------------------
440: 
441:     diagnosis_counts = (
442:         benchmark[
443:             "diagnosis"
444:         ]
445:         .value_counts()
446:         .to_dict()
447:     )
448: 
449:     binary_counts = (
450:         benchmark[
451:             "binary_label"
452:         ]
453:         .value_counts()
454:         .sort_index()
455:         .to_dict()
456:     )
457: 
458:     audit = {
459: 
460:         "dataset": "GSE94533",
461: 
462:         "purpose": (
463:             "Exact 14-miRNA feature matrix corresponding to "
464:             "the historical eLife ovarian-cancer neural-network model."
465:         ),
466: 
467:         "n_patients": int(
468:             len(benchmark)
469:         ),
470: 
471:         "n_features": 14,
472: 
473:         "features": pmml_features,
474: 
475:         "original_model_architecture": {
476:             "input_features": 14,
477:             "hidden_neurons": 7,
478:             "output_classes": 2,
479:             "model_name": (
480:                 "Dataset complete with classificat_MLP 14-7-2"
481:             ),
482:         },
483: 
484:         "feature_provenance": str(
485:             FEATURE_JSON
486:         ),
487: 
488:         "source_expression_matrix": str(
489:             MATRIX_PATH
490:         ),
491: 
492:         "source_metadata": str(
493:             METADATA_PATH
494:         ),
495: 
496:         "diagnosis_counts": diagnosis_counts,
497: 
498:         "binary_counts": binary_counts,
499: 
500:         "binary_endpoint": {
501:             "1": "Cancer",
502:             "0": "Benign, Borderline, or Control",
503:         },
504: 
505:         "validation_checks": {
506:             "patient_count_is_179": (
507:                 len(benchmark) == 179
508:             ),
509:             "exactly_14_features": (
510:                 len(pmml_features) == 14
511:             ),
512:             "all_features_present": (
513:                 len(missing) == 0
514:             ),
515:             "all_values_finite": bool(
516:                 np.isfinite(
517:                     feature_values
518:                 ).all()
519:             ),
520:             "all_values_nonnegative": bool(
521:                 (
522:                     feature_values >= 0
523:                 ).all()
524:             ),
525:             "unique_patient_ids": (
526:                 not benchmark[
527:                     "sample_id"
528:                 ].duplicated().any()
529:             ),
530:             "patient_order_aligned": True,
531:         },
532: 
533:         "important_note": (
534:             "The labels in this benchmark come from the current "
535:             "GSE94533 GEO metadata used to construct the cleaned dataset. "
536:             "They are intentionally not altered to force agreement with "
537:             "historical class counts."
538:         ),
539:     }
540: 
541:     # ------------------------------------------------------------
542:     # Save compact matrix.
543:     # ------------------------------------------------------------
544: 
545:     benchmark.to_csv(
546:         OUTPUT_MATRIX,
547:         index=False,
548:         compression="gzip",
549:     )
550: 
551:     OUTPUT_AUDIT.write_text(
552:         json.dumps(
553:             audit,
554:             indent=2,
555:         )
556:     )
557: 
558:     # ------------------------------------------------------------
559:     # Terminal report.
560:     # ------------------------------------------------------------
561: 
562:     print()
563:     print("=" * 90)
564:     print("14-miRNA BENCHMARK DATASET")
565:     print("=" * 90)
566: 
567:     print()
568:     print(
569:         f"Patients: {len(benchmark)}"
570:     )
571: 
572:     print(
573:         f"Features: {len(pmml_features)}"
574:     )
575: 
576:     print()
577:     print("Feature order:")
578: 
579:     for i, feature in enumerate(
580:         pmml_features,
581:         1,
582:     ):
583: 
584:         print(
585:             f"  {i:2d}. {feature}"
586:         )
587: 
588:     print()
589:     print("Diagnosis counts:")
590: 
591:     for diagnosis, count in diagnosis_counts.items():
592: 
593:         print(
594:             f"  {diagnosis:12s} {count}"
595:         )
596: 
597:     print()
598:     print("Binary counts:")
599: 
600:     print(
601:         f"  0 = non-cancer   "
602:         f"{binary_counts.get(0, 0)}"
603:     )
604: 
605:     print(
606:         f"  1 = cancer       "
607:         f"{binary_counts.get(1, 0)}"
608:     )
609: 
610:     print()
611:     print("Validation:")
612: 
613:     for name, result in audit[
614:         "validation_checks"
615:     ].items():
616: 
617:         print(
618:             f"  {name:30s} {result}"
619:         )
620: 
621:     print()
622:     print("=" * 90)
623:     print("FILES CREATED")
624:     print("=" * 90)
625: 
626:     print()
627:     print(
628:         f"Benchmark matrix:\n"
629:         f"  {OUTPUT_MATRIX}"
630:     )
631: 
632:     print()
633:     print(
634:         f"Audit:\n"
635:         f"  {OUTPUT_AUDIT}"
636:     )
637: 
638:     print()
639:     print("=" * 90)
640:     print("NO MODEL TRAINING PERFORMED")
641:     print("=" * 90)
642: 
643:     print(
644:         """
645: We now have a frozen 14-feature patient-level dataset.
646: 
647: The next scientific step is to reproduce/benchmark the historical
648: classifier and establish a modern baseline using proper patient-level
649: cross-validation.
650: 
651: We should NOT immediately train a large neural network on all 2,582
652: features because 179 patients versus thousands of features creates a
653: serious overfitting risk.
654: """
655:     )
656: 
657: 
658: if __name__ == "__main__":
659:     main()
~~~

## scripts/build_age_192mirna_model.py

~~~python
1: #!/usr/bin/env python3
2: 
3: from pathlib import Path
4: import json
5: import warnings
6: 
7: import numpy as np
8: import pandas as pd
9: 
10: from sklearn.impute import SimpleImputer
11: from sklearn.linear_model import LogisticRegression
12: from sklearn.metrics import (
13:     accuracy_score,
14:     average_precision_score,
15:     confusion_matrix,
16:     roc_auc_score,
17: )
18: from sklearn.model_selection import StratifiedKFold
19: from sklearn.pipeline import Pipeline
20: from sklearn.preprocessing import StandardScaler
21: 
22: warnings.filterwarnings("ignore")
23: 
24: ROOT = Path.home() / "ovarian-ai"
25: 
26: DATA = (
27:     ROOT
28:     / "data"
29:     / "processed"
30:     / "gse94533_14mirna_clinical_benchmark.csv.gz"
31: )
32: 
33: OUTDIR = ROOT / "reports" / "tables"
34: 
35: 
36: def main():
37: 
38:     print("=" * 80)
39:     print("AGE + 192-miRNA MODEL")
40:     print("=" * 80)
41: 
42:     df = pd.read_csv(DATA)
43: 
44:     mirnas = [
45:         c for c in df.columns
46:         if c.startswith("hsa-")
47:     ]
48: 
49:     if len(mirnas) != 14:
50:         raise RuntimeError(
51:             f"Expected 14 miRNA columns in this benchmark, found {len(mirnas)}."
52:         )
53: 
54:     # The current clinical benchmark contains only the exact 14 PMML
55:     # miRNAs. Build the full 192-miRNA matrix from Dataset 2.
56:     xlsx = (
57:         ROOT
58:         / "data"
59:         / "external"
60:         / "elife_28932"
61:         / "elife_poa_e28932_Supplementary_File_6.xlsx"
62:     )
63: 
64:     raw = pd.read_excel(
65:         xlsx,
66:         sheet_name="Supplemental Dataset 2",
67:         header=None,
68:         engine="openpyxl",
69:     )
70: 
71:     feature_names = (
72:         raw.iloc[:, 0]
73:         .astype(str)
74:         .str.strip()
75:     )
76: 
77:     all_features = [
78:         x for x in feature_names
79:         if x.startswith("hsa-")
80:     ]
81: 
82:     # Remove duplicates while preserving order.
83:     all_features = list(dict.fromkeys(all_features))
84: 
85:     if len(all_features) != 192:
86:         raise RuntimeError(
87:             f"Expected 192 miRNA features, found {len(all_features)}."
88:         )
89: 
90:     # The verified PMML cohort establishes the exact 179 column order.
91:     pmml = pd.read_csv(
92:         ROOT
93:         / "reports"
94:         / "tables"
95:         / "original_historical_clean_predictions.csv"
96:     )
97: 
98:     verified_ids = pmml["sample_id"].astype(str).tolist()
99: 
100:     if len(verified_ids) != 179:
101:         raise RuntimeError(
102:             "Verified PMML cohort is not 179 patients."
103:         )
104: 
105:     matrix = {}
106: 
107:     for feature in all_features:
108: 
109:         idx = np.where(
110:             feature_names.values == feature
111:         )[0]
112: 
113:         if len(idx) != 1:
114:             raise RuntimeError(
115:                 f"Could not uniquely locate {feature}."
116:             )
117: 
118:         values = pd.to_numeric(
119:             raw.iloc[idx[0], 1:180],
120:             errors="coerce",
121:         ).to_numpy(float)
122: 
123:         if len(values) != 179:
124:             raise RuntimeError(
125:                 f"{feature}: expected 179 values."
126:             )
127: 
128:         matrix[feature] = values
129: 
130:     Xmirna = pd.DataFrame(matrix)
131:     Xmirna["sample_id"] = verified_ids
132: 
133:     # Join age + verified label.
134:     clinical = df[
135:         ["sample_id", "age", "binary_label"]
136:     ].copy()
137: 
138:     merged = Xmirna.merge(
139:         clinical,
140:         on="sample_id",
141:         how="inner",
142:         validate="one_to_one",
143:     )
144: 
145:     if len(merged) != 179:
146:         raise RuntimeError(
147:             f"Expected 179 merged patients, got {len(merged)}."
148:         )
149: 
150:     features = ["age"] + all_features
151: 
152:     X = merged[features]
153:     y = merged["binary_label"].astype(int)
154: 
155:     print()
156:     print("Patients:", len(merged))
157:     print("Features:", len(features))
158:     print("miRNAs:", len(all_features))
159:     print("Cancer:", int(y.sum()))
160:     print("Non-cancer:", int((y == 0).sum()))
161: 
162:     cv = StratifiedKFold(
163:         n_splits=5,
164:         shuffle=True,
165:         random_state=42,
166:     )
167: 
168:     pipe = Pipeline([
169:         (
170:             "imputer",
171:             SimpleImputer(strategy="median"),
172:         ),
173:         (
174:             "scaler",
175:             StandardScaler(),
176:         ),
177:         (
178:             "model",
179:             LogisticRegression(
180:                 C=0.3,
181:                 penalty="l2",
182:                 solver="liblinear",
183:                 max_iter=5000,
184:             ),
185:         ),
186:     ])
187: 
188:     oof_prob = np.full(len(y), np.nan)
189: 
190:     for fold, (train_idx, test_idx) in enumerate(
191:         cv.split(X, y),
192:         start=1,
193:     ):
194: 
195:         pipe.fit(
196:             X.iloc[train_idx],
197:             y.iloc[train_idx],
198:         )
199: 
200:         oof_prob[test_idx] = pipe.predict_proba(
201:             X.iloc[test_idx]
202:         )[:, 1]
203: 
204:     oof_pred = (
205:         oof_prob >= 0.5
206:     ).astype(int)
207: 
208:     auc = roc_auc_score(y, oof_prob)
209:     pr_auc = average_precision_score(y, oof_prob)
210:     accuracy = accuracy_score(y, oof_pred)
211: 
212:     tn, fp, fn, tp = confusion_matrix(
213:         y,
214:         oof_pred,
215:     ).ravel()
216: 
217:     sensitivity = tp / (tp + fn)
218:     specificity = tn / (tn + fp)
219: 
220:     print()
221:     print("=" * 80)
222:     print("RESULTS")
223:     print("=" * 80)
224: 
225:     print(f"OOF ROC AUC:    {auc:.6f}")
226:     print(f"OOF PR AUC:     {pr_auc:.6f}")
227:     print(f"OOF Accuracy:   {accuracy:.6f}")
228:     print(f"Sensitivity:    {sensitivity:.6f}")
229:     print(f"Specificity:    {specificity:.6f}")
230: 
231:     print()
232:     print("Confusion matrix:")
233:     print(f"  TP={tp} FN={fn}")
234:     print(f"  FP={fp} TN={tn}")
235: 
236:     predictions = pd.DataFrame({
237:         "sample_id": merged["sample_id"],
238:         "true_label": y,
239:         "oof_probability": oof_prob,
240:         "oof_prediction": oof_pred,
241:     })
242: 
243:     pred_path = (
244:         OUTDIR
245:         / "age_plus_192mirna_oof_predictions.csv"
246:     )
247: 
248:     result_path = (
249:         OUTDIR
250:         / "age_plus_192mirna_results.json"
251:     )
252: 
253:     predictions.to_csv(
254:         pred_path,
255:         index=False,
256:     )
257: 
258:     result_path.write_text(
259:         json.dumps(
260:             {
261:                 "model": "Age + 192 miRNA",
262:                 "patients": len(merged),
263:                 "miRNA_features": 192,
264:                 "roc_auc": float(auc),
265:                 "pr_auc": float(pr_auc),
266:                 "accuracy": float(accuracy),
267:                 "sensitivity": float(sensitivity),
268:                 "specificity": float(specificity),
269:                 "tp": int(tp),
270:                 "fn": int(fn),
271:                 "fp": int(fp),
272:                 "tn": int(tn),
273:             },
274:             indent=2,
275:         )
276:     )
277: 
278:     print()
279:     print("Saved:")
280:     print(pred_path)
281:     print(result_path)
282: 
283: 
284: if __name__ == "__main__":
285:     main()
~~~

## scripts/build_age_ca125_192mirna_model.py

~~~python
1: #!/usr/bin/env python3
2: 
3: from pathlib import Path
4: import json
5: import warnings
6: 
7: import numpy as np
8: import pandas as pd
9: 
10: from sklearn.impute import SimpleImputer
11: from sklearn.linear_model import LogisticRegression
12: from sklearn.metrics import (
13:     accuracy_score,
14:     average_precision_score,
15:     confusion_matrix,
16:     roc_auc_score,
17: )
18: from sklearn.model_selection import StratifiedKFold
19: from sklearn.pipeline import Pipeline
20: from sklearn.preprocessing import StandardScaler
21: 
22: warnings.filterwarnings("ignore")
23: 
24: ROOT = Path.home() / "ovarian-ai"
25: 
26: DATA = (
27:     ROOT
28:     / "data"
29:     / "processed"
30:     / "gse94533_14mirna_clinical_benchmark.csv.gz"
31: )
32: 
33: OUTDIR = ROOT / "reports" / "tables"
34: 
35: 
36: def main():
37: 
38:     print("=" * 80)
39:     print("AGE + CA-125 + 192-miRNA MODEL")
40:     print("=" * 80)
41: 
42:     df = pd.read_csv(DATA)
43: 
44:     mirnas = [
45:         c for c in df.columns
46:         if c.startswith("hsa-")
47:     ]
48: 
49:     if len(mirnas) != 14:
50:         raise RuntimeError(
51:             f"Expected 14 miRNA columns in this benchmark, found {len(mirnas)}."
52:         )
53: 
54:     # The current clinical benchmark contains only the exact 14 PMML
55:     # miRNAs. Build the full 192-miRNA matrix from Dataset 2.
56:     xlsx = (
57:         ROOT
58:         / "data"
59:         / "external"
60:         / "elife_28932"
61:         / "elife_poa_e28932_Supplementary_File_6.xlsx"
62:     )
63: 
64:     raw = pd.read_excel(
65:         xlsx,
66:         sheet_name="Supplemental Dataset 2",
67:         header=None,
68:         engine="openpyxl",
69:     )
70: 
71:     feature_names = (
72:         raw.iloc[:, 0]
73:         .astype(str)
74:         .str.strip()
75:     )
76: 
77:     all_features = [
78:         x for x in feature_names
79:         if x.startswith("hsa-")
80:     ]
81: 
82:     # Remove duplicates while preserving order.
83:     all_features = list(dict.fromkeys(all_features))
84: 
85:     if len(all_features) != 192:
86:         raise RuntimeError(
87:             f"Expected 192 miRNA features, found {len(all_features)}."
88:         )
89: 
90:     # The verified PMML cohort establishes the exact 179 column order.
91:     pmml = pd.read_csv(
92:         ROOT
93:         / "reports"
94:         / "tables"
95:         / "original_historical_clean_predictions.csv"
96:     )
97: 
98:     verified_ids = pmml["sample_id"].astype(str).tolist()
99: 
100:     if len(verified_ids) != 179:
101:         raise RuntimeError(
102:             "Verified PMML cohort is not 179 patients."
103:         )
104: 
105:     matrix = {}
106: 
107:     for feature in all_features:
108: 
109:         idx = np.where(
110:             feature_names.values == feature
111:         )[0]
112: 
113:         if len(idx) != 1:
114:             raise RuntimeError(
115:                 f"Could not uniquely locate {feature}."
116:             )
117: 
118:         values = pd.to_numeric(
119:             raw.iloc[idx[0], 1:180],
120:             errors="coerce",
121:         ).to_numpy(float)
122: 
123:         if len(values) != 179:
124:             raise RuntimeError(
125:                 f"{feature}: expected 179 values."
126:             )
127: 
128:         matrix[feature] = values
129: 
130:     Xmirna = pd.DataFrame(matrix)
131:     Xmirna["sample_id"] = verified_ids
132: 
133:     # Join age + verified label.
134:     clinical = df[
135:         ["sample_id", "age", "ca125", "binary_label"]
136:     ].copy()
137: 
138:     merged = Xmirna.merge(
139:         clinical,
140:         on="sample_id",
141:         how="inner",
142:         validate="one_to_one",
143:     )
144: 
145:     if len(merged) != 179:
146:         raise RuntimeError(
147:             f"Expected 179 merged patients, got {len(merged)}."
148:         )
149: 
150:     features = ["age", "ca125"] + all_features
151: 
152:     X = merged[features]
153:     y = merged["binary_label"].astype(int)
154: 
155:     print()
156:     print("Patients:", len(merged))
157:     print("Features:", len(features))
158:     print("miRNAs:", len(all_features))
159:     print("Cancer:", int(y.sum()))
160:     print("Non-cancer:", int((y == 0).sum()))
161: 
162:     cv = StratifiedKFold(
163:         n_splits=5,
164:         shuffle=True,
165:         random_state=42,
166:     )
167: 
168:     pipe = Pipeline([
169:         (
170:             "imputer",
171:             SimpleImputer(strategy="median"),
172:         ),
173:         (
174:             "scaler",
175:             StandardScaler(),
176:         ),
177:         (
178:             "model",
179:             LogisticRegression(
180:                 C=0.3,
181:                 penalty="l2",
182:                 solver="liblinear",
183:                 max_iter=5000,
184:             ),
185:         ),
186:     ])
187: 
188:     oof_prob = np.full(len(y), np.nan)
189: 
190:     for fold, (train_idx, test_idx) in enumerate(
191:         cv.split(X, y),
192:         start=1,
193:     ):
194: 
195:         pipe.fit(
196:             X.iloc[train_idx],
197:             y.iloc[train_idx],
198:         )
199: 
200:         oof_prob[test_idx] = pipe.predict_proba(
201:             X.iloc[test_idx]
202:         )[:, 1]
203: 
204:     oof_pred = (
205:         oof_prob >= 0.5
206:     ).astype(int)
207: 
208:     auc = roc_auc_score(y, oof_prob)
209:     pr_auc = average_precision_score(y, oof_prob)
210:     accuracy = accuracy_score(y, oof_pred)
211: 
212:     tn, fp, fn, tp = confusion_matrix(
213:         y,
214:         oof_pred,
215:     ).ravel()
216: 
217:     sensitivity = tp / (tp + fn)
218:     specificity = tn / (tn + fp)
219: 
220:     print()
221:     print("=" * 80)
222:     print("RESULTS")
223:     print("=" * 80)
224: 
225:     print(f"OOF ROC AUC:    {auc:.6f}")
226:     print(f"OOF PR AUC:     {pr_auc:.6f}")
227:     print(f"OOF Accuracy:   {accuracy:.6f}")
228:     print(f"Sensitivity:    {sensitivity:.6f}")
229:     print(f"Specificity:    {specificity:.6f}")
230: 
231:     print()
232:     print("Confusion matrix:")
233:     print(f"  TP={tp} FN={fn}")
234:     print(f"  FP={fp} TN={tn}")
235: 
236:     predictions = pd.DataFrame({
237:         "sample_id": merged["sample_id"],
238:         "true_label": y,
239:         "oof_probability": oof_prob,
240:         "oof_prediction": oof_pred,
241:     })
242: 
243:     pred_path = (
244:         OUTDIR
245:         / "age_ca125_plus_192mirna_oof_predictions.csv"
246:     )
247: 
248:     result_path = (
249:         OUTDIR
250:         / "age_ca125_plus_192mirna_results.json"
251:     )
252: 
253:     predictions.to_csv(
254:         pred_path,
255:         index=False,
256:     )
257: 
258:     result_path.write_text(
259:         json.dumps(
260:             {
261:                 "model": "Age + CA-125 + 192 miRNA",
262:                 "patients": len(merged),
263:                 "miRNA_features": 192,
264:                 "roc_auc": float(auc),
265:                 "pr_auc": float(pr_auc),
266:                 "accuracy": float(accuracy),
267:                 "sensitivity": float(sensitivity),
268:                 "specificity": float(specificity),
269:                 "tp": int(tp),
270:                 "fn": int(fn),
271:                 "fp": int(fp),
272:                 "tn": int(tn),
273:             },
274:             indent=2,
275:         )
276:     )
277: 
278:     print()
279:     print("Saved:")
280:     print(pred_path)
281:     print(result_path)
282: 
283: 
284: if __name__ == "__main__":
285:     main()
~~~

## scripts/build_age_clinical_models.py

~~~python
1: #!/usr/bin/env python3
2: 
3: from pathlib import Path
4: import json
5: import warnings
6: 
7: import numpy as np
8: import pandas as pd
9: 
10: from sklearn.compose import ColumnTransformer
11: from sklearn.impute import SimpleImputer
12: from sklearn.linear_model import LogisticRegression
13: from sklearn.metrics import (
14:     accuracy_score,
15:     average_precision_score,
16:     confusion_matrix,
17:     roc_auc_score,
18: )
19: from sklearn.model_selection import StratifiedKFold
20: from sklearn.pipeline import Pipeline
21: from sklearn.preprocessing import StandardScaler
22: 
23: 
24: warnings.filterwarnings("ignore")
25: 
26: ROOT = Path(__file__).resolve().parents[1]
27: 
28: DATA = (
29:     ROOT
30:     / "data"
31:     / "processed"
32:     / "gse94533_14mirna_clinical_benchmark.csv.gz"
33: )
34: 
35: OUTDIR = ROOT / "reports" / "tables"
36: 
37: MIRNAS = [
38:     "hsa-miR-23b-3p",
39:     "hsa-miR-29a-3p",
40:     "hsa-miR-32-5p",
41:     "hsa-miR-92a-3p",
42:     "hsa-miR-150-5p",
43:     "hsa-miR-200a-3p",
44:     "hsa-miR-200c-3p",
45:     "hsa-miR-203a",
46:     "hsa-miR-320c",
47:     "hsa-miR-320d",
48:     "hsa-miR-335-5p",
49:     "hsa-miR-450b-5p",
50:     "hsa-miR-1246",
51:     "hsa-miR-1307-5p",
52: ]
53: 
54: 
55: def run_model(df, features, name):
56: 
57:     X = df[features]
58:     y = df["binary_label"].astype(int)
59: 
60:     numeric_pipe = Pipeline([
61:         ("imputer", SimpleImputer(strategy="median")),
62:         ("scaler", StandardScaler()),
63:         (
64:             "model",
65:             LogisticRegression(
66:                 C=1.0,
67:                 max_iter=5000,
68:                 solver="liblinear",
69:             ),
70:         ),
71:     ])
72: 
73:     cv = StratifiedKFold(
74:         n_splits=5,
75:         shuffle=True,
76:         random_state=42,
77:     )
78: 
79:     oof_prob = np.full(len(df), np.nan)
80:     oof_pred = np.full(len(df), -1, dtype=int)
81: 
82:     for train_idx, test_idx in cv.split(X, y):
83: 
84:         model = numeric_pipe
85: 
86:         model.fit(
87:             X.iloc[train_idx],
88:             y.iloc[train_idx],
89:         )
90: 
91:         oof_prob[test_idx] = model.predict_proba(
92:             X.iloc[test_idx]
93:         )[:, 1]
94: 
95:         oof_pred[test_idx] = (
96:             oof_prob[test_idx] >= 0.5
97:         ).astype(int)
98: 
99:     auc = roc_auc_score(y, oof_prob)
100:     pr_auc = average_precision_score(y, oof_prob)
101:     accuracy = accuracy_score(y, oof_pred)
102: 
103:     tn, fp, fn, tp = confusion_matrix(
104:         y,
105:         oof_pred,
106:     ).ravel()
107: 
108:     sensitivity = tp / (tp + fn)
109:     specificity = tn / (tn + fp)
110: 
111:     pred = pd.DataFrame({
112:         "sample_id": df["sample_id"],
113:         "true_label": y,
114:         "oof_probability": oof_prob,
115:         "oof_prediction": oof_pred,
116:     })
117: 
118:     pred_path = (
119:         OUTDIR
120:         / f"{name}_oof_predictions.csv"
121:     )
122: 
123:     result_path = (
124:         OUTDIR
125:         / f"{name}_results.json"
126:     )
127: 
128:     pred.to_csv(
129:         pred_path,
130:         index=False,
131:     )
132: 
133:     results = {
134:         "model": name,
135:         "patients": int(len(df)),
136:         "features": features,
137:         "roc_auc": float(auc),
138:         "pr_auc": float(pr_auc),
139:         "accuracy": float(accuracy),
140:         "sensitivity": float(sensitivity),
141:         "specificity": float(specificity),
142:         "tp": int(tp),
143:         "fn": int(fn),
144:         "fp": int(fp),
145:         "tn": int(tn),
146:     }
147: 
148:     result_path.write_text(
149:         json.dumps(results, indent=2)
150:     )
151: 
152:     return results
153: 
154: 
155: def main():
156: 
157:     print("=" * 80)
158:     print("AGE / CLINICAL OVARIAN CANCER MODELS")
159:     print("=" * 80)
160: 
161:     df = pd.read_csv(DATA)
162: 
163:     if len(df) != 179:
164:         raise RuntimeError(
165:             f"Expected 179 patients, found {len(df)}"
166:         )
167: 
168:     required = {
169:         "sample_id",
170:         "binary_label",
171:         "age",
172:         "ca125",
173:         *MIRNAS,
174:     }
175: 
176:     missing = required - set(df.columns)
177: 
178:     if missing:
179:         raise RuntimeError(
180:             "Missing columns: "
181:             + ", ".join(sorted(missing))
182:         )
183: 
184:     print()
185:     print("Patients:", len(df))
186:     print("Cancer:", int(df["binary_label"].sum()))
187:     print(
188:         "Non-cancer:",
189:         int((df["binary_label"] == 0).sum())
190:     )
191:     print("Age missing:", int(df["age"].isna().sum()))
192:     print("CA-125 missing:", int(df["ca125"].isna().sum()))
193: 
194:     # ---------------------------------------------------------------
195:     # MODEL 1
196:     # ---------------------------------------------------------------
197: 
198:     r1 = run_model(
199:         df,
200:         ["age"] + MIRNAS,
201:         "age_plus_14mirna",
202:     )
203: 
204:     # ---------------------------------------------------------------
205:     # MODEL 2
206:     # ---------------------------------------------------------------
207: 
208:     r2 = run_model(
209:         df,
210:         ["age", "ca125"] + MIRNAS,
211:         "age_ca125_plus_14mirna",
212:     )
213: 
214:     print()
215:     print("=" * 80)
216:     print("FINAL RESULTS")
217:     print("=" * 80)
218: 
219:     for r in [r1, r2]:
220: 
221:         print()
222:         print(r["model"])
223:         print(
224:             f"  ROC AUC:    {r['roc_auc']:.6f}"
225:         )
226:         print(
227:             f"  PR AUC:     {r['pr_auc']:.6f}"
228:         )
229:         print(
230:             f"  Accuracy:   {r['accuracy']:.6f}"
231:         )
232:         print(
233:             f"  Sensitivity:{r['sensitivity']:.6f}"
234:         )
235:         print(
236:             f"  Specificity:{r['specificity']:.6f}"
237:         )
238: 
239:     print()
240:     print("=" * 80)
241:     print("FILES SAVED")
242:     print("=" * 80)
243: 
244:     print(
245:         OUTDIR
246:         / "age_plus_14mirna_oof_predictions.csv"
247:     )
248: 
249:     print(
250:         OUTDIR
251:         / "age_ca125_plus_14mirna_oof_predictions.csv"
252:     )
253: 
254: 
255: if __name__ == "__main__":
256:     main()
~~~

## scripts/build_ca125_baseline.py

~~~python
1: #!/usr/bin/env python3
2: 
3: """
4: CA-125 BASELINE
5: 
6: Evaluates CA-125 alone on the exact 179-patient historical cohort.
7: 
8: This is a NON-ML clinical baseline.
9: 
10: Positive prediction:
11:     higher CA-125 -> greater probability of cancer
12: 
13: Patients:
14:     179
15:     97 cancer
16:     82 non-cancer
17: 
18: CA-125 values are extracted from published eLife Table 2.
19: 
20: Missing CA-125 values are handled using median imputation.
21: 
22: Evaluation:
23:     5-fold patient-level OOF predictions
24: 
25: The purpose is to establish a clinically meaningful baseline
26: before adding miRNA features.
27: """
28: 
29: from __future__ import annotations
30: 
31: import json
32: import re
33: import zipfile
34: from pathlib import Path
35: from xml.etree import ElementTree as ET
36: 
37: import numpy as np
38: import pandas as pd
39: 
40: from sklearn.impute import SimpleImputer
41: from sklearn.metrics import (
42:     accuracy_score,
43:     average_precision_score,
44:     confusion_matrix,
45:     roc_auc_score,
46: )
47: from sklearn.model_selection import StratifiedKFold
48: 
49: 
50: ROOT = Path.home() / "ovarian-ai"
51: 
52: DOCX = (
53:     ROOT
54:     / "data"
55:     / "external"
56:     / "elife_28932"
57:     / "elife_poa_e28932_Supplementary_File_1.docx"
58: )
59: 
60: HISTORICAL_PREDICTIONS = (
61:     ROOT
62:     / "reports"
63:     / "tables"
64:     / "original_historical_clean_predictions.csv"
65: )
66: 
67: OUTPUT_DIR = (
68:     ROOT
69:     / "reports"
70:     / "tables"
71: )
72: 
73: PREDICTIONS_OUTPUT = (
74:     OUTPUT_DIR
75:     / "ca125_baseline_oof_predictions.csv"
76: )
77: 
78: RESULTS_OUTPUT = (
79:     OUTPUT_DIR
80:     / "ca125_baseline_results.json"
81: )
82: 
83: 
84: def norm_id(value):
85: 
86:     s = str(value).strip().lower()
87: 
88:     if re.fullmatch(r"\d+\.0+", s):
89:         s = s.split(".", 1)[0]
90: 
91:     return re.sub(
92:         r"[^a-z0-9]",
93:         "",
94:         s,
95:     )
96: 
97: 
98: def cell_text(cell):
99: 
100:     parts = []
101: 
102:     for node in cell.iter():
103: 
104:         if node.tag.endswith("}t"):
105:             parts.append(node.text or "")
106: 
107:     return "".join(parts).strip()
108: 
109: 
110: def extract_tables(path):
111: 
112:     with zipfile.ZipFile(path) as z:
113:         xml = z.read(
114:             "word/document.xml"
115:         )
116: 
117:     root = ET.fromstring(xml)
118: 
119:     tables = []
120: 
121:     for tbl in root.iter():
122: 
123:         if not tbl.tag.endswith("}tbl"):
124:             continue
125: 
126:         rows = []
127: 
128:         for tr in tbl:
129: 
130:             if not tr.tag.endswith("}tr"):
131:                 continue
132: 
133:             cells = []
134: 
135:             for tc in tr:
136: 
137:                 if tc.tag.endswith("}tc"):
138:                     cells.append(
139:                         cell_text(tc)
140:                     )
141: 
142:             if cells:
143:                 rows.append(cells)
144: 
145:         if rows:
146:             tables.append(rows)
147: 
148:     return tables
149: 
150: 
151: # =====================================================================
152: # LOAD EXACT HISTORICAL COHORT
153: # =====================================================================
154: 
155: print("=" * 90)
156: print("CA-125 BASELINE")
157: print("=" * 90)
158: 
159: pred = pd.read_csv(
160:     HISTORICAL_PREDICTIONS
161: )
162: 
163: if len(pred) != 179:
164:     raise RuntimeError(
165:         f"Expected 179 patients, found {len(pred)}."
166:     )
167: 
168: pred["sample_key"] = (
169:     pred["sample_id"].map(norm_id)
170: )
171: 
172: label_map = (
173:     pred
174:     .set_index("sample_key")["true_label"]
175:     .astype(int)
176: )
177: 
178: print()
179: print("Historical cohort:", len(pred))
180: print(
181:     "Cancer:",
182:     int((pred.true_label == 1).sum())
183: )
184: print(
185:     "Non-cancer:",
186:     int((pred.true_label == 0).sum())
187: )
188: 
189: 
190: # =====================================================================
191: # READ TABLE 2
192: # =====================================================================
193: 
194: print()
195: print("Reading published Table 2...")
196: 
197: tables = extract_tables(DOCX)
198: 
199: if len(tables) < 2:
200:     raise RuntimeError(
201:         "Could not find two supplementary tables."
202:     )
203: 
204: table2 = tables[1]
205: 
206: header = table2[0]
207: 
208: try:
209:     id_col = header.index("ID")
210: except ValueError:
211:     raise RuntimeError(
212:         "Could not locate ID column."
213:     )
214: 
215: try:
216:     ca_col = header.index("CA-125")
217: except ValueError:
218:     raise RuntimeError(
219:         "Could not locate CA-125 column."
220:     )
221: 
222: 
223: # =====================================================================
224: # PARSE PATIENT ROWS
225: # =====================================================================
226: 
227: rows = []
228: 
229: for row in table2[1:]:
230: 
231:     if len(row) <= max(
232:         id_col,
233:         ca_col,
234:     ):
235:         continue
236: 
237:     sid = row[id_col].strip()
238: 
239:     if not sid:
240:         continue
241: 
242:     if sid.lower() == "id":
243:         continue
244: 
245:     ca125 = row[ca_col].strip()
246: 
247:     rows.append(
248:         {
249:             "table2_id": sid,
250:             "ca125_raw": ca125,
251:         }
252:     )
253: 
254: 
255: ca = pd.DataFrame(rows)
256: 
257: print()
258: print(
259:     "Table 2 patient rows:",
260:     len(ca),
261: )
262: 
263: 
264: # =====================================================================
265: # IMPORTANT:
266: # TABLE 2 IDs ARE ACTUAL SAMPLE IDS
267: # =====================================================================
268: 
269: # ---------------------------------------------------------------
270: # Table 2 uses TWO kinds of identifiers:
271: #
272: #   numeric 1..120  = historical patient positions
273: #   e1001...        = actual e-sample IDs
274: #
275: # The historical PMML reproduction established that original
276: # numeric position 56 was excluded because it corresponds to 81195.
277: #
278: # Therefore:
279: #
280: #   Table 2 ID 1..55  -> prediction rows 1..55
281: #   Table 2 ID 57..120 -> prediction rows 56..119
282: #
283: # e-samples retain their explicit IDs.
284: # ---------------------------------------------------------------
285: 
286: prediction_by_key = (
287:     pred
288:     .set_index("sample_key")
289: )
290: 
291: def table2_to_actual_sample_id(value):
292: 
293:     key = norm_id(value)
294: 
295:     # Explicit e-sample.
296:     if key.startswith("e"):
297:         return key
298: 
299:     # Numeric Table-2 position.
300:     try:
301:         position = int(float(value))
302:     except (TypeError, ValueError):
303:         raise RuntimeError(
304:             f"Could not interpret Table 2 ID: {value!r}"
305:         )
306: 
307:     if position == 56:
308:         raise RuntimeError(
309:             "Table 2 unexpectedly contains excluded position 56."
310:         )
311: 
312:     if position < 1 or position > 120:
313:         raise RuntimeError(
314:             f"Unexpected numeric Table 2 ID: {position}"
315:         )
316: 
317:     # Convert original historical position to retained prediction
318:     # row after excluding original position 56.
319:     retained_index = (
320:         position - 1
321:         if position <= 55
322:         else position - 2
323:     )
324: 
325:     return pred.iloc[
326:         retained_index
327:     ]["sample_key"]
328: 
329: 
330: ca["actual_sample_key"] = (
331:     ca["table2_id"].map(
332:         table2_to_actual_sample_id
333:     )
334: )
335: 
336: ca["true_label"] = (
337:     ca["actual_sample_key"].map(
338:         label_map
339:     )
340: )
341: 
342: 
343: if ca["true_label"].isna().any():
344: 
345:     missing = ca.loc[
346:         ca["true_label"].isna(),
347:         "table2_id",
348:     ].tolist()
349: 
350:     raise RuntimeError(
351:         "Missing labels for Table 2 IDs:\n"
352:         + "\n".join(
353:             map(str, missing)
354:         )
355:     )
356: 
357: 
358: # =====================================================================
359: # PARSE CA-125
360: # =====================================================================
361: 
362: def parse_ca125(value):
363: 
364:     s = str(value).strip()
365: 
366:     if (
367:         not s
368:         or s in {"â€¡", "*", "nan"}
369:     ):
370:         return np.nan
371: 
372:     try:
373:         return float(s)
374:     except ValueError:
375:         return np.nan
376: 
377: 
378: ca["ca125"] = (
379:     ca["ca125_raw"].map(
380:         parse_ca125
381:     )
382: )
383: 
384: 
385: # =====================================================================
386: # ALIGN EXACTLY TO HISTORICAL COHORT
387: # =====================================================================
388: 
389: historical = (
390:     pred[
391:         [
392:             "sample_id",
393:             "sample_key",
394:             "true_label",
395:         ]
396:     ]
397:     .merge(
398:         ca[
399:             [
400:                 "actual_sample_key",
401:                 "ca125",
402:             ]
403:         ],
404:         left_on="sample_key",
405:         right_on="actual_sample_key",
406:         how="left",
407:         validate="one_to_one",
408:     )
409: )
410: 
411: 
412: if len(historical) != 179:
413:     raise RuntimeError(
414:         f"Expected 179 aligned patients, got {len(historical)}."
415:     )
416: 
417: 
418: # =====================================================================
419: # VALIDATE
420: # =====================================================================
421: 
422: y = historical[
423:     "true_label"
424: ].to_numpy(dtype=int)
425: 
426: x_raw = historical[
427:     ["ca125"]
428: ].to_numpy(dtype=float)
429: 
430: print()
431: print("=" * 90)
432: print("CA-125 DATA")
433: print("=" * 90)
434: 
435: print()
436: print(
437:     "Patients:",
438:     len(historical),
439: )
440: 
441: print(
442:     "Cancer:",
443:     int((y == 1).sum())
444: )
445: 
446: print(
447:     "Non-cancer:",
448:     int((y == 0).sum())
449: )
450: 
451: print(
452:     "Observed CA-125:",
453:     int(np.isfinite(x_raw).sum())
454: )
455: 
456: print(
457:     "Missing CA-125:",
458:     int(np.isnan(x_raw).sum())
459: )
460: 
461: 
462: # =====================================================================
463: # OOF MEDIAN IMPUTATION
464: # =====================================================================
465: 
466: cv = StratifiedKFold(
467:     n_splits=5,
468:     shuffle=True,
469:     random_state=2026,
470: )
471: 
472: oof = np.full(
473:     len(y),
474:     np.nan,
475: )
476: 
477: 
478: for fold, (
479:     train_idx,
480:     test_idx,
481: ) in enumerate(
482:     cv.split(
483:         x_raw,
484:         y,
485:     ),
486:     start=1,
487: ):
488: 
489:     imputer = SimpleImputer(
490:         strategy="median"
491:     )
492: 
493:     x_train = imputer.fit_transform(
494:         x_raw[train_idx]
495:     )
496: 
497:     x_test = imputer.transform(
498:         x_raw[test_idx]
499:     )
500: 
501:     # CA-125 itself is the score.
502:     # Higher CA-125 => greater cancer likelihood.
503:     oof[test_idx] = (
504:         x_test[:, 0]
505:     )
506: 
507: 
508: # =====================================================================
509: # METRICS
510: # =====================================================================
511: 
512: auc = roc_auc_score(
513:     y,
514:     oof,
515: )
516: 
517: pr_auc = average_precision_score(
518:     y,
519:     oof,
520: )
521: 
522: threshold = 35.0
523: 
524: prediction = (
525:     oof > threshold
526: ).astype(int)
527: 
528: accuracy = accuracy_score(
529:     y,
530:     prediction,
531: )
532: 
533: tn, fp, fn, tp = confusion_matrix(
534:     y,
535:     prediction,
536: ).ravel()
537: 
538: sensitivity = (
539:     tp / (tp + fn)
540: )
541: 
542: specificity = (
543:     tn / (tn + fp)
544: )
545: 
546: 
547: # =====================================================================
548: # SAVE
549: # =====================================================================
550: 
551: historical[
552:     "oof_ca125"
553: ] = oof
554: 
555: historical[
556:     "prediction_ca125_gt35"
557: ] = prediction
558: 
559: OUTPUT_DIR.mkdir(
560:     parents=True,
561:     exist_ok=True,
562: )
563: 
564: historical.to_csv(
565:     PREDICTIONS_OUTPUT,
566:     index=False,
567: )
568: 
569: 
570: results = {
571:     "model": "CA-125 only",
572:     "cohort": {
573:         "patients": 179,
574:         "cancer": int((y == 1).sum()),
575:         "non_cancer": int((y == 0).sum()),
576:     },
577:     "missing_ca125": int(
578:         np.isnan(x_raw).sum()
579:     ),
580:     "evaluation": {
581:         "method": "5-fold patient-level OOF",
582:         "random_state": 2026,
583:     },
584:     "threshold": 35.0,
585:     "metrics": {
586:         "roc_auc": float(auc),
587:         "pr_auc": float(pr_auc),
588:         "accuracy": float(accuracy),
589:         "sensitivity": float(sensitivity),
590:         "specificity": float(specificity),
591:     },
592:     "confusion_matrix": {
593:         "TN": int(tn),
594:         "FP": int(fp),
595:         "FN": int(fn),
596:         "TP": int(tp),
597:     },
598: }
599: 
600: 
601: with RESULTS_OUTPUT.open(
602:     "w"
603: ) as f:
604: 
605:     json.dump(
606:         results,
607:         f,
608:         indent=2,
609:     )
610: 
611: 
612: # =====================================================================
613: # REPORT
614: # =====================================================================
615: 
616: print()
617: print("=" * 90)
618: print("CA-125 BASELINE RESULTS")
619: print("=" * 90)
620: 
621: print()
622: print(f"ROC AUC:      {auc:.6f}")
623: print(f"PR AUC:       {pr_auc:.6f}")
624: print(f"Accuracy:     {accuracy:.6f}")
625: print(f"Sensitivity:  {sensitivity:.6f}")
626: print(f"Specificity:  {specificity:.6f}")
627: 
628: print()
629: print("Confusion matrix:")
630: print(f"  TP={tp} FN={fn}")
631: print(f"  FP={fp} TN={tn}")
632: 
633: print()
634: print("Saved:")
635: print(PREDICTIONS_OUTPUT)
636: print(RESULTS_OUTPUT)
~~~

## scripts/build_ca125_plus_14mirna.py

~~~python
1: #!/usr/bin/env python3
2: 
3: from __future__ import annotations
4: 
5: import json
6: import re
7: import zipfile
8: from pathlib import Path
9: from xml.etree import ElementTree as ET
10: 
11: import numpy as np
12: import pandas as pd
13: 
14: from sklearn.impute import SimpleImputer
15: from sklearn.linear_model import LogisticRegression
16: from sklearn.metrics import (
17:     accuracy_score,
18:     average_precision_score,
19:     confusion_matrix,
20:     roc_auc_score,
21: )
22: from sklearn.model_selection import (
23:     GridSearchCV,
24:     StratifiedKFold,
25: )
26: from sklearn.pipeline import Pipeline
27: from sklearn.preprocessing import StandardScaler
28: 
29: 
30: ROOT = Path.home() / "ovarian-ai"
31: 
32: PREDICTIONS = (
33:     ROOT
34:     / "reports"
35:     / "tables"
36:     / "original_historical_clean_predictions.csv"
37: )
38: 
39: DATASET = (
40:     ROOT
41:     / "data"
42:     / "processed"
43:     / "gse94533_14mirna_benchmark.csv.gz"
44: )
45: 
46: CA125_PREDICTIONS = (
47:     ROOT
48:     / "reports"
49:     / "tables"
50:     / "ca125_baseline_oof_predictions.csv"
51: )
52: 
53: OUTPUT_DIR = (
54:     ROOT
55:     / "reports"
56:     / "tables"
57: )
58: 
59: OUT = (
60:     OUTPUT_DIR
61:     / "ca125_plus_14mirna_oof_predictions.csv"
62: )
63: 
64: RESULTS = (
65:     OUTPUT_DIR
66:     / "ca125_plus_14mirna_results.json"
67: )
68: 
69: FEATURES = [
70:     "hsa-miR-23b-3p",
71:     "hsa-miR-29a-3p",
72:     "hsa-miR-32-5p",
73:     "hsa-miR-92a-3p",
74:     "hsa-miR-150-5p",
75:     "hsa-miR-200a-3p",
76:     "hsa-miR-200c-3p",
77:     "hsa-miR-203a",
78:     "hsa-miR-320c",
79:     "hsa-miR-320d",
80:     "hsa-miR-335-5p",
81:     "hsa-miR-450b-5p",
82:     "hsa-miR-1246",
83:     "hsa-miR-1307-5p",
84: ]
85: 
86: 
87: def norm_id(x):
88:     s = str(x).strip().lower()
89:     if re.fullmatch(r"\d+\.0+", s):
90:         s = s.split(".", 1)[0]
91:     return re.sub(r"[^a-z0-9]", "", s)
92: 
93: 
94: print("=" * 90)
95: print("CA-125 + 14-miRNA MODERN BASELINE")
96: print("=" * 90)
97: 
98: pred = pd.read_csv(PREDICTIONS)
99: 
100: if len(pred) != 179:
101:     raise RuntimeError(f"Expected 179 patients, found {len(pred)}.")
102: 
103: pred["sample_key"] = pred["sample_id"].map(norm_id)
104: 
105: # ---------------------------------------------------------------
106: # Load exact 14-miRNA historical representation.
107: # ---------------------------------------------------------------
108: 
109: data = pd.read_csv(DATASET)
110: 
111: if len(data) != 179:
112:     raise RuntimeError(f"Expected 179 dataset rows, found {len(data)}.")
113: 
114: missing_features = [
115:     f for f in FEATURES
116:     if f not in data.columns
117: ]
118: 
119: if missing_features:
120:     raise RuntimeError(
121:         f"Missing 14-miRNA features: {missing_features}"
122:     )
123: 
124: data["sample_key"] = data["sample_id"].map(norm_id)
125: 
126: data = data[
127:     ["sample_key"] + FEATURES
128: ]
129: 
130: # ---------------------------------------------------------------
131: # Load CA-125.
132: # ---------------------------------------------------------------
133: 
134: ca = pd.read_csv(
135:     CA125_PREDICTIONS
136: )
137: 
138: if len(ca) != 179:
139:     raise RuntimeError(
140:         f"Expected 179 CA-125 rows, found {len(ca)}."
141:     )
142: 
143: ca["sample_key"] = ca["sample_id"].map(norm_id)
144: 
145: if "ca125" not in ca.columns:
146:     raise RuntimeError(
147:         "CA-125 prediction file does not contain 'ca125'."
148:     )
149: 
150: # ---------------------------------------------------------------
151: # Merge everything.
152: # ---------------------------------------------------------------
153: 
154: df = (
155:     pred[
156:         [
157:             "sample_id",
158:             "sample_key",
159:             "true_label",
160:         ]
161:     ]
162:     .merge(
163:         data,
164:         on="sample_key",
165:         how="left",
166:         validate="one_to_one",
167:     )
168:     .merge(
169:         ca[
170:             [
171:                 "sample_key",
172:                 "ca125",
173:             ]
174:         ],
175:         on="sample_key",
176:         how="left",
177:         validate="one_to_one",
178:     )
179: )
180: 
181: if len(df) != 179:
182:     raise RuntimeError(
183:         f"Expected 179 merged patients, got {len(df)}."
184:     )
185: 
186: # ---------------------------------------------------------------
187: # Build matrix.
188: # ---------------------------------------------------------------
189: 
190: X = df[
191:     ["ca125"] + FEATURES
192: ].to_numpy(
193:     dtype=float
194: )
195: 
196: y = df[
197:     "true_label"
198: ].to_numpy(
199:     dtype=int
200: )
201: 
202: if (y == 1).sum() != 97:
203:     raise RuntimeError("Expected 97 cancer patients.")
204: 
205: if (y == 0).sum() != 82:
206:     raise RuntimeError("Expected 82 non-cancer patients.")
207: 
208: print()
209: print("Patients:", len(df))
210: print("Features:", X.shape[1])
211: print("CA-125 + 14 miRNAs")
212: print("Cancer:", int((y == 1).sum()))
213: print("Non-cancer:", int((y == 0).sum()))
214: 
215: # ---------------------------------------------------------------
216: # Nested patient-level CV.
217: # ---------------------------------------------------------------
218: 
219: outer = StratifiedKFold(
220:     n_splits=5,
221:     shuffle=True,
222:     random_state=2026,
223: )
224: 
225: oof_probability = np.full(
226:     len(y),
227:     np.nan,
228: )
229: 
230: oof_prediction = np.full(
231:     len(y),
232:     -1,
233:     dtype=int,
234: )
235: 
236: fold_records = []
237: 
238: for fold, (train_idx, test_idx) in enumerate(
239:     outer.split(X, y),
240:     start=1,
241: ):
242: 
243:     print()
244:     print(f"FOLD {fold}/5")
245: 
246:     inner = StratifiedKFold(
247:         n_splits=3,
248:         shuffle=True,
249:         random_state=2026 + fold,
250:     )
251: 
252:     pipe = Pipeline(
253:         steps=[
254:             (
255:                 "imputer",
256:                 SimpleImputer(strategy="median"),
257:             ),
258:             (
259:                 "scaler",
260:                 StandardScaler(),
261:             ),
262:             (
263:                 "classifier",
264:                 LogisticRegression(
265:                     solver="saga",
266:                     max_iter=10000,
267:                     random_state=2026,
268:                 ),
269:             ),
270:         ]
271:     )
272: 
273:     grid = {
274:         "classifier__C": [
275:             0.01,
276:             0.03,
277:             0.1,
278:             0.3,
279:             1.0,
280:             3.0,
281:         ],
282:         "classifier__l1_ratio": [
283:             0.0,
284:             0.25,
285:             0.5,
286:             0.75,
287:             1.0,
288:         ],
289:         "classifier__penalty": [
290:             "elasticnet",
291:         ],
292:     }
293: 
294:     search = GridSearchCV(
295:         pipe,
296:         grid,
297:         scoring="roc_auc",
298:         cv=inner,
299:         n_jobs=-1,
300:         refit=True,
301:     )
302: 
303:     search.fit(
304:         X[train_idx],
305:         y[train_idx],
306:     )
307: 
308:     probability = search.predict_proba(
309:         X[test_idx]
310:     )[:, 1]
311: 
312:     prediction = (
313:         probability >= 0.5
314:     ).astype(int)
315: 
316:     oof_probability[test_idx] = probability
317:     oof_prediction[test_idx] = prediction
318: 
319:     fold_auc = roc_auc_score(
320:         y[test_idx],
321:         probability,
322:     )
323: 
324:     fold_accuracy = accuracy_score(
325:         y[test_idx],
326:         prediction,
327:     )
328: 
329:     fold_records.append(
330:         {
331:             "fold": fold,
332:             "test_auc": float(fold_auc),
333:             "test_accuracy": float(fold_accuracy),
334:             "inner_cv_auc": float(search.best_score_),
335:             "best_params": search.best_params_,
336:         }
337:     )
338: 
339:     print(
340:         f"  Test AUC: {fold_auc:.6f}"
341:     )
342: 
343:     print(
344:         f"  Accuracy: {fold_accuracy:.6f}"
345:     )
346: 
347:     print(
348:         f"  Best params: {search.best_params_}"
349:     )
350: 
351: # ---------------------------------------------------------------
352: # Final OOF metrics.
353: # ---------------------------------------------------------------
354: 
355: auc = roc_auc_score(
356:     y,
357:     oof_probability,
358: )
359: 
360: pr_auc = average_precision_score(
361:     y,
362:     oof_probability,
363: )
364: 
365: accuracy = accuracy_score(
366:     y,
367:     oof_prediction,
368: )
369: 
370: tn, fp, fn, tp = confusion_matrix(
371:     y,
372:     oof_prediction,
373: ).ravel()
374: 
375: sensitivity = tp / (tp + fn)
376: specificity = tn / (tn + fp)
377: 
378: print()
379: print("=" * 90)
380: print("CA-125 + 14-miRNA RESULTS")
381: print("=" * 90)
382: 
383: print()
384: print(f"OOF ROC AUC:   {auc:.6f}")
385: print(f"OOF PR AUC:    {pr_auc:.6f}")
386: print(f"OOF Accuracy:  {accuracy:.6f}")
387: print(f"Sensitivity:   {sensitivity:.6f}")
388: print(f"Specificity:   {specificity:.6f}")
389: 
390: print()
391: print("Confusion matrix:")
392: print(f"  TP={tp} FN={fn}")
393: print(f"  FP={fp} TN={tn}")
394: 
395: # ---------------------------------------------------------------
396: # Save.
397: # ---------------------------------------------------------------
398: 
399: df["oof_probability"] = oof_probability
400: df["oof_prediction"] = oof_prediction
401: 
402: df.to_csv(
403:     OUT,
404:     index=False,
405: )
406: 
407: results = {
408:     "model": "CA-125 + 14-miRNA elastic-net logistic regression",
409:     "patients": 179,
410:     "cancer": 97,
411:     "non_cancer": 82,
412:     "n_features": 15,
413:     "cross_validation": "5-fold nested patient-level CV",
414:     "metrics": {
415:         "roc_auc": float(auc),
416:         "pr_auc": float(pr_auc),
417:         "accuracy": float(accuracy),
418:         "sensitivity": float(sensitivity),
419:         "specificity": float(specificity),
420:     },
421:     "confusion_matrix": {
422:         "TN": int(tn),
423:         "FP": int(fp),
424:         "FN": int(fn),
425:         "TP": int(tp),
426:     },
427:     "folds": fold_records,
428: }
429: 
430: with RESULTS.open("w") as f:
431:     json.dump(
432:         results,
433:         f,
434:         indent=2,
435:     )
436: 
437: print()
438: print("Saved:")
439: print(OUT)
440: print(RESULTS)
~~~

## scripts/build_cancer_vs_borderline_age_model.py

~~~python
1: #!/usr/bin/env python3
2: 
3: from pathlib import Path
4: import json
5: 
6: import numpy as np
7: import pandas as pd
8: 
9: from sklearn.feature_selection import SelectKBest, f_classif
10: from sklearn.linear_model import LogisticRegression
11: from sklearn.metrics import (
12:     accuracy_score,
13:     average_precision_score,
14:     confusion_matrix,
15:     roc_auc_score,
16: )
17: from sklearn.model_selection import GridSearchCV, StratifiedKFold
18: from sklearn.pipeline import Pipeline
19: from sklearn.preprocessing import StandardScaler
20: 
21: 
22: ROOT = Path.home() / "ovarian-ai"
23: 
24: DATA = (
25:     ROOT
26:     / "data"
27:     / "processed"
28:     / "gse106817_external_validation.csv.gz"
29: )
30: 
31: OUT_PRED = (
32:     ROOT
33:     / "reports"
34:     / "tables"
35:     / "gse106817_cancer_vs_borderline_age_predictions.csv"
36: )
37: 
38: OUT_RESULTS = (
39:     ROOT
40:     / "reports"
41:     / "tables"
42:     / "gse106817_cancer_vs_borderline_age_results.json"
43: )
44: 
45: 
46: print("=" * 80)
47: print("OVARIAN CANCER vs BORDERLINE â€” 170 miRNAs + AGE")
48: print("=" * 80)
49: 
50: df = pd.read_csv(DATA)
51: 
52: df = df[
53:     df["group"].isin([
54:         "cancer",
55:         "borderline",
56:     ])
57: ].copy()
58: 
59: mirna_cols = [
60:     c for c in df.columns
61:     if c.startswith("hsa-")
62: ]
63: 
64: feature_cols = mirna_cols + ["age"]
65: 
66: X = df[feature_cols].apply(
67:     pd.to_numeric,
68:     errors="coerce",
69: )
70: 
71: y = (
72:     df["group"] == "cancer"
73: ).astype(int)
74: 
75: valid = X.notna().all(axis=1)
76: 
77: X = X.loc[valid].reset_index(drop=True)
78: y = y.loc[valid].reset_index(drop=True)
79: meta = df.loc[valid].reset_index(drop=True)
80: 
81: print()
82: print("Patients:", len(X))
83: print("Cancer:", int(y.sum()))
84: print("Borderline:", int((y == 0).sum()))
85: print("Features:", len(feature_cols))
86: print("Age missing:", int(X["age"].isna().sum()))
87: 
88: outer_cv = StratifiedKFold(
89:     n_splits=5,
90:     shuffle=True,
91:     random_state=42,
92: )
93: 
94: inner_cv = StratifiedKFold(
95:     n_splits=4,
96:     shuffle=True,
97:     random_state=43,
98: )
99: 
100: pipeline = Pipeline([
101:     (
102:         "scaler",
103:         StandardScaler(),
104:     ),
105:     (
106:         "select",
107:         SelectKBest(
108:             score_func=f_classif,
109:         ),
110:     ),
111:     (
112:         "classifier",
113:         LogisticRegression(
114:             max_iter=5000,
115:             class_weight="balanced",
116:             solver="liblinear",
117:         ),
118:     ),
119: ])
120: 
121: param_grid = {
122:     "select__k": [
123:         5,
124:         10,
125:         20,
126:         30,
127:         50,
128:         75,
129:         100,
130:         "all",
131:     ],
132:     "classifier__C": [
133:         0.03,
134:         0.1,
135:         0.3,
136:         1.0,
137:         3.0,
138:     ],
139: }
140: 
141: oof_prob = np.full(
142:     len(X),
143:     np.nan,
144: )
145: 
146: oof_pred = np.full(
147:     len(X),
148:     np.nan,
149: )
150: 
151: fold_results = []
152: 
153: print()
154: print("=" * 80)
155: print("NESTED CROSS-VALIDATION")
156: print("=" * 80)
157: 
158: for fold, (train_idx, test_idx) in enumerate(
159:     outer_cv.split(X, y),
160:     start=1,
161: ):
162: 
163:     search = GridSearchCV(
164:         pipeline,
165:         param_grid,
166:         cv=inner_cv,
167:         scoring="roc_auc",
168:         n_jobs=-1,
169:         refit=True,
170:     )
171: 
172:     search.fit(
173:         X.iloc[train_idx],
174:         y.iloc[train_idx],
175:     )
176: 
177:     probability = search.predict_proba(
178:         X.iloc[test_idx]
179:     )[:, 1]
180: 
181:     prediction = (
182:         probability >= 0.5
183:     ).astype(int)
184: 
185:     oof_prob[test_idx] = probability
186:     oof_pred[test_idx] = prediction
187: 
188:     fold_auc = roc_auc_score(
189:         y.iloc[test_idx],
190:         probability,
191:     )
192: 
193:     fold_acc = accuracy_score(
194:         y.iloc[test_idx],
195:         prediction,
196:     )
197: 
198:     fold_results.append({
199:         "fold": fold,
200:         "auc": float(fold_auc),
201:         "accuracy": float(fold_acc),
202:         "best_params": search.best_params_,
203:     })
204: 
205:     print(
206:         f"Fold {fold}: "
207:         f"test={len(test_idx)} "
208:         f"AUC={fold_auc:.6f} "
209:         f"accuracy={fold_acc:.6f} "
210:         f"best={search.best_params_}"
211:     )
212: 
213: auc = roc_auc_score(
214:     y,
215:     oof_prob,
216: )
217: 
218: pr_auc = average_precision_score(
219:     y,
220:     oof_prob,
221: )
222: 
223: accuracy = accuracy_score(
224:     y,
225:     oof_pred,
226: )
227: 
228: tn, fp, fn, tp = confusion_matrix(
229:     y,
230:     oof_pred,
231: ).ravel()
232: 
233: sensitivity = tp / (tp + fn)
234: specificity = tn / (tn + fp)
235: 
236: print()
237: print("=" * 80)
238: print("FINAL RESULTS")
239: print("=" * 80)
240: 
241: print(f"ROC AUC:     {auc:.6f}")
242: print(f"PR AUC:      {pr_auc:.6f}")
243: print(f"Accuracy:    {accuracy:.6f}")
244: print(f"Sensitivity: {sensitivity:.6f}")
245: print(f"Specificity: {specificity:.6f}")
246: 
247: print()
248: print("Confusion matrix:")
249: print(f"  TN={tn} FP={fp}")
250: print(f"  FN={fn} TP={tp}")
251: 
252: pred = meta[
253:     [
254:         "sample_id",
255:         "title",
256:         "group",
257:         "age",
258:     ]
259: ].copy()
260: 
261: pred["true_label"] = y.values
262: pred["cancer_probability"] = oof_prob
263: pred["predicted_label"] = oof_pred.astype(int)
264: 
265: pred.to_csv(
266:     OUT_PRED,
267:     index=False,
268: )
269: 
270: results = {
271:     "dataset": "GSE106817",
272:     "endpoint": "cancer_vs_borderline",
273:     "features": len(feature_cols),
274:     "mirnas": len(mirna_cols),
275:     "age_used": True,
276:     "patients": len(X),
277:     "cancer": int(y.sum()),
278:     "borderline": int((y == 0).sum()),
279:     "roc_auc": float(auc),
280:     "pr_auc": float(pr_auc),
281:     "accuracy": float(accuracy),
282:     "sensitivity": float(sensitivity),
283:     "specificity": float(specificity),
284:     "tn": int(tn),
285:     "fp": int(fp),
286:     "fn": int(fn),
287:     "tp": int(tp),
288:     "outer_cv": "5-fold stratified",
289:     "inner_cv": "4-fold stratified",
290:     "feature_selection_inside_training_folds": True,
291:     "fold_results": fold_results,
292: }
293: 
294: OUT_RESULTS.write_text(
295:     json.dumps(
296:         results,
297:         indent=2,
298:     )
299: )
300: 
301: print()
302: print("Saved:")
303: print(OUT_PRED)
304: print(OUT_RESULTS)
~~~

## scripts/build_cancer_vs_borderline_model.py

~~~python
1: #!/usr/bin/env python3
2: 
3: from pathlib import Path
4: import json
5: 
6: import numpy as np
7: import pandas as pd
8: 
9: from sklearn.feature_selection import SelectKBest, f_classif
10: from sklearn.linear_model import LogisticRegression
11: from sklearn.metrics import (
12:     accuracy_score,
13:     average_precision_score,
14:     confusion_matrix,
15:     roc_auc_score,
16: )
17: from sklearn.model_selection import (
18:     GridSearchCV,
19:     StratifiedKFold,
20: )
21: from sklearn.pipeline import Pipeline
22: from sklearn.preprocessing import StandardScaler
23: 
24: 
25: ROOT = Path.home() / "ovarian-ai"
26: 
27: DATA = (
28:     ROOT
29:     / "data"
30:     / "processed"
31:     / "gse106817_external_validation.csv.gz"
32: )
33: 
34: OUT_PRED = (
35:     ROOT
36:     / "reports"
37:     / "tables"
38:     / "gse106817_cancer_vs_borderline_predictions.csv"
39: )
40: 
41: OUT_RESULTS = (
42:     ROOT
43:     / "reports"
44:     / "tables"
45:     / "gse106817_cancer_vs_borderline_results.json"
46: )
47: 
48: 
49: print("=" * 80)
50: print("OVARIAN CANCER vs BORDERLINE OVARIAN TUMORS")
51: print("=" * 80)
52: 
53: df = pd.read_csv(DATA)
54: 
55: # ---------------------------------------------------------------
56: # Keep ONLY cancer and borderline.
57: # ---------------------------------------------------------------
58: 
59: df = df[
60:     df["group"].isin([
61:         "cancer",
62:         "borderline",
63:     ])
64: ].copy()
65: 
66: feature_cols = [
67:     c for c in df.columns
68:     if c.startswith("hsa-")
69: ]
70: 
71: X = df[feature_cols].apply(
72:     pd.to_numeric,
73:     errors="coerce",
74: )
75: 
76: y = (
77:     df["group"] == "cancer"
78: ).astype(int)
79: 
80: # Remove incomplete rows, if any.
81: valid = X.notna().all(axis=1)
82: 
83: X = X.loc[valid].reset_index(drop=True)
84: y = y.loc[valid].reset_index(drop=True)
85: 
86: meta = df.loc[valid].reset_index(drop=True)
87: 
88: print()
89: print("Patients:", len(X))
90: print(
91:     "Cancer:",
92:     int((y == 1).sum())
93: )
94: print(
95:     "Borderline:",
96:     int((y == 0).sum())
97: )
98: print(
99:     "Features:",
100:     len(feature_cols)
101: )
102: 
103: # ---------------------------------------------------------------
104: # Nested CV
105: #
106: # Feature selection happens INSIDE each training fold.
107: # This prevents information leakage.
108: # ---------------------------------------------------------------
109: 
110: outer_cv = StratifiedKFold(
111:     n_splits=5,
112:     shuffle=True,
113:     random_state=42,
114: )
115: 
116: inner_cv = StratifiedKFold(
117:     n_splits=4,
118:     shuffle=True,
119:     random_state=43,
120: )
121: 
122: pipeline = Pipeline([
123:     (
124:         "scaler",
125:         StandardScaler(),
126:     ),
127:     (
128:         "select",
129:         SelectKBest(
130:             score_func=f_classif,
131:         ),
132:     ),
133:     (
134:         "classifier",
135:         LogisticRegression(
136:             max_iter=5000,
137:             class_weight="balanced",
138:             solver="liblinear",
139:         ),
140:     ),
141: ])
142: 
143: param_grid = {
144:     "select__k": [
145:         5,
146:         10,
147:         20,
148:         30,
149:         50,
150:         "all",
151:     ],
152:     "classifier__C": [
153:         0.03,
154:         0.1,
155:         0.3,
156:         1.0,
157:         3.0,
158:     ],
159: }
160: 
161: oof_prob = np.full(
162:     len(X),
163:     np.nan,
164: )
165: 
166: oof_pred = np.full(
167:     len(X),
168:     np.nan,
169: )
170: 
171: fold_results = []
172: 
173: print()
174: print("=" * 80)
175: print("NESTED 5-FOLD CROSS-VALIDATION")
176: print("=" * 80)
177: 
178: for fold, (train_idx, test_idx) in enumerate(
179:     outer_cv.split(X, y),
180:     start=1,
181: ):
182: 
183:     X_train = X.iloc[train_idx]
184:     X_test = X.iloc[test_idx]
185: 
186:     y_train = y.iloc[train_idx]
187:     y_test = y.iloc[test_idx]
188: 
189:     search = GridSearchCV(
190:         estimator=pipeline,
191:         param_grid=param_grid,
192:         cv=inner_cv,
193:         scoring="roc_auc",
194:         n_jobs=-1,
195:         refit=True,
196:     )
197: 
198:     search.fit(
199:         X_train,
200:         y_train,
201:     )
202: 
203:     probability = search.predict_proba(
204:         X_test
205:     )[:, 1]
206: 
207:     prediction = (
208:         probability >= 0.5
209:     ).astype(int)
210: 
211:     oof_prob[test_idx] = probability
212:     oof_pred[test_idx] = prediction
213: 
214:     fold_auc = roc_auc_score(
215:         y_test,
216:         probability,
217:     )
218: 
219:     fold_acc = accuracy_score(
220:         y_test,
221:         prediction,
222:     )
223: 
224:     fold_results.append({
225:         "fold": fold,
226:         "test_n": int(len(test_idx)),
227:         "auc": float(fold_auc),
228:         "accuracy": float(fold_acc),
229:         "best_params": search.best_params_,
230:     })
231: 
232:     print(
233:         f"Fold {fold}: "
234:         f"test={len(test_idx)} "
235:         f"AUC={fold_auc:.6f} "
236:         f"accuracy={fold_acc:.6f}"
237:     )
238: 
239: # ---------------------------------------------------------------
240: # Overall metrics
241: # ---------------------------------------------------------------
242: 
243: auc = roc_auc_score(
244:     y,
245:     oof_prob,
246: )
247: 
248: pr_auc = average_precision_score(
249:     y,
250:     oof_prob,
251: )
252: 
253: accuracy = accuracy_score(
254:     y,
255:     oof_pred,
256: )
257: 
258: tn, fp, fn, tp = confusion_matrix(
259:     y,
260:     oof_pred,
261: ).ravel()
262: 
263: sensitivity = tp / (tp + fn)
264: specificity = tn / (tn + fp)
265: 
266: print()
267: print("=" * 80)
268: print("FINAL CANCER vs BORDERLINE RESULTS")
269: print("=" * 80)
270: 
271: print(
272:     f"ROC AUC:     {auc:.6f}"
273: )
274: 
275: print(
276:     f"PR AUC:      {pr_auc:.6f}"
277: )
278: 
279: print(
280:     f"Accuracy:    {accuracy:.6f}"
281: )
282: 
283: print(
284:     f"Sensitivity: {sensitivity:.6f}"
285: )
286: 
287: print(
288:     f"Specificity: {specificity:.6f}"
289: )
290: 
291: print()
292: print("Confusion matrix:")
293: print(
294:     f"  TN={tn} FP={fp}"
295: )
296: print(
297:     f"  FN={fn} TP={tp}"
298: )
299: 
300: # ---------------------------------------------------------------
301: # Save predictions
302: # ---------------------------------------------------------------
303: 
304: predictions = meta[
305:     [
306:         "sample_id",
307:         "title",
308:         "group",
309:     ]
310: ].copy()
311: 
312: predictions["true_label"] = y.values
313: predictions["cancer_probability"] = oof_prob
314: predictions["predicted_label"] = oof_pred.astype(int)
315: predictions["correct"] = (
316:     predictions["true_label"]
317:     == predictions["predicted_label"]
318: )
319: 
320: OUT_PRED.parent.mkdir(
321:     parents=True,
322:     exist_ok=True,
323: )
324: 
325: predictions.to_csv(
326:     OUT_PRED,
327:     index=False,
328: )
329: 
330: results = {
331:     "dataset": "GSE106817",
332:     "endpoint": "cancer_vs_borderline",
333:     "patients": int(len(X)),
334:     "cancer": int((y == 1).sum()),
335:     "borderline": int((y == 0).sum()),
336:     "features": int(len(feature_cols)),
337:     "roc_auc": float(auc),
338:     "pr_auc": float(pr_auc),
339:     "accuracy": float(accuracy),
340:     "sensitivity": float(sensitivity),
341:     "specificity": float(specificity),
342:     "tn": int(tn),
343:     "fp": int(fp),
344:     "fn": int(fn),
345:     "tp": int(tp),
346:     "outer_cv": "5-fold stratified",
347:     "inner_cv": "4-fold stratified",
348:     "feature_selection_inside_training_folds": True,
349:     "fold_results": fold_results,
350: }
351: 
352: OUT_RESULTS.write_text(
353:     json.dumps(
354:         results,
355:         indent=2,
356:     )
357: )
358: 
359: print()
360: print("Saved:")
361: print(OUT_PRED)
362: print(OUT_RESULTS)
363: 
364: print()
365: print("IMPORTANT:")
366: print(
367:     "This is an exploratory cancer-vs-borderline model."
368: )
369: print(
370:     "It is NOT yet independent external validation."
371: )
372: print(
373:     "A high AUC here does not by itself prove clinical usefulness."
374: )
~~~

## scripts/build_fair_historical_14mirna_baseline.py

~~~python
1: #!/usr/bin/env python3
2: 
3: """
4: FAIR MODERN 14-miRNA BASELINE
5: 
6: Uses EXACTLY the same historical cohort and Dataset 2 representation
7: used by the reproduced eLife PMML model.
8: 
9: Historical cohort:
10:     179 patients
11:     97 cancer
12:     82 non-cancer
13: 
14: Representation:
15:     eLife Supplementary Dataset 2
16:     batch-adjusted
17:     log10-transformed
18: 
19: Features:
20:     exact 14 PMML miRNAs
21: 
22: Model:
23:     standardized L2 logistic regression
24: 
25: Evaluation:
26:     5-fold stratified patient-level OOF predictions
27: 
28: This is a NEW model.
29: It does NOT use the historical PMML weights.
30: """
31: 
32: from __future__ import annotations
33: 
34: import json
35: from pathlib import Path
36: 
37: import numpy as np
38: import pandas as pd
39: 
40: from sklearn.linear_model import LogisticRegression
41: from sklearn.metrics import (
42:     accuracy_score,
43:     average_precision_score,
44:     confusion_matrix,
45:     roc_auc_score,
46: )
47: from sklearn.model_selection import (
48:     GridSearchCV,
49:     StratifiedKFold,
50: )
51: from sklearn.pipeline import Pipeline
52: from sklearn.preprocessing import StandardScaler
53: 
54: 
55: ROOT = Path.home() / "ovarian-ai"
56: 
57: SUPPLEMENTARY_XLSX = (
58:     ROOT
59:     / "data"
60:     / "external"
61:     / "elife_28932"
62:     / "elife_poa_e28932_Supplementary_File_6.xlsx"
63: )
64: 
65: PREDICTIONS_PATH = (
66:     ROOT
67:     / "reports"
68:     / "tables"
69:     / "original_historical_clean_predictions.csv"
70: )
71: 
72: OUTPUT_DIR = (
73:     ROOT
74:     / "reports"
75:     / "tables"
76: )
77: 
78: PREDICTION_OUTPUT = (
79:     OUTPUT_DIR
80:     / "fair_historical_14mirna_oof_predictions.csv"
81: )
82: 
83: RESULTS_OUTPUT = (
84:     OUTPUT_DIR
85:     / "fair_historical_14mirna_baseline_results.json"
86: )
87: 
88: 
89: FEATURES = [
90:     "hsa-miR-23b-3p",
91:     "hsa-miR-29a-3p",
92:     "hsa-miR-32-5p",
93:     "hsa-miR-92a-3p",
94:     "hsa-miR-150-5p",
95:     "hsa-miR-200a-3p",
96:     "hsa-miR-200c-3p",
97:     "hsa-miR-203a",
98:     "hsa-miR-320c",
99:     "hsa-miR-320d",
100:     "hsa-miR-335-5p",
101:     "hsa-miR-450b-5p",
102:     "hsa-miR-1246",
103:     "hsa-miR-1307-5p",
104: ]
105: 
106: 
107: def norm_id(value):
108: 
109:     s = str(value).strip().lower()
110: 
111:     if s.endswith(".0") and s[:-2].isdigit():
112:         s = s[:-2]
113: 
114:     return s.replace("-", "")
115: 
116: 
117: # =====================================================================
118: # LOAD HISTORICAL SAMPLE ORDER
119: # =====================================================================
120: 
121: print("=" * 90)
122: print("FAIR HISTORICAL 14-miRNA MODERN BASELINE")
123: print("=" * 90)
124: 
125: print()
126: print("Loading historical PMML prediction file...")
127: 
128: pred = pd.read_csv(PREDICTIONS_PATH)
129: 
130: if len(pred) != 179:
131:     raise RuntimeError(
132:         f"Expected 179 historical patients, found {len(pred)}."
133:     )
134: 
135: pred["sample_key"] = pred["sample_id"].map(norm_id)
136: 
137: print(
138:     f"Historical cohort: {len(pred)} patients"
139: )
140: 
141: print(
142:     f"Cancer: {(pred['true_label'] == 1).sum()}"
143: )
144: 
145: print(
146:     f"Non-cancer: {(pred['true_label'] == 0).sum()}"
147: )
148: 
149: 
150: # =====================================================================
151: # READ DATASET 2
152: # =====================================================================
153: 
154: print()
155: print("=" * 90)
156: print("READING eLIFE SUPPLEMENTAL DATASET 2")
157: print("=" * 90)
158: 
159: raw = pd.read_excel(
160:     SUPPLEMENTARY_XLSX,
161:     sheet_name="Supplemental Dataset 2",
162:     header=None,
163:     engine="openpyxl",
164: )
165: 
166: print()
167: print("Raw shape:", raw.shape)
168: 
169: HEADER_ROW = 2
170: 
171: # Dataset 2:
172: #   column 0 = feature name
173: #   columns 1..179 = 179 retained historical samples
174: 
175: header = list(
176:     raw.iloc[
177:         HEADER_ROW,
178:         1:180,
179:     ]
180: )
181: 
182: if len(header) != 179:
183:     raise RuntimeError(
184:         f"Expected 179 Dataset 2 sample columns, found {len(header)}."
185:     )
186: 
187: 
188: # =====================================================================
189: # BUILD THE EXACT HISTORICAL SAMPLE ORDER FROM THE PMML REPRODUCTION
190: # =====================================================================
191: 
192: # The prediction file is already in the exact historical ordering
193: # established by the successful PMML reproduction.
194: 
195: historical_ids = pred["sample_key"].tolist()
196: 
197: if len(set(historical_ids)) != 179:
198:     raise RuntimeError(
199:         "Historical prediction file contains duplicate sample IDs."
200:     )
201: 
202: 
203: # =====================================================================
204: # EXTRACT THE 14 FEATURES
205: # =====================================================================
206: 
207: feature_series = (
208:     raw.iloc[:, 0]
209:     .astype(str)
210:     .str.strip()
211: )
212: 
213: rows = {}
214: 
215: for feature in FEATURES:
216: 
217:     matches = np.where(
218:         feature_series.values == feature
219:     )[0]
220: 
221:     if len(matches) != 1:
222:         raise RuntimeError(
223:             f"Expected exactly one row for {feature}, "
224:             f"found {len(matches)}."
225:         )
226: 
227:     rows[feature] = int(matches[0])
228: 
229: 
230: values = []
231: 
232: for feature in FEATURES:
233: 
234:     row_idx = rows[feature]
235: 
236:     vals = pd.to_numeric(
237:         raw.iloc[
238:             row_idx,
239:             1:180,
240:         ],
241:         errors="coerce",
242:     ).to_numpy(
243:         dtype=np.float64
244:     )
245: 
246:     if len(vals) != 179:
247:         raise RuntimeError(
248:             f"{feature} has {len(vals)} values."
249:         )
250: 
251:     if not np.all(np.isfinite(vals)):
252:         raise RuntimeError(
253:             f"{feature} contains non-finite values."
254:         )
255: 
256:     values.append(vals)
257: 
258: 
259: X = np.vstack(values).T
260: 
261: historical = pd.DataFrame(
262:     X,
263:     columns=FEATURES,
264: )
265: 
266: historical.insert(
267:     0,
268:     "sample_id",
269:     historical_ids,
270: )
271: 
272: 
273: # =====================================================================
274: # ALIGN LABELS
275: # =====================================================================
276: 
277: label_map = (
278:     pred
279:     .set_index("sample_key")["true_label"]
280:     .astype(int)
281: )
282: 
283: historical["sample_key"] = (
284:     historical["sample_id"]
285:     .map(norm_id)
286: )
287: 
288: historical["true_label"] = (
289:     historical["sample_key"]
290:     .map(label_map)
291: )
292: 
293: if historical["true_label"].isna().any():
294:     bad = historical.loc[
295:         historical["true_label"].isna(),
296:         "sample_id",
297:     ].tolist()
298: 
299:     raise RuntimeError(
300:         "Missing labels for historical samples:\n"
301:         + "\n".join(map(str, bad))
302:     )
303: 
304: historical["true_label"] = (
305:     historical["true_label"].astype(int)
306: )
307: 
308: 
309: # =====================================================================
310: # VALIDATION
311: # =====================================================================
312: 
313: y = historical["true_label"].to_numpy()
314: 
315: if len(historical) != 179:
316:     raise RuntimeError(
317:         "Historical dataset is not 179 patients."
318:     )
319: 
320: if X.shape != (179, 14):
321:     raise RuntimeError(
322:         f"Unexpected X shape: {X.shape}"
323:     )
324: 
325: if (y == 1).sum() != 97:
326:     raise RuntimeError(
327:         f"Expected 97 cancer patients, found {(y == 1).sum()}."
328:     )
329: 
330: if (y == 0).sum() != 82:
331:     raise RuntimeError(
332:         f"Expected 82 non-cancer patients, found {(y == 0).sum()}."
333:     )
334: 
335: 
336: print()
337: print("=" * 90)
338: print("HISTORICAL DATA VALIDATION")
339: print("=" * 90)
340: 
341: print()
342: print("Patients:", len(historical))
343: print("Features:", X.shape[1])
344: print("Cancer:", int((y == 1).sum()))
345: print("Non-cancer:", int((y == 0).sum()))
346: 
347: 
348: # =====================================================================
349: # 5-FOLD OOF
350: # =====================================================================
351: 
352: cv = StratifiedKFold(
353:     n_splits=5,
354:     shuffle=True,
355:     random_state=2026,
356: )
357: 
358: oof_probability = np.full(
359:     len(historical),
360:     np.nan,
361: )
362: 
363: oof_prediction = np.full(
364:     len(historical),
365:     -1,
366:     dtype=int,
367: )
368: 
369: fold_records = []
370: 
371: 
372: print()
373: print("=" * 90)
374: print("5-FOLD PATIENT-LEVEL CROSS-VALIDATION")
375: print("=" * 90)
376: 
377: 
378: for fold, (train_idx, test_idx) in enumerate(
379:     cv.split(X, y),
380:     start=1,
381: ):
382: 
383:     model = Pipeline(
384:         steps=[
385:             (
386:                 "scaler",
387:                 StandardScaler(),
388:             ),
389:             (
390:                 "classifier",
391:                 LogisticRegression(
392:                     solver="saga",
393:                     max_iter=10000,
394:                     random_state=2026,
395:                 ),
396:             ),
397:         ]
398:     )
399: 
400:     inner_cv = StratifiedKFold(
401:         n_splits=3,
402:         shuffle=True,
403:         random_state=2026 + fold,
404:     )
405: 
406:     param_grid = {
407:         "classifier__C": [
408:             0.01,
409:             0.03,
410:             0.1,
411:             0.3,
412:             1.0,
413:             3.0,
414:         ],
415:         "classifier__l1_ratio": [
416:             0.0,
417:             0.25,
418:             0.5,
419:             0.75,
420:             1.0,
421:         ],
422:         "classifier__penalty": [
423:             "elasticnet",
424:         ],
425:     }
426: 
427:     search = GridSearchCV(
428:         model,
429:         param_grid=param_grid,
430:         scoring="roc_auc",
431:         cv=inner_cv,
432:         n_jobs=-1,
433:         refit=True,
434:     )
435: 
436:     search.fit(
437:         X[train_idx],
438:         y[train_idx],
439:     )
440: 
441:     probability = search.predict_proba(
442:         X[test_idx]
443:     )[:, 1]
444: 
445:     prediction = (
446:         probability >= 0.5
447:     ).astype(int)
448: 
449:     oof_probability[test_idx] = probability
450:     oof_prediction[test_idx] = prediction
451: 
452:     fold_auc = roc_auc_score(
453:         y[test_idx],
454:         probability,
455:     )
456: 
457:     fold_accuracy = accuracy_score(
458:         y[test_idx],
459:         prediction,
460:     )
461: 
462:     fold_records.append(
463:         {
464:             "fold": fold,
465:             "n_train": int(len(train_idx)),
466:             "n_test": int(len(test_idx)),
467:             "roc_auc": float(fold_auc),
468:             "accuracy": float(fold_accuracy),
469:         }
470:     )
471: 
472:     print(
473:         f"Fold {fold}: "
474:         f"train={len(train_idx)} "
475:         f"test={len(test_idx)} "
476:         f"AUC={fold_auc:.6f} "
477:         f"accuracy={fold_accuracy:.6f}"
478:     )
479: 
480: 
481: # =====================================================================
482: # OOF METRICS
483: # =====================================================================
484: 
485: auc = roc_auc_score(
486:     y,
487:     oof_probability,
488: )
489: 
490: pr_auc = average_precision_score(
491:     y,
492:     oof_probability,
493: )
494: 
495: accuracy = accuracy_score(
496:     y,
497:     oof_prediction,
498: )
499: 
500: tn, fp, fn, tp = confusion_matrix(
501:     y,
502:     oof_prediction,
503: ).ravel()
504: 
505: sensitivity = tp / (tp + fn)
506: 
507: specificity = tn / (tn + fp)
508: 
509: 
510: # =====================================================================
511: # SAVE
512: # =====================================================================
513: 
514: historical["oof_probability_cancer"] = (
515:     oof_probability
516: )
517: 
518: historical["oof_prediction"] = (
519:     oof_prediction
520: )
521: 
522: historical["fold"] = -1
523: 
524: for fold, (_, test_idx) in enumerate(
525:     cv.split(X, y),
526:     start=1,
527: ):
528: 
529:     historical.loc[
530:         test_idx,
531:         "fold",
532:     ] = fold
533: 
534: 
535: OUTPUT_DIR.mkdir(
536:     parents=True,
537:     exist_ok=True,
538: )
539: 
540: historical.to_csv(
541:     PREDICTION_OUTPUT,
542:     index=False,
543: )
544: 
545: 
546: results = {
547:     "model": "standardized L2 logistic regression",
548:     "cohort": {
549:         "patients": 179,
550:         "cancer": 97,
551:         "non_cancer": 82,
552:     },
553:     "representation": {
554:         "source": "eLife Supplementary Dataset 2",
555:         "batch_adjusted": True,
556:         "log10_transformed": True,
557:     },
558:     "features": FEATURES,
559:     "cross_validation": {
560:         "method": "5-fold stratified patient-level OOF",
561:         "random_state": 2026,
562:     },
563:     "metrics": {
564:         "roc_auc": float(auc),
565:         "pr_auc": float(pr_auc),
566:         "accuracy": float(accuracy),
567:         "sensitivity": float(sensitivity),
568:         "specificity": float(specificity),
569:     },
570:     "confusion_matrix": {
571:         "TN": int(tn),
572:         "FP": int(fp),
573:         "FN": int(fn),
574:         "TP": int(tp),
575:     },
576:     "folds": fold_records,
577: }
578: 
579: with RESULTS_OUTPUT.open("w") as f:
580: 
581:     json.dump(
582:         results,
583:         f,
584:         indent=2,
585:     )
586: 
587: 
588: print()
589: print("=" * 90)
590: print("FAIR HISTORICAL 14-miRNA BASELINE RESULTS")
591: print("=" * 90)
592: 
593: print()
594: print(f"Patients:       {len(historical)}")
595: print(f"Cancer:         {(y == 1).sum()}")
596: print(f"Non-cancer:     {(y == 0).sum()}")
597: 
598: print()
599: print(f"OOF ROC AUC:    {auc:.6f}")
600: print(f"OOF PR AUC:     {pr_auc:.6f}")
601: print(f"OOF Accuracy:   {accuracy:.6f}")
602: print(f"Sensitivity:    {sensitivity:.6f}")
603: print(f"Specificity:    {specificity:.6f}")
604: 
605: print()
606: print("Confusion matrix:")
607: print(f"  TP={tp} FN={fn}")
608: print(f"  FP={fp} TN={tn}")
609: 
610: print()
611: print("Saved:")
612: print(PREDICTION_OUTPUT)
613: print(RESULTS_OUTPUT)
~~~

## scripts/build_full_historical_mirna_baseline.py

~~~python
1: #!/usr/bin/env python3
2: 
3: """
4: FULL HISTORICAL miRNA BASELINE
5: 
6: Uses the full miRNA feature set available in eLife Supplemental
7: Dataset 2 on the exact historical 179-patient cohort.
8: 
9: This is a NEW model.
10: 
11: It does NOT use:
12:     - historical PMML weights
13:     - historical neural-network architecture
14: 
15: Model:
16:     StandardScaler
17:         ->
18:     Elastic-net logistic regression
19: 
20: Evaluation:
21:     5-fold outer patient-level stratified CV
22: 
23: Hyperparameter tuning:
24:     3-fold inner CV within each outer training set
25: 
26: Therefore:
27: 
28:     test patients
29:         are never used for
30:         feature/model selection.
31: 
32: This prevents data leakage.
33: 
34: The purpose is to compare:
35: 
36:     historical 14-miRNA neural network
37:                 vs
38:     modern 14-miRNA logistic regression
39:                 vs
40:     modern full-miRNA regularized model
41: """
42: 
43: from __future__ import annotations
44: 
45: import json
46: import re
47: from pathlib import Path
48: 
49: import numpy as np
50: import pandas as pd
51: 
52: from sklearn.linear_model import LogisticRegression
53: from sklearn.metrics import (
54:     accuracy_score,
55:     average_precision_score,
56:     confusion_matrix,
57:     roc_auc_score,
58: )
59: from sklearn.model_selection import (
60:     GridSearchCV,
61:     StratifiedKFold,
62: )
63: from sklearn.pipeline import Pipeline
64: from sklearn.preprocessing import StandardScaler
65: 
66: 
67: # =====================================================================
68: # PATHS
69: # =====================================================================
70: 
71: ROOT = Path.home() / "ovarian-ai"
72: 
73: SUPPLEMENTARY_XLSX = (
74:     ROOT
75:     / "data"
76:     / "external"
77:     / "elife_28932"
78:     / "elife_poa_e28932_Supplementary_File_6.xlsx"
79: )
80: 
81: HISTORICAL_PREDICTIONS = (
82:     ROOT
83:     / "reports"
84:     / "tables"
85:     / "original_historical_clean_predictions.csv"
86: )
87: 
88: OUTPUT_DIR = (
89:     ROOT
90:     / "reports"
91:     / "tables"
92: )
93: 
94: PREDICTIONS_OUTPUT = (
95:     OUTPUT_DIR
96:     / "full_historical_mirna_oof_predictions.csv"
97: )
98: 
99: RESULTS_OUTPUT = (
100:     OUTPUT_DIR
101:     / "full_historical_mirna_baseline_results.json"
102: )
103: 
104: 
105: # =====================================================================
106: # SETTINGS
107: # =====================================================================
108: 
109: OUTER_SPLITS = 5
110: INNER_SPLITS = 3
111: 
112: RANDOM_STATE = 2026
113: 
114: 
115: # =====================================================================
116: # ID NORMALIZATION
117: # =====================================================================
118: 
119: def norm_id(value) -> str:
120: 
121:     s = str(value).strip().lower()
122: 
123:     if re.fullmatch(r"\d+\.0+", s):
124:         s = s.split(".", 1)[0]
125: 
126:     return s.replace("-", "")
127: 
128: 
129: # =====================================================================
130: # IDENTIFY miRNA FEATURE ROWS
131: # =====================================================================
132: 
133: def looks_like_mirna(value) -> bool:
134: 
135:     if pd.isna(value):
136:         return False
137: 
138:     s = str(value).strip()
139: 
140:     return (
141:         s.startswith("hsa-miR-")
142:         or s.startswith("hsa-let-")
143:     )
144: 
145: 
146: # =====================================================================
147: # LOAD HISTORICAL COHORT
148: # =====================================================================
149: 
150: print("=" * 90)
151: print("FULL HISTORICAL miRNA BASELINE")
152: print("=" * 90)
153: 
154: print()
155: print("Loading historical cohort...")
156: 
157: pred = pd.read_csv(
158:     HISTORICAL_PREDICTIONS
159: )
160: 
161: if len(pred) != 179:
162:     raise RuntimeError(
163:         f"Expected 179 historical patients, found {len(pred)}."
164:     )
165: 
166: pred["sample_key"] = (
167:     pred["sample_id"].map(norm_id)
168: )
169: 
170: if pred["sample_key"].duplicated().any():
171:     raise RuntimeError(
172:         "Duplicate historical sample IDs detected."
173:     )
174: 
175: print(
176:     f"Patients: {len(pred)}"
177: )
178: 
179: print(
180:     f"Cancer: {(pred['true_label'] == 1).sum()}"
181: )
182: 
183: print(
184:     f"Non-cancer: {(pred['true_label'] == 0).sum()}"
185: )
186: 
187: 
188: # =====================================================================
189: # READ DATASET 2
190: # =====================================================================
191: 
192: print()
193: print("=" * 90)
194: print("READING eLIFE SUPPLEMENTAL DATASET 2")
195: print("=" * 90)
196: 
197: raw = pd.read_excel(
198:     SUPPLEMENTARY_XLSX,
199:     sheet_name="Supplemental Dataset 2",
200:     header=None,
201:     engine="openpyxl",
202: )
203: 
204: print()
205: print("Raw shape:", raw.shape)
206: 
207: HEADER_ROW = 2
208: 
209: raw_header = list(
210:     raw.iloc[
211:         HEADER_ROW,
212:         1:180,
213:     ]
214: )
215: 
216: if len(raw_header) != 179:
217:     raise RuntimeError(
218:         f"Expected 179 sample columns, found {len(raw_header)}."
219:     )
220: 
221: 
222: # =====================================================================
223: # USE THE SAME HISTORICAL SAMPLE ORDER AS THE PMML REPRODUCTION
224: # =====================================================================
225: 
226: historical_ids = pred["sample_key"].tolist()
227: 
228: # Dataset 2's first 119 numeric columns correspond to the retained
229: # numeric historical samples; the remainder are the 60 e-samples.
230: #
231: # We already independently established the exact historical ordering
232: # in original_historical_clean_predictions.csv.
233: 
234: if len(historical_ids) != 179:
235:     raise RuntimeError(
236:         "Historical ID count is not 179."
237:     )
238: 
239: 
240: # =====================================================================
241: # FIND ALL miRNA ROWS
242: # =====================================================================
243: 
244: feature_names = []
245: 
246: feature_rows = []
247: 
248: for row_idx, value in enumerate(
249:     raw.iloc[:, 0]
250: ):
251: 
252:     if not looks_like_mirna(value):
253:         continue
254: 
255:     feature = str(value).strip()
256: 
257:     # Avoid duplicate feature names.
258:     if feature in feature_names:
259:         continue
260: 
261:     feature_names.append(feature)
262:     feature_rows.append(row_idx)
263: 
264: 
265: if len(feature_names) < 50:
266:     raise RuntimeError(
267:         "Unexpectedly few miRNA features found: "
268:         f"{len(feature_names)}"
269:     )
270: 
271: 
272: print()
273: print("=" * 90)
274: print("FULL miRNA FEATURE SET")
275: print("=" * 90)
276: 
277: print()
278: print(
279:     f"miRNA features detected: {len(feature_names)}"
280: )
281: 
282: print()
283: print("First 20 features:")
284: 
285: for feature in feature_names[:20]:
286:     print(f"  {feature}")
287: 
288: print()
289: print("Last 10 features:")
290: 
291: for feature in feature_names[-10:]:
292:     print(f"  {feature}")
293: 
294: 
295: # =====================================================================
296: # EXTRACT MATRIX
297: # =====================================================================
298: 
299: values = []
300: 
301: for feature, row_idx in zip(
302:     feature_names,
303:     feature_rows,
304: ):
305: 
306:     vals = pd.to_numeric(
307:         raw.iloc[
308:             row_idx,
309:             1:180,
310:         ],
311:         errors="coerce",
312:     ).to_numpy(
313:         dtype=np.float64
314:     )
315: 
316:     if len(vals) != 179:
317:         raise RuntimeError(
318:             f"{feature} has {len(vals)} values."
319:         )
320: 
321:     if not np.all(np.isfinite(vals)):
322:         bad = np.where(
323:             ~np.isfinite(vals)
324:         )[0]
325: 
326:         raise RuntimeError(
327:             f"{feature} contains non-finite values "
328:             f"at columns {bad[:10].tolist()}."
329:         )
330: 
331:     values.append(vals)
332: 
333: 
334: X = np.vstack(values).T
335: 
336: 
337: # =====================================================================
338: # BUILD DATAFRAME
339: # =====================================================================
340: 
341: mirna = pd.DataFrame(
342:     X,
343:     columns=feature_names,
344: )
345: 
346: mirna.insert(
347:     0,
348:     "sample_id",
349:     historical_ids,
350: )
351: 
352: mirna["sample_key"] = (
353:     mirna["sample_id"]
354:     .map(norm_id)
355: )
356: 
357: 
358: # =====================================================================
359: # ALIGN LABELS
360: # =====================================================================
361: 
362: label_map = (
363:     pred
364:     .set_index("sample_key")["true_label"]
365:     .astype(int)
366: )
367: 
368: mirna["true_label"] = (
369:     mirna["sample_key"]
370:     .map(label_map)
371: )
372: 
373: if mirna["true_label"].isna().any():
374: 
375:     bad = mirna.loc[
376:         mirna["true_label"].isna(),
377:         "sample_id",
378:     ].tolist()
379: 
380:     raise RuntimeError(
381:         "Missing labels for:\n"
382:         + "\n".join(
383:             map(str, bad)
384:         )
385:     )
386: 
387: mirna["true_label"] = (
388:     mirna["true_label"].astype(int)
389: )
390: 
391: 
392: # =====================================================================
393: # VALIDATION
394: # =====================================================================
395: 
396: y = mirna["true_label"].to_numpy()
397: 
398: if X.shape[0] != 179:
399:     raise RuntimeError(
400:         f"Expected 179 patients, found {X.shape[0]}."
401:     )
402: 
403: if X.shape[1] != len(feature_names):
404:     raise RuntimeError(
405:         "Feature dimension mismatch."
406:     )
407: 
408: if (y == 1).sum() != 97:
409:     raise RuntimeError(
410:         f"Expected 97 cancer patients, found {(y == 1).sum()}."
411:     )
412: 
413: if (y == 0).sum() != 82:
414:     raise RuntimeError(
415:         f"Expected 82 non-cancer patients, found {(y == 0).sum()}."
416:     )
417: 
418: 
419: print()
420: print("=" * 90)
421: print("DATA VALIDATION")
422: print("=" * 90)
423: 
424: print()
425: print("Patients:", X.shape[0])
426: print("miRNA features:", X.shape[1])
427: print("Cancer:", int((y == 1).sum()))
428: print("Non-cancer:", int((y == 0).sum()))
429: 
430: 
431: # =====================================================================
432: # OUTER CV
433: # =====================================================================
434: 
435: outer_cv = StratifiedKFold(
436:     n_splits=OUTER_SPLITS,
437:     shuffle=True,
438:     random_state=RANDOM_STATE,
439: )
440: 
441: oof_probability = np.full(
442:     len(y),
443:     np.nan,
444: )
445: 
446: oof_prediction = np.full(
447:     len(y),
448:     -1,
449:     dtype=int,
450: )
451: 
452: fold_records = []
453: 
454: 
455: print()
456: print("=" * 90)
457: print("NESTED CROSS-VALIDATION")
458: print("=" * 90)
459: 
460: 
461: for fold, (train_idx, test_idx) in enumerate(
462:     outer_cv.split(X, y),
463:     start=1,
464: ):
465: 
466:     print()
467:     print(
468:         f"OUTER FOLD {fold}/{OUTER_SPLITS}"
469:     )
470: 
471:     X_train = X[train_idx]
472:     X_test = X[test_idx]
473: 
474:     y_train = y[train_idx]
475:     y_test = y[test_idx]
476: 
477: 
478:     # ---------------------------------------------------------------
479:     # Inner CV
480:     # ---------------------------------------------------------------
481: 
482:     inner_cv = StratifiedKFold(
483:         n_splits=INNER_SPLITS,
484:         shuffle=True,
485:         random_state=RANDOM_STATE + fold,
486:     )
487: 
488: 
489:     # ---------------------------------------------------------------
490:     # Pipeline
491:     #
492:     # Scaling occurs separately inside every inner/outer training
493:     # split.
494:     # ---------------------------------------------------------------
495: 
496:     pipeline = Pipeline(
497:         steps=[
498:             (
499:                 "scaler",
500:                 StandardScaler(),
501:             ),
502:             (
503:                 "classifier",
504:                 LogisticRegression(
505:                     solver="saga",
506:                     max_iter=10000,
507:                     random_state=RANDOM_STATE,
508:                 ),
509:             ),
510:         ]
511:     )
512: 
513: 
514:     # ---------------------------------------------------------------
515:     # Hyperparameter grid
516:     # ---------------------------------------------------------------
517: 
518:     param_grid = {
519:         "classifier__C": [
520:             0.01,
521:             0.03,
522:             0.1,
523:             0.3,
524:             1.0,
525:             3.0,
526:         ],
527:         "classifier__l1_ratio": [
528:             0.0,
529:             0.25,
530:             0.5,
531:             0.75,
532:             1.0,
533:         ],
534:         "classifier__penalty": [
535:             "elasticnet",
536:         ],
537:     }
538: 
539: 
540:     search = GridSearchCV(
541:         pipeline,
542:         param_grid=param_grid,
543:         scoring="roc_auc",
544:         cv=inner_cv,
545:         n_jobs=-1,
546:         refit=True,
547:     )
548: 
549: 
550:     search.fit(
551:         X_train,
552:         y_train,
553:     )
554: 
555: 
556:     probability = search.predict_proba(
557:         X_test
558:     )[:, 1]
559: 
560:     prediction = (
561:         probability >= 0.5
562:     ).astype(int)
563: 
564: 
565:     oof_probability[test_idx] = (
566:         probability
567:     )
568: 
569:     oof_prediction[test_idx] = (
570:         prediction
571:     )
572: 
573: 
574:     fold_auc = roc_auc_score(
575:         y_test,
576:         probability,
577:     )
578: 
579:     fold_accuracy = accuracy_score(
580:         y_test,
581:         prediction,
582:     )
583: 
584: 
585:     best_params = {
586:         key: value
587:         for key, value in search.best_params_.items()
588:     }
589: 
590: 
591:     fold_records.append(
592:         {
593:             "fold": fold,
594:             "n_train": int(len(train_idx)),
595:             "n_test": int(len(test_idx)),
596:             "roc_auc": float(fold_auc),
597:             "accuracy": float(fold_accuracy),
598:             "inner_cv_auc": float(
599:                 search.best_score_
600:             ),
601:             "best_params": best_params,
602:         }
603:     )
604: 
605: 
606:     print(
607:         f"  test patients: {len(test_idx)}"
608:     )
609: 
610:     print(
611:         f"  inner-CV AUC: {search.best_score_:.6f}"
612:     )
613: 
614:     print(
615:         f"  test AUC:     {fold_auc:.6f}"
616:     )
617: 
618:     print(
619:         f"  accuracy:     {fold_accuracy:.6f}"
620:     )
621: 
622:     print(
623:         f"  best params:  {best_params}"
624:     )
625: 
626: 
627: # =====================================================================
628: # OOF METRICS
629: # =====================================================================
630: 
631: if np.any(
632:     ~np.isfinite(oof_probability)
633: ):
634:     raise RuntimeError(
635:         "Some patients lack OOF probabilities."
636:     )
637: 
638: if np.any(
639:     oof_prediction < 0
640: ):
641:     raise RuntimeError(
642:         "Some patients lack OOF predictions."
643:     )
644: 
645: 
646: auc = roc_auc_score(
647:     y,
648:     oof_probability,
649: )
650: 
651: pr_auc = average_precision_score(
652:     y,
653:     oof_probability,
654: )
655: 
656: accuracy = accuracy_score(
657:     y,
658:     oof_prediction,
659: )
660: 
661: tn, fp, fn, tp = confusion_matrix(
662:     y,
663:     oof_prediction,
664: ).ravel()
665: 
666: sensitivity = tp / (tp + fn)
667: 
668: specificity = tn / (tn + fp)
669: 
670: 
671: # =====================================================================
672: # SAVE PREDICTIONS
673: # =====================================================================
674: 
675: mirna["oof_probability_cancer"] = (
676:     oof_probability
677: )
678: 
679: mirna["oof_prediction"] = (
680:     oof_prediction
681: )
682: 
683: mirna["fold"] = -1
684: 
685: 
686: for fold, (_, test_idx) in enumerate(
687:     outer_cv.split(X, y),
688:     start=1,
689: ):
690: 
691:     mirna.loc[
692:         test_idx,
693:         "fold",
694:     ] = fold
695: 
696: 
697: OUTPUT_DIR.mkdir(
698:     parents=True,
699:     exist_ok=True,
700: )
701: 
702: mirna.to_csv(
703:     PREDICTIONS_OUTPUT,
704:     index=False,
705: )
706: 
707: 
708: # =====================================================================
709: # SAVE RESULTS
710: # =====================================================================
711: 
712: results = {
713:     "model": "nested-CV elastic-net logistic regression",
714:     "representation": {
715:         "source": "eLife Supplemental Dataset 2",
716:         "batch_adjusted": True,
717:         "log10_transformed": True,
718:     },
719:     "cohort": {
720:         "patients": 179,
721:         "cancer": 97,
722:         "non_cancer": 82,
723:     },
724:     "n_mirna_features": int(
725:         len(feature_names)
726:     ),
727:     "cross_validation": {
728:         "outer": "5-fold stratified patient-level",
729:         "inner": "3-fold stratified",
730:         "random_state": RANDOM_STATE,
731:     },
732:     "metrics": {
733:         "roc_auc": float(auc),
734:         "pr_auc": float(pr_auc),
735:         "accuracy": float(accuracy),
736:         "sensitivity": float(sensitivity),
737:         "specificity": float(specificity),
738:     },
739:     "confusion_matrix": {
740:         "TN": int(tn),
741:         "FP": int(fp),
742:         "FN": int(fn),
743:         "TP": int(tp),
744:     },
745:     "folds": fold_records,
746: }
747: 
748: with RESULTS_OUTPUT.open(
749:     "w"
750: ) as f:
751: 
752:     json.dump(
753:         results,
754:         f,
755:         indent=2,
756:     )
757: 
758: 
759: # =====================================================================
760: # REPORT
761: # =====================================================================
762: 
763: print()
764: print("=" * 90)
765: print("FULL HISTORICAL miRNA BASELINE RESULTS")
766: print("=" * 90)
767: 
768: print()
769: print(f"Patients:       {len(y)}")
770: print(f"miRNA features: {len(feature_names)}")
771: 
772: print()
773: print(f"Cancer:         {(y == 1).sum()}")
774: print(f"Non-cancer:     {(y == 0).sum()}")
775: 
776: print()
777: print(f"OOF ROC AUC:    {auc:.6f}")
778: print(f"OOF PR AUC:     {pr_auc:.6f}")
779: print(f"OOF Accuracy:   {accuracy:.6f}")
780: print(f"Sensitivity:    {sensitivity:.6f}")
781: print(f"Specificity:    {specificity:.6f}")
782: 
783: print()
784: print("Confusion matrix:")
785: print(f"  TP={tp} FN={fn}")
786: print(f"  FP={fp} TN={tn}")
787: 
788: print()
789: print("Saved:")
790: print(PREDICTIONS_OUTPUT)
791: print(RESULTS_OUTPUT)
792: 
793: print()
794: print("=" * 90)
795: print("IMPORTANT")
796: print("=" * 90)
797: 
798: print()
799: print(
800:     "This is a NEW model and is not the historical eLife neural network."
801: )
802: 
803: print(
804:     "All scaling and hyperparameter selection were performed "
805:     "inside the cross-validation procedure."
806: )
~~~

## scripts/build_gse106817_170mirna_baseline.py

~~~python
1: #!/usr/bin/env python3
2: 
3: from pathlib import Path
4: import json
5: 
6: import numpy as np
7: import pandas as pd
8: 
9: from sklearn.linear_model import LogisticRegression
10: from sklearn.metrics import (
11:     accuracy_score,
12:     average_precision_score,
13:     confusion_matrix,
14:     roc_auc_score,
15: )
16: from sklearn.model_selection import StratifiedKFold
17: from sklearn.pipeline import Pipeline
18: from sklearn.preprocessing import StandardScaler
19: 
20: 
21: ROOT = Path.home() / "ovarian-ai"
22: 
23: DATA = (
24:     ROOT
25:     / "data"
26:     / "processed"
27:     / "gse106817_external_validation.csv.gz"
28: )
29: 
30: OUT_PRED = (
31:     ROOT
32:     / "reports"
33:     / "tables"
34:     / "gse106817_170mirna_oof_predictions.csv"
35: )
36: 
37: OUT_RESULTS = (
38:     ROOT
39:     / "reports"
40:     / "tables"
41:     / "gse106817_170mirna_baseline_results.json"
42: )
43: 
44: 
45: print("=" * 80)
46: print("GSE106817 INDEPENDENT 170-miRNA BASELINE")
47: print("=" * 80)
48: 
49: df = pd.read_csv(DATA)
50: 
51: feature_cols = [
52:     c for c in df.columns
53:     if c.startswith("hsa-")
54: ]
55: 
56: X = df[feature_cols].apply(
57:     pd.to_numeric,
58:     errors="coerce",
59: )
60: 
61: y = df["binary_label"].astype(int)
62: 
63: print()
64: print("Patients:", len(df))
65: print("Features:", len(feature_cols))
66: print("Cancer:", int(y.sum()))
67: print("Non-cancer:", int((y == 0).sum()))
68: 
69: # ---------------------------------------------------------------
70: # Remove rows with missing molecular data.
71: # ---------------------------------------------------------------
72: 
73: valid = (
74:     X.notna().all(axis=1)
75:     & y.notna()
76: )
77: 
78: X = X.loc[valid].reset_index(drop=True)
79: y = y.loc[valid].reset_index(drop=True)
80: meta = df.loc[valid].reset_index(drop=True)
81: 
82: print()
83: print("Patients after feature validation:", len(X))
84: print("Missing miRNA rows:", int((~valid).sum()))
85: 
86: # ---------------------------------------------------------------
87: # Model
88: #
89: # Class weighting is important here because GSE106817 has many
90: # non-cancer controls.
91: # ---------------------------------------------------------------
92: 
93: model = Pipeline([
94:     (
95:         "scaler",
96:         StandardScaler(),
97:     ),
98:     (
99:         "classifier",
100:         LogisticRegression(
101:             max_iter=5000,
102:             class_weight="balanced",
103:             solver="liblinear",
104:             random_state=42,
105:         ),
106:     ),
107: ])
108: 
109: cv = StratifiedKFold(
110:     n_splits=5,
111:     shuffle=True,
112:     random_state=42,
113: )
114: 
115: oof_prob = np.full(
116:     len(X),
117:     np.nan,
118: )
119: 
120: oof_pred = np.full(
121:     len(X),
122:     np.nan,
123: )
124: 
125: print()
126: print("=" * 80)
127: print("5-FOLD PATIENT-LEVEL OOF VALIDATION")
128: print("=" * 80)
129: 
130: for fold, (train_idx, test_idx) in enumerate(
131:     cv.split(X, y),
132:     start=1,
133: ):
134: 
135:     X_train = X.iloc[train_idx]
136:     X_test = X.iloc[test_idx]
137: 
138:     y_train = y.iloc[train_idx]
139: 
140:     model.fit(
141:         X_train,
142:         y_train,
143:     )
144: 
145:     prob = model.predict_proba(
146:         X_test
147:     )[:, 1]
148: 
149:     pred = (
150:         prob >= 0.5
151:     ).astype(int)
152: 
153:     oof_prob[test_idx] = prob
154:     oof_pred[test_idx] = pred
155: 
156:     auc = roc_auc_score(
157:         y.iloc[test_idx],
158:         prob,
159:     )
160: 
161:     acc = accuracy_score(
162:         y.iloc[test_idx],
163:         pred,
164:     )
165: 
166:     print(
167:         f"Fold {fold}: "
168:         f"train={len(train_idx)} "
169:         f"test={len(test_idx)} "
170:         f"AUC={auc:.6f} "
171:         f"accuracy={acc:.6f}"
172:     )
173: 
174: # ---------------------------------------------------------------
175: # Overall OOF metrics
176: # ---------------------------------------------------------------
177: 
178: auc = roc_auc_score(
179:     y,
180:     oof_prob,
181: )
182: 
183: pr_auc = average_precision_score(
184:     y,
185:     oof_prob,
186: )
187: 
188: accuracy = accuracy_score(
189:     y,
190:     oof_pred,
191: )
192: 
193: cm = confusion_matrix(
194:     y,
195:     oof_pred,
196: )
197: 
198: tn, fp, fn, tp = cm.ravel()
199: 
200: sensitivity = (
201:     tp / (tp + fn)
202:     if (tp + fn)
203:     else np.nan
204: )
205: 
206: specificity = (
207:     tn / (tn + fp)
208:     if (tn + fp)
209:     else np.nan
210: )
211: 
212: print()
213: print("=" * 80)
214: print("RESULTS")
215: print("=" * 80)
216: 
217: print(
218:     f"OOF ROC AUC:    {auc:.6f}"
219: )
220: 
221: print(
222:     f"OOF PR AUC:     {pr_auc:.6f}"
223: )
224: 
225: print(
226:     f"OOF Accuracy:   {accuracy:.6f}"
227: )
228: 
229: print(
230:     f"Sensitivity:    {sensitivity:.6f}"
231: )
232: 
233: print(
234:     f"Specificity:    {specificity:.6f}"
235: )
236: 
237: print()
238: print("Confusion matrix:")
239: print(
240:     f"  TP={tp} FN={fn}"
241: )
242: print(
243:     f"  FP={fp} TN={tn}"
244: )
245: 
246: # ---------------------------------------------------------------
247: # Save OOF predictions
248: # ---------------------------------------------------------------
249: 
250: predictions = meta[
251:     [
252:         "sample_id",
253:         "title",
254:         "group",
255:         "age",
256:         "binary_label",
257:     ]
258: ].copy()
259: 
260: predictions["prediction_probability"] = (
261:     oof_prob
262: )
263: 
264: predictions["prediction"] = (
265:     oof_pred.astype(int)
266: )
267: 
268: predictions["correct"] = (
269:     predictions["prediction"]
270:     == predictions["binary_label"]
271: )
272: 
273: OUT_PRED.parent.mkdir(
274:     parents=True,
275:     exist_ok=True,
276: )
277: 
278: predictions.to_csv(
279:     OUT_PRED,
280:     index=False,
281: )
282: 
283: results = {
284:     "dataset": "GSE106817",
285:     "patients": int(len(X)),
286:     "features": int(len(feature_cols)),
287:     "cancer": int(y.sum()),
288:     "non_cancer": int((y == 0).sum()),
289:     "roc_auc": float(auc),
290:     "pr_auc": float(pr_auc),
291:     "accuracy": float(accuracy),
292:     "sensitivity": float(sensitivity),
293:     "specificity": float(specificity),
294:     "tp": int(tp),
295:     "fn": int(fn),
296:     "fp": int(fp),
297:     "tn": int(tn),
298:     "model": "balanced logistic regression",
299:     "evaluation": "5-fold stratified patient-level OOF",
300: }
301: 
302: OUT_RESULTS.write_text(
303:     json.dumps(
304:         results,
305:         indent=2,
306:     )
307: )
308: 
309: print()
310: print("=" * 80)
311: print("FILES")
312: print("=" * 80)
313: 
314: print(OUT_PRED)
315: print(OUT_RESULTS)
316: 
317: print()
318: print("IMPORTANT:")
319: print(
320:     "This is an independent GSE106817 baseline."
321: )
322: print(
323:     "It is NOT the frozen GSE94533 model."
324: )
325: print(
326:     "Age is not used because most GSE106817 controls lack age."
327: )
~~~

## scripts/build_gse106817_allmirna_borderline_discovery.py

~~~python
1: #!/usr/bin/env python3
2: 
3: from pathlib import Path
4: import json
5: 
6: import numpy as np
7: import pandas as pd
8: 
9: from sklearn.feature_selection import SelectKBest, f_classif
10: from sklearn.linear_model import LogisticRegression
11: from sklearn.metrics import (
12:     accuracy_score,
13:     average_precision_score,
14:     confusion_matrix,
15:     roc_auc_score,
16: )
17: from sklearn.model_selection import GridSearchCV, StratifiedKFold
18: from sklearn.pipeline import Pipeline
19: from sklearn.preprocessing import StandardScaler
20: 
21: 
22: ROOT = Path.home() / "ovarian-ai"
23: 
24: DATA = (
25:     ROOT
26:     / "data"
27:     / "processed"
28:     / "gse106817_external_validation.csv.gz"
29: )
30: 
31: OUT_PRED = (
32:     ROOT
33:     / "reports"
34:     / "tables"
35:     / "gse106817_allmirna_borderline_oof_predictions.csv"
36: )
37: 
38: OUT_RESULTS = (
39:     ROOT
40:     / "reports"
41:     / "tables"
42:     / "gse106817_allmirna_borderline_results.json"
43: )
44: 
45: 
46: print("=" * 80)
47: print("GSE106817 ALL-miRNA CANCER vs BORDERLINE DISCOVERY")
48: print("=" * 80)
49: 
50: df = pd.read_csv(DATA)
51: 
52: # ---------------------------------------------------------------
53: # Only cancer and borderline.
54: # ---------------------------------------------------------------
55: 
56: df = df[
57:     df["group"].isin([
58:         "cancer",
59:         "borderline",
60:     ])
61: ].copy()
62: 
63: # ---------------------------------------------------------------
64: # Every available miRNA in GSE106817.
65: # ---------------------------------------------------------------
66: 
67: mirna_cols = [
68:     c for c in df.columns
69:     if c.startswith("hsa-")
70: ]
71: 
72: feature_cols = mirna_cols + ["age"]
73: 
74: X = df[feature_cols].apply(
75:     pd.to_numeric,
76:     errors="coerce",
77: )
78: 
79: y = (
80:     df["group"] == "cancer"
81: ).astype(int)
82: 
83: # ---------------------------------------------------------------
84: # Remove incomplete rows.
85: # ---------------------------------------------------------------
86: 
87: valid = X.notna().all(axis=1)
88: 
89: X = X.loc[valid].reset_index(drop=True)
90: y = y.loc[valid].reset_index(drop=True)
91: meta = df.loc[valid].reset_index(drop=True)
92: 
93: print()
94: print("Patients:", len(X))
95: print("Cancer:", int(y.sum()))
96: print("Borderline:", int((y == 0).sum()))
97: print("miRNAs:", len(mirna_cols))
98: print("Total features including age:", len(feature_cols))
99: 
100: # ---------------------------------------------------------------
101: # Nested CV
102: #
103: # Feature selection happens ONLY inside the training folds.
104: # ---------------------------------------------------------------
105: 
106: outer_cv = StratifiedKFold(
107:     n_splits=5,
108:     shuffle=True,
109:     random_state=42,
110: )
111: 
112: inner_cv = StratifiedKFold(
113:     n_splits=4,
114:     shuffle=True,
115:     random_state=43,
116: )
117: 
118: pipeline = Pipeline([
119:     (
120:         "scale",
121:         StandardScaler(),
122:     ),
123:     (
124:         "select",
125:         SelectKBest(
126:             score_func=f_classif,
127:         ),
128:     ),
129:     (
130:         "clf",
131:         LogisticRegression(
132:             max_iter=5000,
133:             class_weight="balanced",
134:             solver="liblinear",
135:         ),
136:     ),
137: ])
138: 
139: param_grid = {
140:     "select__k": [
141:         5,
142:         10,
143:         20,
144:         30,
145:         50,
146:         75,
147:         100,
148:         150,
149:         200,
150:         300,
151:     ],
152:     "clf__C": [
153:         0.01,
154:         0.03,
155:         0.1,
156:         0.3,
157:         1.0,
158:         3.0,
159:     ],
160: }
161: 
162: oof_prob = np.full(
163:     len(X),
164:     np.nan,
165: )
166: 
167: oof_pred = np.full(
168:     len(X),
169:     np.nan,
170: )
171: 
172: fold_results = []
173: 
174: for fold, (train_idx, test_idx) in enumerate(
175:     outer_cv.split(X, y),
176:     start=1,
177: ):
178: 
179:     print()
180:     print("=" * 80)
181:     print(f"OUTER FOLD {fold}/5")
182:     print("=" * 80)
183: 
184:     search = GridSearchCV(
185:         estimator=pipeline,
186:         param_grid=param_grid,
187:         cv=inner_cv,
188:         scoring="roc_auc",
189:         n_jobs=-1,
190:         refit=True,
191:     )
192: 
193:     search.fit(
194:         X.iloc[train_idx],
195:         y.iloc[train_idx],
196:     )
197: 
198:     probability = search.predict_proba(
199:         X.iloc[test_idx]
200:     )[:, 1]
201: 
202:     prediction = (
203:         probability >= 0.5
204:     ).astype(int)
205: 
206:     oof_prob[test_idx] = probability
207:     oof_pred[test_idx] = prediction
208: 
209:     auc = roc_auc_score(
210:         y.iloc[test_idx],
211:         probability,
212:     )
213: 
214:     accuracy = accuracy_score(
215:         y.iloc[test_idx],
216:         prediction,
217:     )
218: 
219:     fold_results.append({
220:         "fold": fold,
221:         "test_n": int(len(test_idx)),
222:         "auc": float(auc),
223:         "accuracy": float(accuracy),
224:         "best_params": search.best_params_,
225:     })
226: 
227:     print(
228:         f"test={len(test_idx)} "
229:         f"AUC={auc:.6f} "
230:         f"accuracy={accuracy:.6f}"
231:     )
232: 
233:     print(
234:         "best:",
235:         search.best_params_,
236:     )
237: 
238: # ---------------------------------------------------------------
239: # Overall OOF metrics
240: # ---------------------------------------------------------------
241: 
242: auc = roc_auc_score(
243:     y,
244:     oof_prob,
245: )
246: 
247: pr_auc = average_precision_score(
248:     y,
249:     oof_prob,
250: )
251: 
252: accuracy = accuracy_score(
253:     y,
254:     oof_pred,
255: )
256: 
257: tn, fp, fn, tp = confusion_matrix(
258:     y,
259:     oof_pred,
260: ).ravel()
261: 
262: sensitivity = tp / (tp + fn)
263: specificity = tn / (tn + fp)
264: 
265: print()
266: print("=" * 80)
267: print("ALL-miRNA DISCOVERY RESULTS")
268: print("=" * 80)
269: 
270: print(f"ROC AUC:     {auc:.6f}")
271: print(f"PR AUC:      {pr_auc:.6f}")
272: print(f"Accuracy:    {accuracy:.6f}")
273: print(f"Sensitivity: {sensitivity:.6f}")
274: print(f"Specificity: {specificity:.6f}")
275: 
276: print()
277: print("Confusion matrix:")
278: print(f"  TN={tn} FP={fp}")
279: print(f"  FN={fn} TP={tp}")
280: 
281: # ---------------------------------------------------------------
282: # Save predictions
283: # ---------------------------------------------------------------
284: 
285: pred = meta[
286:     [
287:         "sample_id",
288:         "title",
289:         "group",
290:         "age",
291:     ]
292: ].copy()
293: 
294: pred["true_label"] = y.values
295: pred["cancer_probability"] = oof_prob
296: pred["predicted_label"] = oof_pred.astype(int)
297: 
298: OUT_PRED.parent.mkdir(
299:     parents=True,
300:     exist_ok=True,
301: )
302: 
303: pred.to_csv(
304:     OUT_PRED,
305:     index=False,
306: )
307: 
308: results = {
309:     "dataset": "GSE106817",
310:     "endpoint": "cancer_vs_borderline",
311:     "patients": int(len(X)),
312:     "cancer": int(y.sum()),
313:     "borderline": int((y == 0).sum()),
314:     "mirnas": int(len(mirna_cols)),
315:     "features_with_age": int(len(feature_cols)),
316:     "roc_auc": float(auc),
317:     "pr_auc": float(pr_auc),
318:     "accuracy": float(accuracy),
319:     "sensitivity": float(sensitivity),
320:     "specificity": float(specificity),
321:     "tn": int(tn),
322:     "fp": int(fp),
323:     "fn": int(fn),
324:     "tp": int(tp),
325:     "feature_selection_inside_training_folds": True,
326:     "outer_cv": "5-fold stratified",
327:     "inner_cv": "4-fold stratified",
328:     "fold_results": fold_results,
329: }
330: 
331: OUT_RESULTS.write_text(
332:     json.dumps(
333:         results,
334:         indent=2,
335:     )
336: )
337: 
338: print()
339: print("Saved:")
340: print(OUT_PRED)
341: print(OUT_RESULTS)
342: 
343: print()
344: print("IMPORTANT:")
345: print(
346:     "This is an exploratory discovery model."
347: )
348: print(
349:     "It is not external validation."
350: )
351: print(
352:     "Feature selection was performed inside training folds."
353: )
~~~

## scripts/build_gse8841_borderline_dataset.py

~~~python
1: #!/usr/bin/env python3
2: 
3: from pathlib import Path
4: import csv
5: import gzip
6: import re
7: 
8: import numpy as np
9: import pandas as pd
10: 
11: 
12: ROOT = Path.home() / "ovarian-ai"
13: 
14: INPUT = (
15:     ROOT
16:     / "data"
17:     / "external"
18:     / "gse8841"
19:     / "GSE8841_series_matrix.txt.gz"
20: )
21: 
22: OUTPUT = (
23:     ROOT
24:     / "data"
25:     / "processed"
26:     / "gse8841_borderline_expression.csv.gz"
27: )
28: 
29: print("=" * 80)
30: print("BUILDING GSE8841 BORDERLINE EXPRESSION DATASET")
31: print("=" * 80)
32: 
33: # ------------------------------------------------------------------
34: # Read sample metadata without losing repeated characteristics rows
35: # ------------------------------------------------------------------
36: 
37: metadata = {}
38: 
39: with gzip.open(
40:     INPUT,
41:     "rt",
42:     encoding="utf-8",
43:     errors="replace",
44: ) as f:
45: 
46:     for line in f:
47: 
48:         if line.startswith(
49:             "!series_matrix_table_begin"
50:         ):
51:             break
52: 
53:         if line.startswith("!Sample_"):
54: 
55:             parts = next(
56:                 csv.reader(
57:                     [line.rstrip("\n")],
58:                     delimiter="\t",
59:                     quotechar='"',
60:                 )
61:             )
62: 
63:             metadata.setdefault(
64:                 parts[0],
65:                 []
66:             ).append(parts[1:])
67: 
68: 
69: accessions = metadata[
70:     "!Sample_geo_accession"
71: ][0]
72: 
73: titles = metadata[
74:     "!Sample_title"
75: ][0]
76: 
77: char_rows = metadata[
78:     "!Sample_characteristics_ch1"
79: ]
80: 
81: n_samples = len(accessions)
82: 
83: print()
84: print("Samples:", n_samples)
85: 
86: # ------------------------------------------------------------------
87: # Build per-sample metadata
88: # ------------------------------------------------------------------
89: 
90: sample_info = [
91:     {}
92:     for _ in range(n_samples)
93: ]
94: 
95: for row in char_rows:
96: 
97:     for i, value in enumerate(row):
98: 
99:         if ":" not in value:
100:             continue
101: 
102:         key, val = value.split(
103:             ":",
104:             1,
105:         )
106: 
107:         sample_info[i][
108:             key.strip().lower()
109:         ] = val.strip()
110: 
111: 
112: records = []
113: 
114: for i in range(n_samples):
115: 
116:     info = sample_info[i]
117: 
118:     stage = info.get(
119:         "figo stage",
120:         "",
121:     )
122: 
123:     substage = info.get(
124:         "figo substage",
125:         "",
126:     )
127: 
128:     grade = info.get(
129:         "grade",
130:         "",
131:     )
132: 
133:     age = info.get(
134:         "age(years)",
135:         "",
136:     )
137: 
138:     os = info.get(
139:         "overall survival(days)",
140:         "",
141:     )
142: 
143:     pfs = info.get(
144:         "progression free survival(days)",
145:         "",
146:     )
147: 
148:     histotype = info.get(
149:         "histotype",
150:         "",
151:     )
152: 
153:     relapse = info.get(
154:         "relapsed",
155:         "",
156:     )
157: 
158:     status = info.get(
159:         "status",
160:         "",
161:     )
162: 
163:     # --------------------------------------------------------------
164:     # Define endpoint.
165:     #
166:     # Borderline is explicitly defined by Grade:Borderline.
167:     # Cancer is FIGO Stage:1 with non-borderline grade.
168:     # --------------------------------------------------------------
169: 
170:     if grade.strip().lower() == "borderline":
171: 
172:         label = "borderline"
173: 
174:     elif stage.strip() == "1":
175: 
176:         label = "stage_I_cancer"
177: 
178:     else:
179: 
180:         raise RuntimeError(
181:             f"Unexpected sample classification: "
182:             f"{accessions[i]} | {info}"
183:         )
184: 
185:     records.append({
186:         "sample_id": accessions[i],
187:         "title": titles[i],
188:         "label": label,
189:         "figo_stage": stage,
190:         "figo_substage": substage,
191:         "histotype": histotype,
192:         "grade": grade,
193:         "age": pd.to_numeric(
194:             age,
195:             errors="coerce",
196:         ),
197:         "overall_survival_days": pd.to_numeric(
198:             os,
199:             errors="coerce",
200:         ),
201:         "progression_free_survival_days": pd.to_numeric(
202:             pfs,
203:             errors="coerce",
204:         ),
205:         "relapsed": relapse,
206:         "status": status,
207:     })
208: 
209: 
210: clinical = pd.DataFrame(records)
211: 
212: # ------------------------------------------------------------------
213: # Parse expression matrix
214: # ------------------------------------------------------------------
215: 
216: print()
217: print("=" * 80)
218: print("READING EXPRESSION MATRIX")
219: print("=" * 80)
220: 
221: matrix_lines = []
222: 
223: with gzip.open(
224:     INPUT,
225:     "rt",
226:     encoding="utf-8",
227:     errors="replace",
228: ) as f:
229: 
230:     inside = False
231: 
232:     for line in f:
233: 
234:         if line.startswith(
235:             "!series_matrix_table_begin"
236:         ):
237:             inside = True
238:             continue
239: 
240:         if line.startswith(
241:             "!series_matrix_table_end"
242:         ):
243:             break
244: 
245:         if inside and line.strip():
246:             matrix_lines.append(
247:                 line
248:             )
249: 
250: if not matrix_lines:
251:     raise RuntimeError(
252:         "Expression matrix was not found."
253:     )
254: 
255: from io import StringIO
256: 
257: expr = pd.read_csv(
258:     StringIO(
259:         "".join(matrix_lines)
260:     ),
261:     sep="\t",
262:     quotechar='"',
263: )
264: 
265: print(
266:     "Expression rows:",
267:     len(expr)
268: )
269: 
270: print(
271:     "Expression columns:",
272:     len(expr.columns)
273: )
274: 
275: # ------------------------------------------------------------------
276: # Verify sample columns
277: # ------------------------------------------------------------------
278: 
279: sample_cols = [
280:     c
281:     for c in expr.columns
282:     if c in set(accessions)
283: ]
284: 
285: if len(sample_cols) != n_samples:
286: 
287:     missing = sorted(
288:         set(accessions)
289:         - set(sample_cols)
290:     )
291: 
292:     raise RuntimeError(
293:         "Expression/sample mismatch. "
294:         f"Missing sample columns: {missing}"
295:     )
296: 
297: print(
298:     "Matched sample columns:",
299:     len(sample_cols)
300: )
301: 
302: # ------------------------------------------------------------------
303: # Convert expression values to numeric
304: # ------------------------------------------------------------------
305: 
306: feature_id_col = expr.columns[0]
307: 
308: expr[feature_id_col] = (
309:     expr[feature_id_col]
310:     .astype(str)
311: )
312: 
313: for c in sample_cols:
314: 
315:     expr[c] = pd.to_numeric(
316:         expr[c],
317:         errors="coerce",
318:     )
319: 
320: # ------------------------------------------------------------------
321: # Transpose: patients x expression features
322: # ------------------------------------------------------------------
323: 
324: expr_patient = expr[
325:     [feature_id_col] + sample_cols
326: ].set_index(
327:     feature_id_col
328: ).T.reset_index()
329: 
330: expr_patient = expr_patient.rename(
331:     columns={
332:         "index": "sample_id"
333:     }
334: )
335: 
336: expr_patient["sample_id"] = (
337:     expr_patient["sample_id"]
338:     .astype(str)
339: )
340: 
341: # ------------------------------------------------------------------
342: # Remove duplicate probe IDs if any
343: # ------------------------------------------------------------------
344: 
345: if expr_patient["sample_id"].duplicated().any():
346: 
347:     raise RuntimeError(
348:         "Duplicate sample IDs found."
349:     )
350: 
351: # ------------------------------------------------------------------
352: # Merge clinical + expression
353: # ------------------------------------------------------------------
354: 
355: df = clinical.merge(
356:     expr_patient,
357:     on="sample_id",
358:     how="inner",
359:     validate="one_to_one",
360: )
361: 
362: # ------------------------------------------------------------------
363: # Validate
364: # ------------------------------------------------------------------
365: 
366: expression_cols = [
367:     c
368:     for c in expr_patient.columns
369:     if c != "sample_id"
370: ]
371: 
372: print()
373: print("=" * 80)
374: print("FINAL DATASET VALIDATION")
375: print("=" * 80)
376: 
377: print(
378:     "Patients:",
379:     len(df)
380: )
381: 
382: print(
383:     "Stage-I cancer:",
384:     int(
385:         (
386:             df["label"]
387:             == "stage_I_cancer"
388:         ).sum()
389:     )
390: )
391: 
392: print(
393:     "Borderline:",
394:     int(
395:         (
396:             df["label"]
397:             == "borderline"
398:         ).sum()
399:     )
400: )
401: 
402: print(
403:     "Expression features:",
404:     len(expression_cols)
405: )
406: 
407: print(
408:     "Age observed:",
409:     int(
410:         df["age"].notna().sum()
411:     )
412: )
413: 
414: print(
415:     "Expression missing values:",
416:     int(
417:         df[expression_cols]
418:         .isna()
419:         .sum()
420:         .sum()
421:     )
422: )
423: 
424: print()
425: print("Histotype counts:")
426: print(
427:     pd.crosstab(
428:         df["label"],
429:         df["histotype"],
430:     ).to_string()
431: )
432: 
433: # ------------------------------------------------------------------
434: # Save
435: # ------------------------------------------------------------------
436: 
437: OUTPUT.parent.mkdir(
438:     parents=True,
439:     exist_ok=True,
440: )
441: 
442: df.to_csv(
443:     OUTPUT,
444:     index=False,
445:     compression="gzip",
446: )
447: 
448: print()
449: print("=" * 80)
450: print("SAVED")
451: print("=" * 80)
452: 
453: print(OUTPUT)
454: 
455: print()
456: print("IMPORTANT:")
457: print(
458:     "FIGO stage and grade are retained as metadata,"
459: )
460: print(
461:     "but must NOT be used as predictive features."
462: )
463: print(
464:     "This dataset is tissue-expression data,"
465: )
466: print(
467:     "not serum-miRNA data."
468: )
~~~

## scripts/build_gse8841_strict_borderline_model.py

~~~python
1: #!/usr/bin/env python3
2: 
3: from pathlib import Path
4: import json
5: 
6: import numpy as np
7: import pandas as pd
8: 
9: from sklearn.feature_selection import SelectKBest, f_classif
10: from sklearn.linear_model import LogisticRegression
11: from sklearn.metrics import (
12:     accuracy_score,
13:     balanced_accuracy_score,
14:     confusion_matrix,
15:     roc_auc_score,
16: )
17: from sklearn.model_selection import StratifiedKFold
18: from sklearn.pipeline import Pipeline
19: from sklearn.preprocessing import StandardScaler
20: 
21: 
22: ROOT = Path.home() / "ovarian-ai"
23: 
24: DATA = ROOT / "data/processed/gse8841_borderline_expression.csv.gz"
25: 
26: OUT_PRED = ROOT / "reports/tables/gse8841_strict_borderline_predictions.csv"
27: 
28: OUT_RESULTS = ROOT / "reports/tables/gse8841_strict_borderline_results.json"
29: 
30: 
31: print("=" * 80)
32: print("GSE8841 STRICT BORDERLINE MODEL")
33: print("=" * 80)
34: 
35: df = pd.read_csv(DATA)
36: 
37: # ---------------------------------------------------------------
38: # Restrict to serous + mucinous.
39: # Histotype is used ONLY for stratification, never as a feature.
40: # ---------------------------------------------------------------
41: 
42: df = df[
43:     df["histotype"].str.lower().isin(["serous", "mucinous"])
44: ].copy()
45: 
46: df["histotype"] = df["histotype"].str.lower()
47: 
48: y = (
49:     df["label"] == "stage_I_cancer"
50: ).astype(int)
51: 
52: # Create strata so each fold preserves BOTH diagnosis and histotype.
53: strata = (
54:     df["histotype"]
55:     + "_"
56:     + df["label"]
57: )
58: 
59: metadata = [
60:     "sample_id",
61:     "title",
62:     "label",
63:     "figo_stage",
64:     "figo_substage",
65:     "histotype",
66:     "grade",
67:     "age",
68:     "overall_survival_days",
69:     "progression_free_survival_days",
70:     "relapsed",
71:     "status",
72: ]
73: 
74: expression_cols = [
75:     c for c in df.columns
76:     if c not in metadata
77: ]
78: 
79: X = df[expression_cols].apply(
80:     pd.to_numeric,
81:     errors="coerce",
82: )
83: 
84: print()
85: print("Patients:", len(df))
86: print("Stage-I cancer:", int((y == 1).sum()))
87: print("Borderline:", int((y == 0).sum()))
88: print("Expression features:", X.shape[1])
89: 
90: print()
91: print("Histotype Ã— diagnosis:")
92: print(
93:     pd.crosstab(
94:         df["histotype"],
95:         df["label"],
96:     ).to_string()
97: )
98: 
99: # ---------------------------------------------------------------
100: # Remove expression features with missing values.
101: # ---------------------------------------------------------------
102: 
103: missing_fraction = X.isna().mean()
104: 
105: keep = missing_fraction == 0
106: 
107: X = X.loc[:, keep]
108: 
109: print()
110: print(
111:     "Features removed for missingness:",
112:     int((~keep).sum())
113: )
114: 
115: print(
116:     "Features retained:",
117:     X.shape[1]
118: )
119: 
120: if X.isna().any().any():
121:     raise RuntimeError("Missing values remain.")
122: 
123: # ---------------------------------------------------------------
124: # Pre-specified feature count.
125: #
126: # We are NOT searching for the best k after seeing the results.
127: # k=10 is fixed before evaluation.
128: # ---------------------------------------------------------------
129: 
130: K = min(10, X.shape[1])
131: 
132: N_REPEATS = 20
133: 
134: repeat_results = []
135: prediction_rows = []
136: 
137: for repeat in range(1, N_REPEATS + 1):
138: 
139:     cv = StratifiedKFold(
140:         n_splits=5,
141:         shuffle=True,
142:         random_state=1000 + repeat,
143:     )
144: 
145:     repeat_true = []
146:     repeat_prob = []
147:     repeat_pred = []
148: 
149:     for fold, (train_idx, test_idx) in enumerate(
150:         cv.split(X, strata),
151:         start=1,
152:     ):
153: 
154:         X_train = X.iloc[train_idx]
155:         X_test = X.iloc[test_idx]
156: 
157:         y_train = y.iloc[train_idx]
158:         y_test = y.iloc[test_idx]
159: 
160:         model = Pipeline([
161:             (
162:                 "scale",
163:                 StandardScaler(),
164:             ),
165:             (
166:                 "select",
167:                 SelectKBest(
168:                     score_func=f_classif,
169:                     k=K,
170:                 ),
171:             ),
172:             (
173:                 "clf",
174:                 LogisticRegression(
175:                     max_iter=5000,
176:                     class_weight="balanced",
177:                     C=0.1,
178:                     solver="liblinear",
179:                 ),
180:             ),
181:         ])
182: 
183:         model.fit(
184:             X_train,
185:             y_train,
186:         )
187: 
188:         prob = model.predict_proba(
189:             X_test
190:         )[:, 1]
191: 
192:         pred = (
193:             prob >= 0.5
194:         ).astype(int)
195: 
196:         repeat_true.extend(
197:             y_test.tolist()
198:         )
199: 
200:         repeat_prob.extend(
201:             prob.tolist()
202:         )
203: 
204:         repeat_pred.extend(
205:             pred.tolist()
206:         )
207: 
208:         for idx, p, pr in zip(
209:             test_idx,
210:             prob,
211:             pred,
212:         ):
213: 
214:             prediction_rows.append({
215:                 "repeat": repeat,
216:                 "fold": fold,
217:                 "sample_id": df.iloc[idx]["sample_id"],
218:                 "true_label": int(y.iloc[idx]),
219:                 "cancer_probability": float(p),
220:                 "predicted_label": int(pr),
221:             })
222: 
223:     repeat_true = np.asarray(
224:         repeat_true
225:     )
226: 
227:     repeat_prob = np.asarray(
228:         repeat_prob
229:     )
230: 
231:     repeat_pred = np.asarray(
232:         repeat_pred
233:     )
234: 
235:     auc = roc_auc_score(
236:         repeat_true,
237:         repeat_prob,
238:     )
239: 
240:     accuracy = accuracy_score(
241:         repeat_true,
242:         repeat_pred,
243:     )
244: 
245:     balanced_acc = balanced_accuracy_score(
246:         repeat_true,
247:         repeat_pred,
248:     )
249: 
250:     tn, fp, fn, tp = confusion_matrix(
251:         repeat_true,
252:         repeat_pred,
253:     ).ravel()
254: 
255:     sensitivity = (
256:         tp / (tp + fn)
257:         if tp + fn
258:         else np.nan
259:     )
260: 
261:     specificity = (
262:         tn / (tn + fp)
263:         if tn + fp
264:         else np.nan
265:     )
266: 
267:     repeat_results.append({
268:         "repeat": repeat,
269:         "auc": float(auc),
270:         "accuracy": float(accuracy),
271:         "balanced_accuracy": float(balanced_acc),
272:         "sensitivity": float(sensitivity),
273:         "specificity": float(specificity),
274:         "tn": int(tn),
275:         "fp": int(fp),
276:         "fn": int(fn),
277:         "tp": int(tp),
278:     })
279: 
280:     print(
281:         f"Repeat {repeat:02d}/20: "
282:         f"AUC={auc:.6f} "
283:         f"balanced_acc={balanced_acc:.6f}"
284:     )
285: 
286: # ---------------------------------------------------------------
287: # Summary across complete repeats
288: # ---------------------------------------------------------------
289: 
290: aucs = np.array([
291:     r["auc"]
292:     for r in repeat_results
293: ])
294: 
295: balanced_accs = np.array([
296:     r["balanced_accuracy"]
297:     for r in repeat_results
298: ])
299: 
300: sensitivities = np.array([
301:     r["sensitivity"]
302:     for r in repeat_results
303: ])
304: 
305: specificities = np.array([
306:     r["specificity"]
307:     for r in repeat_results
308: ])
309: 
310: print()
311: print("=" * 80)
312: print("FINAL GSE8841 STRICT BORDERLINE RESULTS")
313: print("=" * 80)
314: 
315: print(
316:     f"Patients: {len(df)}"
317: )
318: 
319: print(
320:     f"Stage-I cancer: {int((y == 1).sum())}"
321: )
322: 
323: print(
324:     f"Borderline: {int((y == 0).sum())}"
325: )
326: 
327: print()
328: print(
329:     f"Mean 5-fold AUC:       {aucs.mean():.6f}"
330: )
331: 
332: print(
333:     f"SD of AUC:             {aucs.std(ddof=1):.6f}"
334: )
335: 
336: print(
337:     f"Median AUC:            {np.median(aucs):.6f}"
338: )
339: 
340: print(
341:     f"Min AUC:               {aucs.min():.6f}"
342: )
343: 
344: print(
345:     f"Max AUC:               {aucs.max():.6f}"
346: )
347: 
348: print()
349: print(
350:     f"Mean balanced accuracy:{balanced_accs.mean():.6f}"
351: )
352: 
353: print(
354:     f"Mean sensitivity:      {sensitivities.mean():.6f}"
355: )
356: 
357: print(
358:     f"Mean specificity:      {specificities.mean():.6f}"
359: )
360: 
361: print()
362: print("=" * 80)
363: print("INTERPRETATION")
364: print("=" * 80)
365: 
366: if aucs.mean() < 0.60:
367:     print(
368:         "Weak evidence of reproducible discrimination."
369:     )
370: elif aucs.mean() < 0.70:
371:     print(
372:         "Modest exploratory discrimination."
373:     )
374: elif aucs.mean() < 0.80:
375:     print(
376:         "Potentially meaningful exploratory discrimination."
377:     )
378: else:
379:     print(
380:         "Strong exploratory discrimination; independent "
381:         "replication is still required."
382:     )
383: 
384: print()
385: print("=" * 80)
386: print("SAVING")
387: print("=" * 80)
388: 
389: pred_df = pd.DataFrame(
390:     prediction_rows
391: )
392: 
393: pred_df.to_csv(
394:     OUT_PRED,
395:     index=False,
396: )
397: 
398: results = {
399:     "dataset": "GSE8841",
400:     "subset": "serous + mucinous",
401:     "patients": int(len(df)),
402:     "stage_I_cancer": int((y == 1).sum()),
403:     "borderline": int((y == 0).sum()),
404:     "expression_features_input": int(len(expression_cols)),
405:     "expression_features_retained": int(X.shape[1]),
406:     "features_selected_per_fold": int(K),
407:     "repeats": int(N_REPEATS),
408:     "folds_per_repeat": 5,
409:     "mean_auc": float(aucs.mean()),
410:     "sd_auc": float(aucs.std(ddof=1)),
411:     "median_auc": float(np.median(aucs)),
412:     "min_auc": float(aucs.min()),
413:     "max_auc": float(aucs.max()),
414:     "mean_balanced_accuracy": float(
415:         balanced_accs.mean()
416:     ),
417:     "mean_sensitivity": float(
418:         sensitivities.mean()
419:     ),
420:     "mean_specificity": float(
421:         specificities.mean()
422:     ),
423:     "repeat_results": repeat_results,
424: }
425: 
426: OUT_RESULTS.write_text(
427:     json.dumps(
428:         results,
429:         indent=2,
430:     )
431: )
432: 
433: print()
434: print(OUT_PRED)
435: print(OUT_RESULTS)
436: 
437: print()
438: print("NO CLINICAL CLAIM.")
439: print(
440:     "This is exploratory tissue-expression research."
441: )
442: print(
443:     "GSE8841 is not serum miRNA validation."
444: )
~~~

## scripts/build_gse94533_clean_dataset.py

~~~python
1: #!/usr/bin/env python3
2: 
3: from __future__ import annotations
4: 
5: import gzip
6: import json
7: import re
8: from pathlib import Path
9: 
10: import numpy as np
11: import pandas as pd
12: 
13: 
14: # ================================================================
15: # PROJECT PATHS
16: # ================================================================
17: 
18: ROOT = Path(__file__).resolve().parents[1]
19: 
20: RAW_DIR = ROOT / "data" / "raw"
21: INTERIM_DIR = ROOT / "data" / "interim"
22: PROCESSED_DIR = ROOT / "data" / "processed"
23: REPORT_DIR = ROOT / "reports" / "tables"
24: 
25: METADATA_INPUT = INTERIM_DIR / "gse94533_geo_metadata_full.csv"
26: 
27: METADATA_OUTPUT = PROCESSED_DIR / "gse94533_metadata.csv"
28: MATRIX_OUTPUT = PROCESSED_DIR / "gse94533_mirna_tpm.csv.gz"
29: AUDIT_OUTPUT = REPORT_DIR / "gse94533_clean_dataset_audit.json"
30: 
31: 
32: # ================================================================
33: # SAMPLE-ID NORMALIZATION
34: # ================================================================
35: 
36: def normalize_sample_id(value: object) -> str:
37: 
38:     if pd.isna(value):
39:         return ""
40: 
41:     s = str(value).strip().lower()
42: 
43:     # Correct Excel-style numeric IDs.
44:     #
45:     # Examples:
46:     #   50057.0 -> 50057
47:     #   50217.0 -> 50217
48:     #
49:     if re.fullmatch(r"\d+\.0+", s):
50:         s = s.split(".", 1)[0]
51: 
52:     # Remove punctuation.
53:     #
54:     # E-1001 -> e1001
55:     # 50057  -> 50057
56:     #
57:     s = re.sub(r"[^a-z0-9]", "", s)
58: 
59:     return s
60: 
61: 
62: # ================================================================
63: # FIND INPUT WORKBOOKS
64: # ================================================================
65: 
66: def find_workbooks() -> list[Path]:
67: 
68:     files = sorted(
69:         RAW_DIR.glob("GSE94533_Processed_file_Cohort*.xlsx")
70:     )
71: 
72:     if len(files) != 2:
73: 
74:         raise RuntimeError(
75:             "Expected exactly two processed GSE94533 workbooks.\n"
76:             f"Found {len(files)}:\n"
77:             + "\n".join(str(x) for x in files)
78:         )
79: 
80:     return files
81: 
82: 
83: # ================================================================
84: # LOAD GEO METADATA
85: # ================================================================
86: 
87: def load_metadata() -> pd.DataFrame:
88: 
89:     if not METADATA_INPUT.exists():
90: 
91:         raise RuntimeError(
92:             f"Could not find GEO metadata at:\n{METADATA_INPUT}\n\n"
93:             "Run scripts/investigate_gse94533_labels.py first."
94:         )
95: 
96:     meta = pd.read_csv(
97:         METADATA_INPUT,
98:         dtype=str,
99:     )
100: 
101:     required = {
102:         "geo_accession",
103:         "title",
104:         "diagnosis",
105:     }
106: 
107:     missing = required - set(meta.columns)
108: 
109:     if missing:
110: 
111:         raise RuntimeError(
112:             f"GEO metadata is missing required columns: {missing}"
113:         )
114: 
115:     meta = meta.copy()
116: 
117:     meta["sample_id_normalized"] = meta["title"].map(
118:         normalize_sample_id
119:     )
120: 
121:     # ------------------------------------------------------------
122:     # Normalize diagnosis
123:     # ------------------------------------------------------------
124: 
125:     def normalize_diagnosis(value: object) -> str:
126: 
127:         if pd.isna(value):
128:             return ""
129: 
130:         s = str(value).strip().lower()
131: 
132:         if s == "cancer":
133:             return "Cancer"
134: 
135:         if s == "benign":
136:             return "Benign"
137: 
138:         if s == "borderline":
139:             return "Borderline"
140: 
141:         if s == "control":
142:             return "Control"
143: 
144:         return str(value).strip()
145: 
146:     meta["diagnosis"] = meta["diagnosis"].map(
147:         normalize_diagnosis
148:     )
149: 
150:     # ------------------------------------------------------------
151:     # Binary endpoint
152:     # ------------------------------------------------------------
153: 
154:     valid = {
155:         "Cancer",
156:         "Benign",
157:         "Borderline",
158:         "Control",
159:     }
160: 
161:     bad = sorted(
162:         set(meta["diagnosis"].dropna())
163:         - valid
164:     )
165: 
166:     if bad:
167: 
168:         raise RuntimeError(
169:             "Unexpected GEO diagnosis values:\n"
170:             + "\n".join(map(str, bad))
171:         )
172: 
173:     meta["binary_label"] = (
174:         meta["diagnosis"] == "Cancer"
175:     ).astype(int)
176: 
177:     # ------------------------------------------------------------
178:     # Explicit exclusion
179:     # ------------------------------------------------------------
180: 
181:     meta["excluded_from_analysis"] = (
182:         meta["sample_id_normalized"] == "e1056"
183:     )
184: 
185:     return meta
186: 
187: 
188: # ================================================================
189: # FIND SAMPLE COLUMNS IN EXPRESSION WORKBOOK
190: # ================================================================
191: 
192: 
193: def find_sample_columns(
194:     workbook: Path,
195:     valid_ids: set[str],
196: ) -> list[tuple[object, str]]:
197: 
198:     """
199:     Identify expression-matrix sample columns without changing the
200:     original pandas column type.
201: 
202:     This matters because Cohort 1 uses identifiers such as:
203: 
204:         e1036
205: 
206:     while Cohort 2 uses numeric Excel headers such as:
207: 
208:         50217.0
209: 
210:     We must preserve the original column object for dataframe indexing
211:     and only use the normalized string for matching to GEO metadata.
212: 
213:     Returns
214:     -------
215:     list of:
216:         (original_excel_column_object, normalized_sample_id)
217:     """
218: 
219:     header = pd.read_excel(
220:         workbook,
221:         sheet_name="mirna_tpm",
222:         header=2,
223:         nrows=0,
224:         engine="openpyxl",
225:     )
226: 
227:     results = []
228: 
229:     for column in header.columns:
230: 
231:         normalized = normalize_sample_id(column)
232: 
233:         if normalized in valid_ids:
234: 
235:             # IMPORTANT:
236:             # Keep `column` exactly as pandas read it.
237:             #
238:             # Do NOT use str(column) here.
239:             #
240:             results.append(
241:                 (column, normalized)
242:             )
243: 
244:     return results
245: 
246: 
247: def read_expression_workbook(
248:     workbook: Path,
249:     valid_ids: set[str],
250: ) -> tuple[pd.DataFrame, pd.DataFrame]:
251: 
252:     print()
253:     print("=" * 90)
254:     print(f"READING {workbook.name}")
255:     print("=" * 90)
256: 
257:     sample_columns = find_sample_columns(
258:         workbook,
259:         valid_ids,
260:     )
261: 
262:     if not sample_columns:
263: 
264:         raise RuntimeError(
265:             f"No GEO sample IDs matched {workbook.name}."
266:         )
267: 
268:     print()
269:     print(f"Matched sample columns: {len(sample_columns)}")
270: 
271:     # ------------------------------------------------------------
272:     # Read the actual matrix.
273:     # ------------------------------------------------------------
274: 
275:     df = pd.read_excel(
276:         workbook,
277:         sheet_name="mirna_tpm",
278:         header=2,
279:         engine="openpyxl",
280:     )
281: 
282:     columns = list(df.columns)
283: 
284:     feature_column = columns[0]
285: 
286:     print(
287:         f"Feature column: {feature_column!r}"
288:     )
289: 
290:     # ------------------------------------------------------------
291:     # Select sample columns using the ORIGINAL pandas objects.
292:     # ------------------------------------------------------------
293:     #
294:     # Cohort 2 contains headers such as 50217.0.
295:     #
296:     # We must select:
297:     #
298:     #     50217.0
299:     #
300:     # rather than:
301:     #
302:     #     "50217"
303:     #
304:     # because those are different pandas column labels.
305: 
306:     selected_columns = [
307:         original
308:         for original, _ in sample_columns
309:     ]
310: 
311:     expression = df[
312:         [feature_column] + selected_columns
313:     ].copy()
314: 
315:     expression = expression.rename(
316:         columns={
317:             feature_column: "feature_name"
318:         }
319:     )
320: 
321:     # ------------------------------------------------------------
322:     # Clean feature names.
323:     # ------------------------------------------------------------
324: 
325:     expression["feature_name"] = (
326:         expression["feature_name"]
327:         .astype("string")
328:         .str.strip()
329:     )
330: 
331:     expression = expression[
332:         expression["feature_name"].notna()
333:         & (expression["feature_name"] != "")
334:     ].copy()
335: 
336:     # ------------------------------------------------------------
337:     # Convert expression values to numeric.
338:     # ------------------------------------------------------------
339: 
340:     for column in selected_columns:
341: 
342:         expression[column] = pd.to_numeric(
343:             expression[column],
344:             errors="coerce",
345:         )
346: 
347:     expression[selected_columns] = (
348:         expression[selected_columns]
349:         .replace([float("inf"), float("-inf")], float("nan"))
350:         .fillna(0.0)
351:     )
352: 
353:     # ------------------------------------------------------------
354:     # Rename sample columns AFTER selecting them.
355:     # ------------------------------------------------------------
356: 
357:     rename_map = {}
358: 
359:     mapping_rows = []
360: 
361:     for original, normalized in sample_columns:
362: 
363:         rename_map[original] = normalized
364: 
365:         mapping_rows.append(
366:             {
367:                 "source_workbook": workbook.name,
368:                 "excel_column": str(original),
369:                 "sample_id_normalized": normalized,
370:             }
371:         )
372: 
373:     expression = expression.rename(
374:         columns=rename_map
375:     )
376: 
377:     normalized_sample_columns = [
378:         normalized
379:         for _, normalized in sample_columns
380:     ]
381: 
382:     # ------------------------------------------------------------
383:     # Aggregate duplicate miRNA feature names.
384:     # ------------------------------------------------------------
385: 
386:     grouped = (
387:         expression
388:         .groupby(
389:             "feature_name",
390:             sort=True,
391:         )[normalized_sample_columns]
392:         .sum()
393:     )
394: 
395:     print(
396:         f"Features after aggregation: {grouped.shape[0]}"
397:     )
398: 
399:     # ------------------------------------------------------------
400:     # Transpose:
401:     #
402:     # rows    = patients
403:     # columns = miRNAs
404:     # ------------------------------------------------------------
405: 
406:     patient_matrix = grouped.T.copy()
407: 
408:     patient_matrix.index = [
409:         normalize_sample_id(x)
410:         for x in patient_matrix.index
411:     ]
412: 
413:     patient_matrix.index.name = (
414:         "sample_id_normalized"
415:     )
416: 
417:     mapping = pd.DataFrame(
418:         mapping_rows
419:     )
420: 
421:     # ------------------------------------------------------------
422:     # Safety checks.
423:     # ------------------------------------------------------------
424: 
425:     if len(patient_matrix.index) != len(set(patient_matrix.index)):
426: 
427:         raise RuntimeError(
428:             f"Duplicate normalized patient IDs detected in "
429:             f"{workbook.name} after renaming."
430:         )
431: 
432:     return patient_matrix, mapping
433: 
434: def build_expression_matrix(
435:     metadata: pd.DataFrame,
436: ) -> tuple[pd.DataFrame, pd.DataFrame]:
437: 
438:     workbooks = find_workbooks()
439: 
440:     valid_geo_ids = set(
441:         metadata["sample_id_normalized"]
442:     )
443: 
444:     matrices = []
445:     mappings = []
446: 
447:     for workbook in workbooks:
448: 
449:         matrix, mapping = read_expression_workbook(
450:             workbook=workbook,
451:             valid_ids=valid_geo_ids,
452:         )
453: 
454:         matrices.append(matrix)
455:         mappings.append(mapping)
456: 
457:     # ------------------------------------------------------------
458:     # Make sure sample IDs aren't duplicated across workbooks.
459:     # ------------------------------------------------------------
460: 
461:     all_ids = [
462:         sample_id
463:         for matrix in matrices
464:         for sample_id in matrix.index
465:     ]
466: 
467:     duplicates = sorted(
468:         {
469:             x
470:             for x in all_ids
471:             if all_ids.count(x) > 1
472:         }
473:     )
474: 
475:     if duplicates:
476: 
477:         raise RuntimeError(
478:             "Duplicate patients across expression workbooks:\n"
479:             + "\n".join(duplicates)
480:         )
481: 
482:     # ------------------------------------------------------------
483:     # Combine all patients.
484:     # ------------------------------------------------------------
485: 
486:     combined = pd.concat(
487:         matrices,
488:         axis=0,
489:         join="outer",
490:     )
491: 
492:     combined = combined.fillna(0.0)
493: 
494:     # ------------------------------------------------------------
495:     # Remove E-1056.
496:     # ------------------------------------------------------------
497: 
498:     if "e1056" in combined.index:
499: 
500:         print()
501:         print(
502:             "Removing E-1056 according to the original study."
503:         )
504: 
505:         combined = combined.drop(
506:             index="e1056"
507:         )
508: 
509:     # ------------------------------------------------------------
510:     # Restrict GEO metadata to expression patients.
511:     # ------------------------------------------------------------
512: 
513:     metadata = metadata[
514:         ~metadata["excluded_from_analysis"]
515:     ].copy()
516: 
517:     metadata_ids = set(
518:         metadata["sample_id_normalized"]
519:     )
520: 
521:     expression_ids = set(
522:         combined.index
523:     )
524: 
525:     missing_expression = sorted(
526:         metadata_ids - expression_ids
527:     )
528: 
529:     missing_metadata = sorted(
530:         expression_ids - metadata_ids
531:     )
532: 
533:     if missing_expression:
534: 
535:         print()
536:         print(
537:             "WARNING: GEO metadata samples missing from expression:"
538:         )
539: 
540:         for sample_id in missing_expression:
541:             print(" ", sample_id)
542: 
543:     if missing_metadata:
544: 
545:         print()
546:         print(
547:             "WARNING: Expression samples missing from GEO metadata:"
548:         )
549: 
550:         for sample_id in missing_metadata:
551:             print(" ", sample_id)
552: 
553:     # Only retain patients for which we have both expression and labels.
554:     common = [
555:         sample_id
556:         for sample_id in metadata[
557:             "sample_id_normalized"
558:         ]
559:         if sample_id in combined.index
560:     ]
561: 
562:     metadata = (
563:         metadata
564:         .set_index("sample_id_normalized")
565:         .loc[common]
566:         .reset_index()
567:     )
568: 
569:     combined = combined.loc[
570:         common
571:     ]
572: 
573:     return combined, metadata
574: 
575: 
576: # ================================================================
577: # DATASET AUDIT
578: # ================================================================
579: 
580: def audit_dataset(
581:     matrix: pd.DataFrame,
582:     metadata: pd.DataFrame,
583: ) -> dict:
584: 
585:     print()
586:     print("=" * 90)
587:     print("FINAL DATASET AUDIT")
588:     print("=" * 90)
589: 
590:     n_patients = int(matrix.shape[0])
591:     n_features = int(matrix.shape[1])
592: 
593:     diagnosis_counts = (
594:         metadata["diagnosis"]
595:         .value_counts()
596:         .to_dict()
597:     )
598: 
599:     cancer_count = int(
600:         (metadata["binary_label"] == 1).sum()
601:     )
602: 
603:     noncancer_count = int(
604:         (metadata["binary_label"] == 0).sum()
605:     )
606: 
607:     finite = bool(
608:         np.isfinite(
609:             matrix.to_numpy(dtype=float)
610:         ).all()
611:     )
612: 
613:     nonnegative = bool(
614:         (
615:             matrix.to_numpy(dtype=float) >= 0
616:         ).all()
617:     )
618: 
619:     unique_patients = (
620:         not matrix.index.duplicated().any()
621:     )
622: 
623:     aligned = (
624:         list(matrix.index)
625:         == list(metadata["sample_id_normalized"])
626:     )
627: 
628:     audit = {
629:         "dataset": "GSE94533",
630:         "n_patients": n_patients,
631:         "n_mirna_features": n_features,
632:         "diagnosis_counts": diagnosis_counts,
633:         "binary_counts": {
634:             "cancer": cancer_count,
635:             "non_cancer": noncancer_count,
636:         },
637:         "excluded_samples": [
638:             "E-1056",
639:         ],
640:         "matrix_checks": {
641:             "all_finite": finite,
642:             "all_nonnegative": nonnegative,
643:             "unique_patient_ids": unique_patients,
644:             "metadata_expression_aligned": aligned,
645:         },
646:         "endpoint_definition": {
647:             "1": "GEO diagnosis = Cancer",
648:             "0": (
649:                 "GEO diagnosis = Benign, Borderline, or Control"
650:             ),
651:         },
652:         "important_provenance_note": (
653:             "Current GEO diagnosis annotations are preserved. "
654:             "They are not overwritten to force agreement with "
655:             "the diagnosis counts reported in the 2017 eLife paper."
656:         ),
657:     }
658: 
659:     print()
660:     print(f"Patients:          {n_patients}")
661:     print(f"miRNA features:    {n_features}")
662: 
663:     print()
664:     print("Diagnosis counts:")
665:     for key, value in diagnosis_counts.items():
666:         print(
667:             f"  {key:12s} {value}"
668:         )
669: 
670:     print()
671:     print("Binary counts:")
672:     print(
673:         f"  Cancer:        {cancer_count}"
674:     )
675:     print(
676:         f"  Non-cancer:    {noncancer_count}"
677:     )
678: 
679:     print()
680:     print("Matrix checks:")
681:     print(
682:         f"  Finite values:        {finite}"
683:     )
684:     print(
685:         f"  Non-negative values:   {nonnegative}"
686:     )
687:     print(
688:         f"  Unique patient IDs:    {unique_patients}"
689:     )
690:     print(
691:         f"  Metadata aligned:      {aligned}"
692:     )
693: 
694:     # ------------------------------------------------------------
695:     # Hard integrity requirements.
696:     # ------------------------------------------------------------
697: 
698:     if n_patients != 179:
699: 
700:         raise RuntimeError(
701:             f"Expected 179 eligible patients after E-1056 exclusion, "
702:             f"but found {n_patients}."
703:         )
704: 
705:     if not finite:
706: 
707:         raise RuntimeError(
708:             "Expression matrix contains NaN or infinite values."
709:         )
710: 
711:     if not nonnegative:
712: 
713:         raise RuntimeError(
714:             "Expression matrix contains negative TPM values."
715:         )
716: 
717:     if not unique_patients:
718: 
719:         raise RuntimeError(
720:             "Duplicate patient IDs remain."
721:         )
722: 
723:     if not aligned:
724: 
725:         raise RuntimeError(
726:             "Metadata and expression matrix are not aligned."
727:         )
728: 
729:     # We intentionally do NOT impose the historical 98/21/45/15
730:     # count constraint here, because the current GEO metadata are
731:     # the labels being used for this dataset.
732: 
733:     return audit
734: 
735: 
736: # ================================================================
737: # SAVE DATA
738: # ================================================================
739: 
740: def save_outputs(
741:     matrix: pd.DataFrame,
742:     metadata: pd.DataFrame,
743:     audit: dict,
744: ) -> None:
745: 
746:     INTERIM_DIR.mkdir(
747:         parents=True,
748:         exist_ok=True,
749:     )
750: 
751:     PROCESSED_DIR.mkdir(
752:         parents=True,
753:         exist_ok=True,
754:     )
755: 
756:     REPORT_DIR.mkdir(
757:         parents=True,
758:         exist_ok=True,
759:     )
760: 
761:     # ------------------------------------------------------------
762:     # Metadata
763:     # ------------------------------------------------------------
764: 
765:     metadata.to_csv(
766:         METADATA_OUTPUT,
767:         index=False,
768:     )
769: 
770:     # ------------------------------------------------------------
771:     # Matrix
772:     # ------------------------------------------------------------
773: 
774:     matrix.to_csv(
775:         MATRIX_OUTPUT,
776:         compression="gzip",
777:     )
778: 
779:     # ------------------------------------------------------------
780:     # Audit
781:     # ------------------------------------------------------------
782: 
783:     AUDIT_OUTPUT.write_text(
784:         json.dumps(
785:             audit,
786:             indent=2,
787:         )
788:     )
789: 
790:     print()
791:     print("=" * 90)
792:     print("OUTPUT FILES")
793:     print("=" * 90)
794: 
795:     print()
796:     print(
797:         f"Metadata:\n  {METADATA_OUTPUT}"
798:     )
799: 
800:     print()
801:     print(
802:         f"Expression matrix:\n  {MATRIX_OUTPUT}"
803:     )
804: 
805:     print()
806:     print(
807:         f"Audit:\n  {AUDIT_OUTPUT}"
808:     )
809: 
810: 
811: # ================================================================
812: # MAIN
813: # ================================================================
814: 
815: def main() -> None:
816: 
817:     print("=" * 90)
818:     print("GSE94533 CLEAN PATIENT-LEVEL DATASET BUILDER")
819:     print("=" * 90)
820: 
821:     # ------------------------------------------------------------
822:     # Load metadata.
823:     # ------------------------------------------------------------
824: 
825:     metadata = load_metadata()
826: 
827:     print()
828:     print(
829:         f"GEO metadata records: {len(metadata)}"
830:     )
831: 
832:     print()
833:     print("Current GEO diagnosis distribution:")
834: 
835:     print(
836:         metadata["diagnosis"]
837:         .value_counts()
838:         .to_string()
839:     )
840: 
841:     # ------------------------------------------------------------
842:     # Build expression matrix.
843:     # ------------------------------------------------------------
844: 
845:     matrix, metadata = build_expression_matrix(
846:         metadata
847:     )
848: 
849:     # ------------------------------------------------------------
850:     # Audit.
851:     # ------------------------------------------------------------
852: 
853:     audit = audit_dataset(
854:         matrix=matrix,
855:         metadata=metadata,
856:     )
857: 
858:     # ------------------------------------------------------------
859:     # Save.
860:     # ------------------------------------------------------------
861: 
862:     save_outputs(
863:         matrix=matrix,
864:         metadata=metadata,
865:         audit=audit,
866:     )
867: 
868:     print()
869:     print("=" * 90)
870:     print("SUCCESS")
871:     print("=" * 90)
872: 
873:     print(
874:         """
875: The clean patient-level GSE94533 dataset is now ready.
876: 
877: No machine-learning model has been trained.
878: 
879: The dataset now has:
880: 
881:     patient
882:        |
883:        +---- GEO diagnosis
884:        +---- binary cancer endpoint
885:        +---- stage
886:        +---- histology
887:        |
888:        +---- serum miRNA TPM features
889: 
890: Next we can build:
891: 
892:     1. A faithful 14-miRNA historical benchmark
893:     2. A modern regularized miRNA baseline
894:     3. Stage-I/II analysis
895:     4. Benign-vs-cancer analysis
896:     5. Later, the multimodal architecture
897:     """
898:     )
899: 
900: 
901: if __name__ == "__main__":
902:     main()
~~~

## scripts/build_gse94533_clinical_fusion.py

~~~python
1: #!/usr/bin/env python3
2: 
3: """Attach eLife Age and CA-125 values to the verified GSE94533 patients."""
4: 
5: import re
6: import zipfile
7: from pathlib import Path
8: from xml.etree import ElementTree as ET
9: 
10: import numpy as np
11: import pandas as pd
12: 
13: 
14: PROJECT = Path(__file__).resolve().parents[1]
15: MATRIX = PROJECT / "data" / "processed" / "gse94533_mirna_tpm.csv.gz"
16: METADATA = PROJECT / "data" / "processed" / "gse94533_metadata.csv"
17: SUPPLEMENT = PROJECT / "data" / "external" / "elife_28932" / "elife_poa_e28932_Supplementary_File_6.xlsx"
18: TABLE_DOCX = PROJECT / "data" / "external" / "elife_28932" / "elife_poa_e28932_Supplementary_File_1.docx"
19: OUTPUT = PROJECT / "data" / "processed" / "gse94533_14mirna_clinical.csv.gz"
20: 
21: FEATURES = [
22:     "hsa-miR-23b-3p", "hsa-miR-29a-3p", "hsa-miR-32-5p",
23:     "hsa-miR-92a-3p", "hsa-miR-150-5p", "hsa-miR-200a-3p",
24:     "hsa-miR-200c-3p", "hsa-miR-203a", "hsa-miR-320c",
25:     "hsa-miR-320d", "hsa-miR-335-5p", "hsa-miR-450b-5p",
26:     "hsa-miR-1246", "hsa-miR-1307-5p",
27: ]
28: 
29: 
30: def normalize(value):
31:     text = str(value).strip().lower()
32:     if re.fullmatch(r"\d+\.0+", text):
33:         text = text.split(".", 1)[0]
34:     return re.sub(r"[^a-z0-9]", "", text)
35: 
36: 
37: def cell_text(cell):
38:     return "".join(node.text or "" for node in cell.iter() if node.tag.endswith("}t")).strip()
39: 
40: 
41: def docx_tables(path):
42:     with zipfile.ZipFile(path) as archive:
43:         root = ET.fromstring(archive.read("word/document.xml"))
44:     tables = []
45:     for table in root.iter():
46:         if not table.tag.endswith("}tbl"):
47:             continue
48:         rows = []
49:         for tr in table:
50:             if tr.tag.endswith("}tr"):
51:                 rows.append([cell_text(tc) for tc in tr if tc.tag.endswith("}tc")])
52:         if rows:
53:             tables.append(rows)
54:     return tables
55: 
56: 
57: def number(value):
58:     text = str(value).strip()
59:     if text in {"", "*", "â€¡", "nan"}:
60:         return np.nan
61:     try:
62:         return float(text)
63:     except ValueError:
64:         return np.nan
65: 
66: 
67: def subject_by_number():
68:     sheet = pd.read_excel(SUPPLEMENT, sheet_name="Supplemental Dataset 6", header=None)
69:     labels = sheet.iloc[:, 0].astype(str).str.strip()
70:     subject_row = sheet.index[labels.eq("Subject ID")]
71:     number_row = sheet.index[labels.eq("Sample Number")]
72:     if len(subject_row) != 1 or len(number_row) != 1:
73:         raise RuntimeError("Could not identify Subject ID/Sample Number rows.")
74:     result = {}
75:     for sample_number, subject in zip(
76:         sheet.iloc[number_row[0], 1:], sheet.iloc[subject_row[0], 1:]
77:     ):
78:         if pd.notna(sample_number) and pd.notna(subject):
79:             result[int(float(sample_number))] = normalize(subject)
80:     return result
81: 
82: 
83: def clinical_by_sample():
84:     tables = docx_tables(TABLE_DOCX)
85:     table = next(
86:         table for table in tables
87:         if table and {str(x).strip() for x in table[0]} >= {"ID", "Age", "CA-125"}
88:     )
89:     header = table[0]
90:     id_col, age_col, ca_col = header.index("ID"), header.index("Age"), header.index("CA-125")
91:     number_map = subject_by_number()
92:     result = {}
93:     for row in table[1:]:
94:         if len(row) <= max(id_col, age_col, ca_col) or not row[id_col].strip():
95:             continue
96:         raw_id = row[id_col].strip()
97:         key = normalize(raw_id)
98:         if raw_id.replace(".", "", 1).isdigit():
99:             key = number_map.get(int(float(raw_id)), "")
100:         if key:
101:             result[key] = {"age": number(row[age_col]), "ca125": number(row[ca_col])}
102:     return result
103: 
104: 
105: def main():
106:     metadata = pd.read_csv(METADATA)
107:     matrix = pd.read_csv(MATRIX)
108:     missing = [feature for feature in FEATURES if feature not in matrix]
109:     if missing:
110:         raise RuntimeError(f"Missing historical miRNAs: {missing}")
111:     clinical = clinical_by_sample()
112:     base = metadata[["sample_id_normalized", "diagnosis", "binary_label"]].merge(
113:         matrix[["sample_id_normalized"] + FEATURES],
114:         on="sample_id_normalized", how="inner", validate="one_to_one",
115:     )
116:     clinical_frame = pd.DataFrame.from_dict(clinical, orient="index")
117:     clinical_frame.index.name = "sample_id_normalized"
118:     clinical_frame = clinical_frame.reset_index()
119:     output = base.merge(clinical_frame, on="sample_id_normalized", how="left", validate="one_to_one")
120:     if len(output) != 179:
121:         raise RuntimeError(f"Expected 179 patients, found {len(output)}")
122:     # CA-125 is legitimately missing for a subset of the published cohort;
123:     # the modeling pipeline imputes it inside each training fold.
124:     if output["age"].notna().sum() < 170 or output["ca125"].notna().sum() < 100:
125:         raise RuntimeError("Too few matched clinical values; refusing an uncertain join.")
126:     OUTPUT.parent.mkdir(parents=True, exist_ok=True)
127:     output.to_csv(OUTPUT, index=False, compression="gzip")
128:     print(f"Patients: {len(output)}")
129:     print(f"Age observed: {int(output['age'].notna().sum())}")
130:     print(f"CA-125 observed: {int(output['ca125'].notna().sum())}")
131:     print(f"Saved: {OUTPUT}")
132: 
133: 
134: if __name__ == "__main__":
135:     main()
~~~

## scripts/build_gse94533_dataset.py

~~~python
1: #!/usr/bin/env python3
2: 
3: """
4: Build a clean, patient-level GSE94533 serum-miRNA dataset.
5: 
6: Purpose
7: -------
8: GSE94533 contains processed serum miRNA expression matrices in two Excel
9: workbooks. The Excel workbooks contain sample IDs and molecular measurements,
10: but the diagnostic labels are stored in GEO sample metadata.
11: 
12: This script:
13: 
14: 1. Downloads the authoritative GEO family SOFT metadata.
15: 2. Extracts sample title, diagnosis, stage, histology, and GEO accession.
16: 3. Normalizes sample IDs so Excel and GEO identifiers can be joined safely.
17: 4. Reads the two processed GSE94533 miRNA TPM workbooks.
18: 5. Identifies the sample columns automatically.
19: 6. Excludes E-1056, which the original GSE94533 study excluded.
20: 7. Converts the GEO diagnosis into:
21:        1 = invasive ovarian cancer
22:        0 = non-cancer comparison group
23: 8. Preserves the original diagnosis as well.
24: 9. Builds a patient x miRNA matrix.
25: 10. Performs strict integrity checks.
26: 11. Writes:
27:        data/interim/gse94533_metadata.csv
28:        data/processed/gse94533_mirna_tpm.csv.gz
29:        data/processed/gse94533_dataset_manifest.json
30: 
31: IMPORTANT
32: ---------
33: This script deliberately DOES NOT train an AI model.
34: 
35: We first establish a correct ground-truth dataset.
36: """
37: 
38: from __future__ import annotations
39: 
40: import gzip
41: import json
42: import re
43: import sys
44: import urllib.request
45: from collections import Counter, defaultdict
46: from pathlib import Path
47: 
48: import numpy as np
49: import pandas as pd
50: 
51: 
52: # ---------------------------------------------------------------------
53: # Configuration
54: # ---------------------------------------------------------------------
55: 
56: ROOT = Path(__file__).resolve().parents[1]
57: 
58: RAW_DIR = ROOT / "data" / "raw"
59: INTERIM_DIR = ROOT / "data" / "interim"
60: PROCESSED_DIR = ROOT / "data" / "processed"
61: REPORT_DIR = ROOT / "reports" / "tables"
62: SCRIPT_DIR = ROOT / "scripts"
63: 
64: SOFT_URL = (
65:     "https://ftp.ncbi.nlm.nih.gov/geo/series/"
66:     "GSE94nnn/GSE94533/soft/GSE94533_family.soft.gz"
67: )
68: 
69: SOFT_PATH = INTERIM_DIR / "GSE94533_family.soft.gz"
70: 
71: METADATA_PATH = INTERIM_DIR / "gse94533_metadata.csv"
72: MATRIX_PATH = PROCESSED_DIR / "gse94533_mirna_tpm.csv.gz"
73: MANIFEST_PATH = PROCESSED_DIR / "gse94533_dataset_manifest.json"
74: 
75: 
76: # ---------------------------------------------------------------------
77: # Utility functions
78: # ---------------------------------------------------------------------
79: 
80: def normalize_sample_id(value: object) -> str:
81:     """
82:     Convert different representations of the same sample identifier
83:     into one canonical form.
84: 
85:     Examples
86:     --------
87:     E-1036   -> e1036
88:     e1036    -> e1036
89:     50217    -> 50217
90:     50217.0  -> 50217
91:     """
92: 
93:     if pd.isna(value):
94:         return ""
95: 
96:     s = str(value).strip().lower()
97: 
98:     # Excel often turns integer-looking IDs into floats.
99:     if re.fullmatch(r"\d+\.0+", s):
100:         s = s.split(".", 1)[0]
101: 
102:     # Remove spaces and punctuation.
103:     s = re.sub(r"[^a-z0-9]", "", s)
104: 
105:     return s
106: 
107: 
108: def classify_diagnosis(raw_value: str) -> tuple[str | None, int | None]:
109:     """
110:     Convert the original GEO diagnosis into the study endpoint.
111: 
112:     Returns
113:     -------
114:     diagnosis_group:
115:         cancer / borderline / benign / healthy / None
116: 
117:     binary_label:
118:         1 for invasive ovarian cancer
119:         0 for comparison groups
120:         None for unknown
121:     """
122: 
123:     if raw_value is None:
124:         return None, None
125: 
126:     s = str(raw_value).strip().lower()
127: 
128:     if not s:
129:         return None, None
130: 
131:     # IMPORTANT:
132:     # Check non-cancer categories first.
133:     # This prevents a future label such as
134:     # "borderline tumor" from being accidentally classified as cancer
135:     # just because it contains the word "tumor".
136: 
137:     if "borderline" in s:
138:         return "borderline", 0
139: 
140:     if "benign" in s:
141:         return "benign", 0
142: 
143:     if "healthy" in s:
144:         return "healthy", 0
145: 
146:     if "control" in s and "cancer" not in s:
147:         return "healthy", 0
148: 
149:     # Positive endpoint.
150:     if "invasive" in s:
151:         return "cancer", 1
152: 
153:     if "malignant" in s:
154:         return "cancer", 1
155: 
156:     if "cancer" in s:
157:         return "cancer", 1
158: 
159:     return None, None
160: 
161: 
162: def download_file(url: str, destination: Path) -> None:
163:     """
164:     Download a file from NCBI with a normal HTTP User-Agent.
165:     """
166: 
167:     destination.parent.mkdir(parents=True, exist_ok=True)
168: 
169:     if destination.exists() and destination.stat().st_size > 0:
170:         print(f"[download] Already exists: {destination}")
171:         return
172: 
173:     print(f"[download] Downloading:")
174:     print(f"           {url}")
175:     print(f"           -> {destination}")
176: 
177:     request = urllib.request.Request(
178:         url,
179:         headers={
180:             "User-Agent": "ovarian-ai-research-prototype/1.0"
181:         },
182:     )
183: 
184:     try:
185:         with urllib.request.urlopen(request, timeout=120) as response:
186:             data = response.read()
187:     except Exception as exc:
188:         raise RuntimeError(
189:             f"Could not download GEO metadata.\n"
190:             f"URL: {url}\n"
191:             f"Error: {exc}"
192:         ) from exc
193: 
194:     destination.write_bytes(data)
195: 
196:     print(
197:         f"[download] Complete: {destination.stat().st_size:,} bytes"
198:     )
199: 
200: 
201: # ---------------------------------------------------------------------
202: # GEO SOFT parser
203: # ---------------------------------------------------------------------
204: 
205: def parse_geo_soft(path: Path) -> pd.DataFrame:
206:     """
207:     Parse the GSE94533 family SOFT file.
208: 
209:     We only extract sample-level metadata.
210:     """
211: 
212:     print()
213:     print("=" * 90)
214:     print("PARSING GEO SAMPLE METADATA")
215:     print("=" * 90)
216: 
217:     records: list[dict] = []
218: 
219:     current: dict | None = None
220: 
221:     with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
222: 
223:         for raw_line in fh:
224: 
225:             line = raw_line.rstrip("\n")
226: 
227:             # A new GEO sample record.
228:             if line.startswith("^SAMPLE = "):
229: 
230:                 if current is not None:
231:                     records.append(current)
232: 
233:                 accession = line.split("=", 1)[1].strip()
234: 
235:                 current = {
236:                     "geo_accession": accession,
237:                     "title": None,
238:                     "characteristics": {},
239:                 }
240: 
241:                 continue
242: 
243:             if current is None:
244:                 continue
245: 
246:             if line.startswith("!Sample_title = "):
247: 
248:                 title = line.split("=", 1)[1].strip()
249:                 current["title"] = title
250:                 continue
251: 
252:             if line.startswith("!Sample_characteristics_ch1 = "):
253: 
254:                 value = line.split("=", 1)[1].strip()
255: 
256:                 if ":" in value:
257:                     key, val = value.split(":", 1)
258:                     key = key.strip().lower()
259:                     val = val.strip()
260:                 else:
261:                     key = "other"
262:                     val = value.strip()
263: 
264:                 current["characteristics"][key] = val
265: 
266:                 continue
267: 
268:     if current is not None:
269:         records.append(current)
270: 
271:     rows = []
272: 
273:     for record in records:
274: 
275:         chars = record["characteristics"]
276: 
277:         diagnosis_raw = chars.get("diagnosis")
278: 
279:         diagnosis_group, binary_label = classify_diagnosis(
280:             diagnosis_raw
281:         )
282: 
283:         rows.append(
284:             {
285:                 "geo_accession": record["geo_accession"],
286:                 "sample_id": record["title"],
287:                 "sample_id_normalized": normalize_sample_id(
288:                     record["title"]
289:                 ),
290:                 "diagnosis_raw": diagnosis_raw,
291:                 "diagnosis_group": diagnosis_group,
292:                 "binary_label": binary_label,
293:                 "stage": chars.get("stage"),
294:                 "histology": chars.get("histology"),
295:             }
296:         )
297: 
298:     metadata = pd.DataFrame(rows)
299: 
300:     print()
301:     print(f"GEO sample records found: {len(metadata)}")
302: 
303:     if len(metadata) != 180:
304:         print(
305:             "WARNING: expected 180 GEO sample records according to "
306:             "the GSE94533 record."
307:         )
308: 
309:     return metadata
310: 
311: 
312: # ---------------------------------------------------------------------
313: # Metadata validation
314: # ---------------------------------------------------------------------
315: 
316: def validate_geo_metadata(metadata: pd.DataFrame) -> None:
317: 
318:     print()
319:     print("=" * 90)
320:     print("VALIDATING GEO METADATA")
321:     print("=" * 90)
322: 
323:     print()
324:     print("Diagnosis raw values:")
325:     print(
326:         metadata["diagnosis_raw"]
327:         .fillna("<MISSING>")
328:         .value_counts(dropna=False)
329:         .to_string()
330:     )
331: 
332:     print()
333:     print("Normalized diagnosis groups:")
334: 
335:     print(
336:         metadata["diagnosis_group"]
337:         .fillna("<UNKNOWN>")
338:         .value_counts(dropna=False)
339:         .to_string()
340:     )
341: 
342:     unknown = metadata[metadata["diagnosis_group"].isna()]
343: 
344:     if len(unknown) > 0:
345: 
346:         print()
347:         print("WARNING: Unknown diagnosis records:")
348:         print(
349:             unknown[
350:                 [
351:                     "geo_accession",
352:                     "sample_id",
353:                     "diagnosis_raw",
354:                 ]
355:             ].to_string(index=False)
356:         )
357: 
358:         raise RuntimeError(
359:             "GEO diagnosis values could not all be mapped safely. "
360:             "The script refuses to guess their labels."
361:         )
362: 
363:     # E-1056 is intentionally still present at this point. The original
364:     # publication's table counts do not exactly match the current GEO
365:     # annotation, so the source metadata is retained as authoritative and
366:     # the excluded sample is recorded explicitly below.
367:     print()
368:     print("Expected post-exclusion sample count: 179")
369:     print(f"Current GEO metadata count: {len(metadata)}")
370: 
371:     # Verify the published diagnosis composition among the 179
372:     # eligible samples after excluding E-1056.
373:     eligible = metadata[
374:         metadata["sample_id_normalized"] != "e1056"
375:     ].copy()
376: 
377:     observed = (
378:         eligible["diagnosis_group"]
379:         .value_counts()
380:         .to_dict()
381:     )
382: 
383:     print()
384:     print("Observed eligible diagnosis counts:")
385:     print(json.dumps(observed, indent=2))
386: 
387:     if len(eligible) != 179:
388:         raise RuntimeError(
389:             f"Eligible patient count is {len(eligible)}, "
390:             "but expected 179 after excluding E-1056."
391:         )
392: 
393:     # Check sample ID uniqueness.
394:     duplicate_ids = (
395:         eligible["sample_id_normalized"]
396:         .loc[
397:             eligible["sample_id_normalized"].duplicated(
398:                 keep=False
399:             )
400:         ]
401:         .unique()
402:         .tolist()
403:     )
404: 
405:     if duplicate_ids:
406: 
407:         raise RuntimeError(
408:             "Duplicate normalized sample IDs found: "
409:             + ", ".join(duplicate_ids)
410:         )
411: 
412:     print()
413:     print("GEO metadata validation PASSED.")
414: 
415: 
416: # ---------------------------------------------------------------------
417: # Excel loading
418: # ---------------------------------------------------------------------
419: 
420: def find_workbooks() -> list[Path]:
421: 
422:     candidates = sorted(
423:         RAW_DIR.glob("GSE94533_Processed_file_Cohort*.xlsx")
424:     )
425: 
426:     if len(candidates) != 2:
427: 
428:         raise RuntimeError(
429:             "Expected exactly two GSE94533 processed Excel workbooks, "
430:             f"but found {len(candidates)}:\n"
431:             + "\n".join(str(p) for p in candidates)
432:         )
433: 
434:     return candidates
435: 
436: 
437: def load_mirna_matrix_from_workbook(
438:     workbook: Path,
439:     metadata_ids: set[str],
440: ) -> tuple[pd.DataFrame, pd.DataFrame]:
441: 
442:     print()
443:     print("=" * 90)
444:     print(f"READING: {workbook.name}")
445:     print("=" * 90)
446: 
447:     print()
448:     print("Reading sheet 'mirna_tpm' with header row 3...")
449: 
450:     df = pd.read_excel(
451:         workbook,
452:         sheet_name="mirna_tpm",
453:         header=2,
454:         engine="openpyxl",
455:     )
456: 
457:     print(f"Raw shape: {df.shape}")
458: 
459:     if df.shape[1] < 5:
460:         raise RuntimeError(
461:             f"{workbook.name} does not look like the expected "
462:             "GSE94533 miRNA TPM matrix."
463:         )
464: 
465:     # Identify the name columns.
466:     columns = list(df.columns)
467: 
468:     feature_name_col = columns[0]
469: 
470:     # The standard processed file contains:
471:     # Name, Sequence, Name count, Identical miRNAs, samples...
472:     identical_col = None
473: 
474:     for col in columns[:6]:
475: 
476:         if str(col).strip().lower() in {
477:             "identical mirnas",
478:             "identical miRNAs".lower(),
479:         }:
480:             identical_col = col
481:             break
482: 
483:     if identical_col is None:
484:         # Fall back to the fourth column used by the processed file.
485:         if len(columns) >= 4:
486:             identical_col = columns[3]
487: 
488:     print()
489:     print(f"Feature-name column: {feature_name_col!r}")
490:     print(f"Identical-miRNA column: {identical_col!r}")
491: 
492:     sample_columns = []
493:     sample_map_rows = []
494: 
495:     for col in columns:
496: 
497:         normalized = normalize_sample_id(col)
498: 
499:         if normalized in metadata_ids:
500: 
501:             sample_columns.append(col)
502: 
503:             sample_map_rows.append(
504:                 {
505:                     "excel_column": str(col),
506:                     "sample_id_normalized": normalized,
507:                     "source_workbook": workbook.name,
508:                 }
509:             )
510: 
511:     print()
512:     print(f"Matched sample columns: {len(sample_columns)}")
513: 
514:     if not sample_columns:
515: 
516:         raise RuntimeError(
517:             f"No sample columns from {workbook.name} matched GEO metadata."
518:         )
519: 
520:     # Show the first several matches for human inspection.
521:     print()
522:     print("First sample-column matches:")
523: 
524:     for row in sample_map_rows[:20]:
525:         print(
526:             f"  {row['excel_column']!r}"
527:             f" -> {row['sample_id_normalized']}"
528:         )
529: 
530:     # Extract only the feature identifier and sample measurements.
531:     expression = df[
532:         [feature_name_col, identical_col] + sample_columns
533:     ].copy()
534: 
535:     expression = expression.rename(
536:         columns={
537:             feature_name_col: "feature_name",
538:             identical_col: "canonical_mirna",
539:         }
540:     )
541: 
542:     # Use the canonical miRNA name when available.
543:     expression["feature_id"] = (
544:         expression["canonical_mirna"]
545:         .astype("string")
546:         .fillna(
547:             expression["feature_name"].astype("string")
548:         )
549:         .str.strip()
550:     )
551: 
552:     # Remove empty feature identifiers.
553:     expression = expression[
554:         expression["feature_id"].notna()
555:         & (expression["feature_id"] != "")
556:     ].copy()
557: 
558:     # Convert expression columns to numeric.
559:     for col in sample_columns:
560:         expression[col] = pd.to_numeric(
561:             expression[col],
562:             errors="coerce",
563:         )
564: 
565:     # Replace missing expression measurements with zero only
566:     # after verifying that the measurement table is numeric.
567:     expression[sample_columns] = (
568:         expression[sample_columns]
569:         .replace([np.inf, -np.inf], np.nan)
570:         .fillna(0.0)
571:     )
572: 
573:     # Collapse identical canonical miRNAs.
574:     #
575:     # Multiple sequence representations can sometimes correspond to
576:     # the same canonical miRNA. For TPM, summing these representations
577:     # is more appropriate than silently dropping one.
578:     grouped = (
579:         expression
580:         .groupby("feature_id", sort=True)[sample_columns]
581:         .sum()
582:     )
583: 
584:     # Transpose so rows = patients and columns = miRNAs.
585:     patient_matrix = grouped.T.copy()
586: 
587:     # Normalize sample names.
588:     patient_matrix.index = [
589:         normalize_sample_id(x)
590:         for x in patient_matrix.index
591:     ]
592: 
593:     patient_matrix.index.name = "sample_id_normalized"
594: 
595:     sample_map = pd.DataFrame(sample_map_rows)
596: 
597:     return patient_matrix, sample_map
598: 
599: 
600: # ---------------------------------------------------------------------
601: # Combine both cohorts
602: # ---------------------------------------------------------------------
603: 
604: def build_patient_level_dataset(
605:     metadata: pd.DataFrame,
606: ) -> tuple[pd.DataFrame, pd.DataFrame]:
607: 
608:     workbooks = find_workbooks()
609: 
610:     metadata_ids = set(
611:         metadata["sample_id_normalized"].astype(str)
612:     )
613: 
614:     matrices = []
615:     mappings = []
616: 
617:     for workbook in workbooks:
618: 
619:         matrix, sample_map = load_mirna_matrix_from_workbook(
620:             workbook=workbook,
621:             metadata_ids=metadata_ids,
622:         )
623: 
624:         matrices.append(matrix)
625:         mappings.append(sample_map)
626: 
627:     # Check whether the two processed workbooks accidentally contain
628:     # overlapping patient IDs.
629:     all_indices = [
630:         sample_id
631:         for matrix in matrices
632:         for sample_id in matrix.index
633:     ]
634: 
635:     duplicates = [
636:         sample_id
637:         for sample_id, count in Counter(all_indices).items()
638:         if count > 1
639:     ]
640: 
641:     if duplicates:
642: 
643:         raise RuntimeError(
644:             "The processed workbooks share sample IDs:\n"
645:             + "\n".join(duplicates)
646:         )
647: 
648:     # Align features across the two workbooks.
649:     combined_matrix = pd.concat(
650:         matrices,
651:         axis=0,
652:         join="outer",
653:     )
654: 
655:     combined_matrix = combined_matrix.fillna(0.0)
656: 
657:     # Restrict metadata to samples actually present in expression data.
658:     eligible_meta = metadata[
659:         metadata["sample_id_normalized"].isin(
660:             combined_matrix.index
661:         )
662:         & (
663:             metadata["sample_id_normalized"] != "e1056"
664:         )
665:     ].copy()
666: 
667:     # Exclude E-1056 explicitly.
668:     combined_matrix = combined_matrix[
669:         combined_matrix.index != "e1056"
670:     ].copy()
671: 
672:     # Verify one row per patient.
673:     if combined_matrix.index.duplicated().any():
674:         raise RuntimeError(
675:             "Duplicate patient rows remain after combining workbooks."
676:         )
677: 
678:     # Determine sample source workbook.
679:     mapping = pd.concat(
680:         mappings,
681:         ignore_index=True,
682:     )
683: 
684:     source_lookup = dict(
685:         zip(
686:             mapping["sample_id_normalized"],
687:             mapping["source_workbook"],
688:         )
689:     )
690: 
691:     eligible_meta["source_workbook"] = (
692:         eligible_meta["sample_id_normalized"]
693:         .map(source_lookup)
694:     )
695: 
696:     # Keep metadata and expression in exactly the same order.
697:     common_ids = [
698:         sample_id
699:         for sample_id in eligible_meta["sample_id_normalized"]
700:         if sample_id in combined_matrix.index
701:     ]
702: 
703:     eligible_meta = (
704:         eligible_meta
705:         .set_index("sample_id_normalized")
706:         .loc[common_ids]
707:         .reset_index()
708:     )
709: 
710:     combined_matrix = combined_matrix.loc[
711:         common_ids
712:     ]
713: 
714:     return eligible_meta, combined_matrix
715: 
716: 
717: # ---------------------------------------------------------------------
718: # Final integrity checks
719: # ---------------------------------------------------------------------
720: 
721: def final_integrity_checks(
722:     metadata: pd.DataFrame,
723:     matrix: pd.DataFrame,
724: ) -> None:
725: 
726:     print()
727:     print("=" * 90)
728:     print("FINAL DATASET INTEGRITY CHECK")
729:     print("=" * 90)
730: 
731:     print()
732:     print(f"Patients: {matrix.shape[0]}")
733:     print(f"miRNA features: {matrix.shape[1]}")
734: 
735:     if matrix.shape[0] != 179:
736:         raise RuntimeError(
737:             f"Expected 179 eligible patients, found {matrix.shape[0]}."
738:         )
739: 
740:     if len(metadata) != matrix.shape[0]:
741:         raise RuntimeError(
742:             "Metadata and expression matrix patient counts do not match."
743:         )
744: 
745:     if list(metadata["sample_id_normalized"]) != list(matrix.index):
746:         raise RuntimeError(
747:             "Metadata patient order and matrix patient order do not match."
748:         )
749: 
750:     if matrix.index.duplicated().any():
751:         raise RuntimeError(
752:             "Duplicate patient IDs detected in expression matrix."
753:         )
754: 
755:     if not np.isfinite(matrix.to_numpy(dtype=float)).all():
756:         raise RuntimeError(
757:             "Expression matrix contains NaN or infinite values."
758:         )
759: 
760:     # Expression cannot be negative.
761:     if (matrix.to_numpy(dtype=float) < 0).any():
762:         raise RuntimeError(
763:             "Negative TPM values detected."
764:         )
765: 
766:     print()
767:     print("Diagnosis distribution:")
768:     print(
769:         metadata["diagnosis_group"]
770:         .value_counts()
771:         .reindex(
772:             ["cancer", "borderline", "benign", "healthy"]
773:         )
774:         .fillna(0)
775:         .astype(int)
776:         .to_string()
777:     )
778: 
779:     print()
780:     print("Binary endpoint distribution:")
781:     print(
782:         metadata["binary_label"]
783:         .value_counts()
784:         .sort_index()
785:         .to_string()
786:     )
787: 
788:     print()
789:     print("Patients by source workbook:")
790:     print(
791:         metadata["source_workbook"]
792:         .value_counts()
793:         .to_string()
794:     )
795: 
796:     print()
797:     print("Patients by diagnosis and source workbook:")
798:     print(
799:         pd.crosstab(
800:             metadata["source_workbook"],
801:             metadata["diagnosis_group"],
802:         )
803:         .reindex(
804:             columns=[
805:                 "cancer",
806:                 "borderline",
807:                 "benign",
808:                 "healthy",
809:             ],
810:             fill_value=0,
811:         )
812:         .to_string()
813:     )
814: 
815:     print()
816:     print("Stage distribution among cancer cases:")
817: 
818:     cancer_stage = metadata.loc[
819:         metadata["diagnosis_group"] == "cancer",
820:         "stage",
821:     ]
822: 
823:     print(
824:         cancer_stage
825:         .fillna("<MISSING>")
826:         .value_counts(dropna=False)
827:         .to_string()
828:     )
829: 
830:     print()
831:     print("Histology distribution:")
832:     print(
833:         metadata["histology"]
834:         .fillna("<MISSING>")
835:         .value_counts(dropna=False)
836:         .head(30)
837:         .to_string()
838:     )
839: 
840:     print()
841:     print("FINAL INTEGRITY CHECK PASSED.")
842: 
843: 
844: # ---------------------------------------------------------------------
845: # Write outputs
846: # ---------------------------------------------------------------------
847: 
848: def write_outputs(
849:     metadata: pd.DataFrame,
850:     matrix: pd.DataFrame,
851: ) -> None:
852: 
853:     INTERIM_DIR.mkdir(parents=True, exist_ok=True)
854:     PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
855:     REPORT_DIR.mkdir(parents=True, exist_ok=True)
856: 
857:     metadata.to_csv(
858:         METADATA_PATH,
859:         index=False,
860:     )
861: 
862:     matrix.to_csv(
863:         MATRIX_PATH,
864:         compression="gzip",
865:     )
866: 
867:     manifest = {
868:         "dataset": "GSE94533",
869:         "description": (
870:             "Patient-level serum miRNA TPM matrix aligned to GEO "
871:             "diagnosis metadata."
872:         ),
873:         "source": (
874:             "NCBI Gene Expression Omnibus GSE94533"
875:         ),
876:         "geo_metadata_source": SOFT_URL,
877:         "expression_sources": [
878:             str(p)
879:             for p in find_workbooks()
880:         ],
881:         "excluded_samples": [
882:             {
883:                 "sample_id": "E-1056",
884:                 "reason": (
885:                     "Excluded by original study because of "
886:                     "abnormally high miR-122 believed to be "
887:                     "associated with a recent myocardial infarction."
888:                 ),
889:             }
890:         ],
891:         "n_patients": int(matrix.shape[0]),
892:         "n_features": int(matrix.shape[1]),
893:         "diagnosis_counts": (
894:             metadata["diagnosis_group"]
895:             .value_counts()
896:             .to_dict()
897:         ),
898:         "binary_label_definition": {
899:             "1": "invasive ovarian cancer",
900:             "0": (
901:                 "borderline tumor, benign ovarian disease, "
902:                 "or healthy control"
903:             ),
904:         },
905:         "expression_measurement": "TPM",
906:         "aggregation": (
907:             "Multiple sequence representations mapping to the same "
908:             "canonical miRNA were summed."
909:         ),
910:         "note": (
911:             "Research dataset. This is not a validated clinical "
912:             "diagnostic test."
913:         ),
914:     }
915: 
916:     MANIFEST_PATH.write_text(
917:         json.dumps(
918:             manifest,
919:             indent=2,
920:         )
921:     )
922: 
923:     print()
924:     print("=" * 90)
925:     print("OUTPUTS WRITTEN")
926:     print("=" * 90)
927: 
928:     print()
929:     print(f"Metadata:")
930:     print(f"  {METADATA_PATH}")
931: 
932:     print()
933:     print(f"Expression matrix:")
934:     print(f"  {MATRIX_PATH}")
935: 
936:     print()
937:     print(f"Manifest:")
938:     print(f"  {MANIFEST_PATH}")
939: 
940: 
941: # ---------------------------------------------------------------------
942: # Main
943: # ---------------------------------------------------------------------
944: 
945: def main() -> int:
946: 
947:     print("=" * 90)
948:     print("GSE94533 PATIENT-LEVEL DATASET BUILDER")
949:     print("=" * 90)
950: 
951:     print()
952:     print(f"Project root: {ROOT}")
953: 
954:     RAW_DIR.mkdir(parents=True, exist_ok=True)
955:     INTERIM_DIR.mkdir(parents=True, exist_ok=True)
956:     PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
957:     REPORT_DIR.mkdir(parents=True, exist_ok=True)
958: 
959:     # -------------------------------------------------------------
960:     # Step 1: obtain official GEO metadata
961:     # -------------------------------------------------------------
962: 
963:     download_file(
964:         SOFT_URL,
965:         SOFT_PATH,
966:     )
967: 
968:     metadata = parse_geo_soft(SOFT_PATH)
969: 
970:     if metadata.empty:
971:         raise RuntimeError(
972:             "No GEO sample metadata were extracted."
973:         )
974: 
975:     validate_geo_metadata(metadata)
976: 
977:     # -------------------------------------------------------------
978:     # Step 2: build patient-level expression matrix
979:     # -------------------------------------------------------------
980: 
981:     final_metadata, matrix = build_patient_level_dataset(
982:         metadata
983:     )
984: 
985:     # -------------------------------------------------------------
986:     # Step 3: strict integrity checks
987:     # -------------------------------------------------------------
988: 
989:     final_integrity_checks(
990:         final_metadata,
991:         matrix,
992:     )
993: 
994:     # -------------------------------------------------------------
995:     # Step 4: write clean data products
996:     # -------------------------------------------------------------
997: 
998:     write_outputs(
999:         final_metadata,
1000:         matrix,
1001:     )
1002: 
1003:     print()
1004:     print("=" * 90)
1005:     print("GSE94533 DATASET BUILD COMPLETE")
1006:     print("=" * 90)
1007: 
1008:     print()
1009:     print("We have NOT trained a model yet.")
1010:     print()
1011:     print(
1012:         "The next stage will use this frozen patient-level dataset "
1013:         "to build the first miRNA diagnostic baseline."
1014:     )
1015: 
1016:     return 0
1017: 
1018: 
1019: if __name__ == "__main__":
1020:     try:
1021:         raise SystemExit(main())
1022:     except KeyboardInterrupt:
1023:         print("\nInterrupted by user.")
1024:         raise SystemExit(130)
1025:     except Exception as exc:
1026:         print()
1027:         print("=" * 90)
1028:         print("DATASET BUILD FAILED")
1029:         print("=" * 90)
1030:         print()
1031:         print(f"{type(exc).__name__}: {exc}")
1032:         print()
1033:         print(
1034:             "No model training was performed."
1035:         )
1036:         raise SystemExit(1)
~~~

## scripts/build_modern_14mirna_baseline.py

~~~python
1: #!/usr/bin/env python3
2: 
3: """
4: MODERN 14-miRNA BASELINE
5: ========================
6: 
7: Build a leakage-safe modern baseline using the exact 14 miRNAs from
8: the historical eLife PMML model.
9: 
10: This script DOES NOT use the historical neural-network weights.
11: 
12: It trains a new regularized logistic-regression model using:
13: 
14:     exact 14 historical miRNAs
15:         ->
16:     StandardScaler
17:         ->
18:     LogisticRegression
19: 
20: Evaluation:
21: 
22:     5-fold stratified patient-level cross-validation
23: 
24: All preprocessing occurs INSIDE each fold.
25: 
26: The output is an out-of-fold prediction for every patient.
27: 
28: This provides a fair modern baseline against which the fixed
29: historical PMML network can later be compared.
30: 
31: No neural network is trained here.
32: """
33: 
34: from __future__ import annotations
35: 
36: import json
37: from pathlib import Path
38: 
39: import numpy as np
40: import pandas as pd
41: 
42: from sklearn.linear_model import LogisticRegression
43: from sklearn.metrics import (
44:     accuracy_score,
45:     average_precision_score,
46:     confusion_matrix,
47:     roc_auc_score,
48: )
49: from sklearn.model_selection import StratifiedKFold
50: from sklearn.pipeline import Pipeline
51: from sklearn.preprocessing import StandardScaler
52: 
53: 
54: # =====================================================================
55: # PATHS
56: # =====================================================================
57: 
58: ROOT = Path.home() / "ovarian-ai"
59: 
60: DATA_PATH = (
61:     ROOT
62:     / "data"
63:     / "processed"
64:     / "gse94533_14mirna_benchmark.csv.gz"
65: )
66: 
67: OUTPUT_DIR = (
68:     ROOT
69:     / "reports"
70:     / "tables"
71: )
72: 
73: PREDICTIONS_OUTPUT = (
74:     OUTPUT_DIR
75:     / "modern_14mirna_oof_predictions.csv"
76: )
77: 
78: RESULTS_OUTPUT = (
79:     OUTPUT_DIR
80:     / "modern_14mirna_baseline_results.json"
81: )
82: 
83: 
84: # =====================================================================
85: # EXACT HISTORICAL FEATURES
86: # =====================================================================
87: 
88: FEATURES = [
89:     "hsa-miR-23b-3p",
90:     "hsa-miR-29a-3p",
91:     "hsa-miR-32-5p",
92:     "hsa-miR-92a-3p",
93:     "hsa-miR-150-5p",
94:     "hsa-miR-200a-3p",
95:     "hsa-miR-200c-3p",
96:     "hsa-miR-203a",
97:     "hsa-miR-320c",
98:     "hsa-miR-320d",
99:     "hsa-miR-335-5p",
100:     "hsa-miR-450b-5p",
101:     "hsa-miR-1246",
102:     "hsa-miR-1307-5p",
103: ]
104: 
105: 
106: # =====================================================================
107: # LOAD DATA
108: # =====================================================================
109: 
110: print("=" * 90)
111: print("MODERN 14-miRNA BASELINE")
112: print("=" * 90)
113: 
114: print()
115: print("Loading:")
116: print(DATA_PATH)
117: 
118: if not DATA_PATH.exists():
119:     raise RuntimeError(
120:         f"Dataset not found:\n{DATA_PATH}"
121:     )
122: 
123: df = pd.read_csv(DATA_PATH)
124: 
125: print()
126: print("Dataset shape:", df.shape)
127: 
128: print()
129: print("Columns:")
130: print(df.columns.tolist())
131: 
132: 
133: # =====================================================================
134: # IDENTIFY LABEL COLUMN
135: # =====================================================================
136: 
137: label_candidates = [
138:     "true_label",
139:     "binary_label",
140:     "y",
141:     "label",
142: ]
143: 
144: label_column = None
145: 
146: for candidate in label_candidates:
147: 
148:     if candidate in df.columns:
149:         label_column = candidate
150:         break
151: 
152: if label_column is None:
153:     raise RuntimeError(
154:         "Could not find binary label column. "
155:         f"Available columns: {df.columns.tolist()}"
156:     )
157: 
158: 
159: # =====================================================================
160: # VALIDATE FEATURES
161: # =====================================================================
162: 
163: missing_features = [
164:     feature
165:     for feature in FEATURES
166:     if feature not in df.columns
167: ]
168: 
169: if missing_features:
170:     raise RuntimeError(
171:         "Missing required 14-miRNA features:\n"
172:         + "\n".join(missing_features)
173:     )
174: 
175: 
176: # =====================================================================
177: # BUILD X / y
178: # =====================================================================
179: 
180: X = df[FEATURES].to_numpy(
181:     dtype=np.float64
182: )
183: 
184: y = df[label_column].to_numpy(
185:     dtype=int
186: )
187: 
188: if len(df) != 179:
189:     raise RuntimeError(
190:         f"Expected 179 patients, found {len(df)}."
191:     )
192: 
193: if len(np.unique(y)) != 2:
194:     raise RuntimeError(
195:         "Expected exactly two classes."
196:     )
197: 
198: if not np.all(np.isfinite(X)):
199:     raise RuntimeError(
200:         "Expression matrix contains non-finite values."
201:     )
202: 
203: print()
204: print("=" * 90)
205: print("DATA VALIDATION")
206: print("=" * 90)
207: 
208: print()
209: print(f"Patients: {len(df)}")
210: print(f"Features: {X.shape[1]}")
211: print(f"Cancer: {(y == 1).sum()}")
212: print(f"Non-cancer: {(y == 0).sum()}")
213: 
214: 
215: # =====================================================================
216: # CROSS-VALIDATION
217: # =====================================================================
218: 
219: cv = StratifiedKFold(
220:     n_splits=5,
221:     shuffle=True,
222:     random_state=2026,
223: )
224: 
225: oof_probability = np.full(
226:     len(df),
227:     np.nan,
228: )
229: 
230: oof_prediction = np.full(
231:     len(df),
232:     -1,
233:     dtype=int,
234: )
235: 
236: fold_records = []
237: 
238: 
239: print()
240: print("=" * 90)
241: print("5-FOLD STRATIFIED PATIENT-LEVEL CROSS-VALIDATION")
242: print("=" * 90)
243: 
244: 
245: for fold, (train_idx, test_idx) in enumerate(
246:     cv.split(X, y),
247:     start=1,
248: ):
249: 
250:     X_train = X[train_idx]
251:     X_test = X[test_idx]
252: 
253:     y_train = y[train_idx]
254:     y_test = y[test_idx]
255: 
256:     # ---------------------------------------------------------------
257:     # Pipeline prevents leakage:
258:     #
259:     # scaler is fitted ONLY on training data.
260:     # ---------------------------------------------------------------
261: 
262:     model = Pipeline(
263:         steps=[
264:             (
265:                 "scaler",
266:                 StandardScaler(),
267:             ),
268:             (
269:                 "classifier",
270:                 LogisticRegression(
271:                     penalty="l2",
272:                     C=1.0,
273:                     solver="liblinear",
274:                     max_iter=5000,
275:                     random_state=2026,
276:                 ),
277:             ),
278:         ]
279:     )
280: 
281:     model.fit(
282:         X_train,
283:         y_train,
284:     )
285: 
286:     probability = model.predict_proba(
287:         X_test
288:     )[:, 1]
289: 
290:     prediction = (
291:         probability >= 0.5
292:     ).astype(int)
293: 
294:     oof_probability[test_idx] = probability
295:     oof_prediction[test_idx] = prediction
296: 
297:     fold_auc = roc_auc_score(
298:         y_test,
299:         probability,
300:     )
301: 
302:     fold_accuracy = accuracy_score(
303:         y_test,
304:         prediction,
305:     )
306: 
307:     fold_records.append(
308:         {
309:             "fold": fold,
310:             "n_train": len(train_idx),
311:             "n_test": len(test_idx),
312:             "cancer_train": int((y_train == 1).sum()),
313:             "noncancer_train": int((y_train == 0).sum()),
314:             "cancer_test": int((y_test == 1).sum()),
315:             "noncancer_test": int((y_test == 0).sum()),
316:             "roc_auc": float(fold_auc),
317:             "accuracy": float(fold_accuracy),
318:         }
319:     )
320: 
321:     print()
322:     print(
323:         f"Fold {fold}: "
324:         f"train={len(train_idx)}, "
325:         f"test={len(test_idx)}, "
326:         f"AUC={fold_auc:.6f}, "
327:         f"accuracy={fold_accuracy:.6f}"
328:     )
329: 
330: 
331: # =====================================================================
332: # OUT-OF-FOLD METRICS
333: # =====================================================================
334: 
335: if np.any(~np.isfinite(oof_probability)):
336:     raise RuntimeError(
337:         "Some patients did not receive an out-of-fold prediction."
338:     )
339: 
340: if np.any(oof_prediction < 0):
341:     raise RuntimeError(
342:         "Some patients did not receive an out-of-fold class prediction."
343:     )
344: 
345: 
346: oof_auc = roc_auc_score(
347:     y,
348:     oof_probability,
349: )
350: 
351: oof_pr_auc = average_precision_score(
352:     y,
353:     oof_probability,
354: )
355: 
356: oof_accuracy = accuracy_score(
357:     y,
358:     oof_prediction,
359: )
360: 
361: cm = confusion_matrix(
362:     y,
363:     oof_prediction,
364: )
365: 
366: tn, fp, fn, tp = cm.ravel()
367: 
368: sensitivity = (
369:     tp / (tp + fn)
370:     if (tp + fn) > 0
371:     else float("nan")
372: )
373: 
374: specificity = (
375:     tn / (tn + fp)
376:     if (tn + fp) > 0
377:     else float("nan")
378: )
379: 
380: 
381: # =====================================================================
382: # SAVE PREDICTIONS
383: # =====================================================================
384: 
385: predictions = df.copy()
386: 
387: predictions["oof_probability_cancer"] = (
388:     oof_probability
389: )
390: 
391: predictions["oof_prediction"] = (
392:     oof_prediction
393: )
394: 
395: predictions["fold"] = -1
396: 
397: for fold, (_, test_idx) in enumerate(
398:     cv.split(X, y),
399:     start=1,
400: ):
401:     predictions.loc[
402:         test_idx,
403:         "fold",
404:     ] = fold
405: 
406: 
407: OUTPUT_DIR.mkdir(
408:     parents=True,
409:     exist_ok=True,
410: )
411: 
412: predictions.to_csv(
413:     PREDICTIONS_OUTPUT,
414:     index=False,
415: )
416: 
417: 
418: # =====================================================================
419: # SAVE RESULTS
420: # =====================================================================
421: 
422: results = {
423:     "model": "standardized L2 logistic regression",
424:     "features": FEATURES,
425:     "n_patients": int(len(df)),
426:     "n_features": int(len(FEATURES)),
427:     "cancer_count": int((y == 1).sum()),
428:     "noncancer_count": int((y == 0).sum()),
429:     "cross_validation": {
430:         "method": "5-fold stratified patient-level CV",
431:         "shuffle": True,
432:         "random_state": 2026,
433:     },
434:     "metrics": {
435:         "roc_auc": float(oof_auc),
436:         "pr_auc": float(oof_pr_auc),
437:         "accuracy": float(oof_accuracy),
438:         "sensitivity_at_0.5": float(sensitivity),
439:         "specificity_at_0.5": float(specificity),
440:     },
441:     "confusion_matrix": {
442:         "TN": int(tn),
443:         "FP": int(fp),
444:         "FN": int(fn),
445:         "TP": int(tp),
446:     },
447:     "folds": fold_records,
448: }
449: 
450: with RESULTS_OUTPUT.open(
451:     "w"
452: ) as f:
453: 
454:     json.dump(
455:         results,
456:         f,
457:         indent=2,
458:     )
459: 
460: 
461: # =====================================================================
462: # REPORT
463: # =====================================================================
464: 
465: print()
466: print("=" * 90)
467: print("MODERN 14-miRNA BASELINE RESULTS")
468: print("=" * 90)
469: 
470: print()
471: print(f"Patients:       {len(df)}")
472: print(f"Cancer:         {(y == 1).sum()}")
473: print(f"Non-cancer:     {(y == 0).sum()}")
474: 
475: print()
476: print(f"OOF ROC AUC:    {oof_auc:.6f}")
477: print(f"OOF PR AUC:     {oof_pr_auc:.6f}")
478: print(f"OOF Accuracy:   {oof_accuracy:.6f}")
479: print(f"Sensitivity:    {sensitivity:.6f}")
480: print(f"Specificity:    {specificity:.6f}")
481: 
482: print()
483: print("Confusion matrix:")
484: print(
485:     f"  TP={tp} FN={fn}"
486: )
487: print(
488:     f"  FP={fp} TN={tn}"
489: )
490: 
491: print()
492: print("=" * 90)
493: print("FILES")
494: print("=" * 90)
495: 
496: print()
497: print("OOF predictions:")
498: print(PREDICTIONS_OUTPUT)
499: 
500: print()
501: print("Results:")
502: print(RESULTS_OUTPUT)
503: 
504: print()
505: print("=" * 90)
506: print("IMPORTANT")
507: print("=" * 90)
508: 
509: print()
510: print(
511:     "This is a NEW modern model. It is NOT the historical "
512:     "eLife neural network."
513: )
514: 
515: print(
516:     "All scaling was performed separately inside each training fold."
517: )
518: 
519: print(
520:     "The out-of-fold AUC is therefore a patient-level cross-validated "
521:     "estimate rather than a resubstitution/full-cohort score."
522: )
~~~
