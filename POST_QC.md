# Post-Pipeline Quality Checks

Runs once per input CSV **after** every per-target Parquet has been written by Step 8. Its job is to catch regressions in pipeline-produced columns — *not* to re-validate the raw input (that's [QUALITY_CHECKS.md](QUALITY_CHECKS.md)).

This document covers when post-QC runs, what each check verifies, and the supplementary files it produces. For the entry point in the code, see [src/post_quality_check.py](src/post_quality_check.py).

## When it runs

Inside `process_csv_files`, right after the per-target Steps 3–9 loop ends. Gated by `--end-at >= 8` (Step 8 must have produced Parquet output). If the gate isn't met — or the `Step8_FullColumns/` folder is missing/empty — post-QC writes a one-line *"skipped: no Parquet files found"* log and returns. It is **best-effort** and never blocks the rest of the pipeline.

If **input QC fails** for a given file, the pipeline skips it entirely (including post-QC) — so the absence of a `PostQClog_*` file when QC failed is expected.

## Inputs

- The directory `ProcessedData_<csv_basename>/Step8_FullColumns/` (one `.parquet` per target).
- These are concatenated row-wise into a single DataFrame, and the checks run against the whole thing.

## Outputs

Two files written next to the input-QC logs, in `ProcessedData_<csv_basename>/`:

- **Plain text log** — `PostQClog_<YYYYMMDD>_<csv_basename>.log`, grouped by section, same layout as the input-QC log (sections → numbered checks → Overall → Statistics Summary).
- **Excel companion** — `PostQClog_<YYYYMMDD>_<csv_basename>.xlsx`, one row per check with columns `Section`, `Check #`, `Criteria`, `Status`, `Detail`. Rows are color-coded (green = PASS, red = FAIL, yellow = WARN) and the header row is frozen. A second sheet, `Statistics`, mirrors the stats block from the log.

## Checks

Total: **23 checks across 5 sections.** All reuse the helpers from [src/quality_check.py](src/quality_check.py).

### Label Checks

1. **`LABEL` only contains `{0, 1}`** — FAIL. Set-membership against `{0, 1}`. Catches the silent-zero regression we hit when the binary-label convention changed (see [PIPELINE.md → Step 7](PIPELINE.md)).
2. **`AIRCHECK_LABEL` only contains `{-2, -1, 0, 1, 2, 3, 4}`** — FAIL. The full label set produced by Step 6.

### Score Range Checks

3. **`PVALUE` is in `[0, 1]`** — FAIL. Inclusive on both ends; Welch's t-test should never produce a value outside this range.
4. **`TARGET_INTENSITY_VALUE >= 0`** — FAIL.
5. **`NONTARGET_INTENSITY_VALUE >= 0`** — FAIL.
6. **`SELECTIVE_VALUE >= 0`** — FAIL.
7. **`NTC_VALUE >= 0`** — FAIL.
8. **`ENRICHMENT >= 0`** — FAIL.
9. **`EASMS_ENRICHMENT >= 0`** — FAIL.
10. **`SELECTIVE_ENRICHMENT >= 0`** — FAIL.

### Molecular Property Checks

11. **`MW` (molecular weight) is positive (> 0)** — FAIL. A real molecule can't have zero or negative mass.
12. **`ALOGP` is in `[-5, 10]` (typical range)** — WARN. Outliers (e.g. extremely lipophilic / hydrophilic compounds) can be legitimate, so this surfaces as informational rather than fatal. Bounds are `ALOGP_MIN` / `ALOGP_MAX` at the top of [src/post_quality_check.py](src/post_quality_check.py).

### Flag Column Checks

13. **`MassSpec_Detected` only contains `{Y, N}`** — FAIL. Set-membership.
14. **`HAD_DUPLICATE_INTENSITY` only contains `{Y, N}`** — FAIL. Set-membership.

### Fingerprint Length Checks

Each row's fingerprint vector must have exactly the expected number of bits/features. The check works on both fingerprint storage formats — numpy arrays (`fp_format="array"`, the default) and legacy comma-separated strings (`fp_format="string"`).

15. **`ECFP4` length = 2048** — FAIL.
16. **`ECFP6` length = 2048** — FAIL.
17. **`FCFP4` length = 2048** — FAIL.
18. **`FCFP6` length = 2048** — FAIL.
19. **`MACCS` length = 167** — FAIL.
20. **`RDK` length = 2048** — FAIL.
21. **`AVALON` length = 2048** — FAIL.
22. **`TOPTOR` length = 2048** — FAIL.
23. **`ATOMPAIR` length = 2048** — FAIL.

Expected lengths come from the `_dimension` attribute of each `HitGen*FPFunc` class in [src/fingerprints.py](src/fingerprints.py). If a future fingerprint changes dimension, update `FP_EXPECTED_LENGTHS` at the top of [src/post_quality_check.py](src/post_quality_check.py).

## What is deliberately not included

- **Recomputation checks** like *"verify `TARGET_INTENSITY_VALUE == mean(POS_INT_REP1..3)`"*. Those would re-run the pipeline's math against the same code that just produced the value — the check would only fail if the same code disagreed with itself a second time.
- **Re-validation of raw input columns.** That happens before Step 1 in [QUALITY_CHECKS.md](QUALITY_CHECKS.md).

## Tuning

Constants at the top of [src/post_quality_check.py](src/post_quality_check.py):

- `LABEL_VALUES`, `AIRCHECK_LABEL_VALUES` — discrete allowed sets.
- `PVALUE_MIN`, `PVALUE_MAX` — p-value bounds (inclusive).
- `ALOGP_MIN`, `ALOGP_MAX` — typical ALOGP range (used for the WARN check).
- `FP_EXPECTED_LENGTHS` — fingerprint vector dimensions, one entry per fingerprint.
- `YN_COLUMNS` — list of `{Y, N}` flag columns.

Adding a new column-content check is the same pattern as the input-QC `SECTIONS` — add a `(description, function)` tuple to the relevant section in `SECTIONS` at the bottom of the file. The orchestrator picks it up on the next run; no other wiring needed.
