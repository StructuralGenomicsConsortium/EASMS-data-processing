# Data Processing Pipeline

What happens to a raw ASMS results CSV once Quality Check passes. Steps 1–8 below run sequentially per input file; each step's output is persisted on disk so any step can be resumed from a saved checkpoint.

For QC details (Step 0), see [QUALITY_CHECKS.md](QUALITY_CHECKS.md). For run instructions and the `--start-from` / `--end-at` flags, see [Running a Subset of Steps](#running-a-subset-of-steps) below or [USAGE.md](USAGE.md).

## Pipeline Steps

The entry point is [src/Main.py](src/Main.py). It processes the raw CSV(s) you pass via `--input-file` (one file) or `--input-dir` (a folder of CSVs) — local or `gs://` — running Step 0 (QC) first, then the steps below for files that pass QC.

### 1. Split by target — [`separate_protein_files.split_protein_data`](src/separate_protein_files.py)

Groups rows in the raw CSV by `TARGET_ID` and writes one CSV per target into `ProcessedData_<csv_basename>/Step1_Separated/`. Requires `PROTEIN_NUMBER`, `ASMS_BATCH_NAME`, and `TARGET_ID` columns.

### 2. Compute scores — [`add_scores.compute_and_add_scores`](src/add_scores.py)

For each per-target file, this step uses the *other* per-target files in the batch as a non-target reference and adds the following columns:

- **`TARGET_VALUE`** — mean of the three replicate intensities for the current target:
  `TARGET_VALUE = mean(POS_INT_REP1, POS_INT_REP2, POS_INT_REP3)` *(skipping NaN)*.

- **`SELECTIVE_VALUE`** — per compound, the **maximum** of `TARGET_VALUE` across all *other* targets in the batch (i.e., the strongest off-target signal observed for this compound).

- **`NTC_VALUE`** — per compound, the **minimum** of `TARGET_VALUE` across all other targets (acts as a no-target-control reference: the weakest off-target signal).

- **`ENRICHMENT`** — current target's signal over the weakest off-target signal:
  `ENRICHMENT = TARGET_VALUE / NTC_VALUE`.

- **`SELECTIVE_ENRICHMENT`** — current target's signal over the strongest off-target signal:
  `SELECTIVE_ENRICHMENT = TARGET_VALUE / SELECTIVE_VALUE`. Values ≥ 1 indicate selectivity for this target.

- **`MEAN_NONTARGET_VALUES`** — for the current compound, take each *other* target file, compute the mean of its three replicates, then average those means across files. This is the "average off-target signal" for the compound.

- **`EASMS_ENRICHMENT`** — current target's signal over the mean off-target signal:
  `EASMS_ENRICHMENT = TARGET_VALUE / MEAN_NONTARGET_VALUES`.

- **`PVALUE`** — two-sample Welch's t-test ([`scipy.stats.ttest_ind`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html) with `equal_var=False`) comparing:
  - the current target's three replicates (`POS_INT_REP{1,2,3}`)
  - vs. **all** pooled replicates from the other targets for the same compound.

  Returns `None` if either side has insufficient samples (current < 1 value, or other < 3 values), and `1.0` if either side has zero variance.

Requires the three replicate columns `POS_INT_REP1`, `POS_INT_REP2`, `POS_INT_REP3` and a `COMPOUND_ID` column on every input file.

### 3. Filter anomalies — [`anomaly_selection.filter_anomalous_data`](src/anomaly_selection.py)

Resolves duplicate-SMILES rows with conflicting enrichment:
- All `ENRICHMENT < 1` → keep the smallest
- All `ENRICHMENT > 10` → keep the largest
- Mixed values → drop all of them (ambiguous)

Anomalies and removals are logged for audit.

### 4. Handle isomers — [`isomer_handling.handle_isomers`](src/isomer_handling.py)

Splits rows whose `SMILES` contains multiple isomers separated by `;` into one row per isomer, and records the original group in a new `ISOMERS` column.

### 5. Add negative samples — [`add_negatives.add_negative_samples_from_masterlist`](src/add_negatives.py)

Resolves the master-list file for the current raw CSV from its `LIBRARY_NAME` column — loads `<LIBRARY_NAME>.xlsx` or `<LIBRARY_NAME>.csv` from the master-lists folder (`--masterlists-dir`, default `<repo>/MasterLists`) — and adds its compounds as negative samples. Each master list must contain a `SMILES` column; `formula`, `COMPOUND_ID`/`SGC ID for Component`, and `SGC ID for Pool` are copied onto the negatives when present.

### 6. Generate ML labels — [`produce_ml_labels.generate_ml_labels`](src/produce_ml_labels.py)

Assigns an `AIRCHECK_LABEL` integer based on `EASMS_ENRICHMENT`, `PVALUE`, `ISOMERS`, and `HAD_DUPLICATE_INTENSITY` (values range from −2 to 4 — see the module's docstring for the exact rules). This is the last CSV-format step.

### 7. Extract fingerprints + rename + binary label — [`fingerprint_extraction.extract_fingerprints`](src/fingerprint_extraction.py)

Three transformations applied together:

- **Fingerprints / descriptors** (via [src/fingerprints.py](src/fingerprints.py), [src/utils.py](src/utils.py)): `MW`, `ALOGP`, and `ECFP4`, `ECFP6`, `FCFP4`, `FCFP6`, `MACCS`, `RDK`, `AVALON`, `TOPTOR`, `ATOMPAIR`.
- **Column renames** (inline in `Main.py`): `TARGET_VALUE` → `TARGET_INTENSITY_VALUE`, `MEAN_NONTARGET_VALUES` → `NONTARGET_INTENSITY_VALUE`.
- **Binary label**: `LABEL = 1 if BINARY_LABEL == "Y" else 0`.

This is the first Parquet-format step (CSV is dropped because the wide fingerprint columns make it slow and huge).

#### Fingerprint storage format — `TypeOfFp`

Each fingerprint column (`ECFP4`, `MACCS`, …) can be stored in one of two formats. Set the `TypeOfFp` constant near the top of the `__main__` block in [src/Main.py](src/Main.py):

| `TypeOfFp` | Stored as | When to use |
|---|---|---|
| `"array"` *(default)* | `numpy.float32` array per row — ready to feed directly into a model | New runs. Downstream code can read the column without any string-to-array conversion. |
| `"string"` | Comma-separated string per row (e.g. `"1,0,0,1,…"`) | Legacy format from earlier pipeline versions. Use if downstream code expects the string layout (e.g. `np.fromstring(x, sep=',', dtype=np.float32)`). |

Both formats survive a Parquet round-trip; `"array"` is just one step closer to the form a model consumes.

### 8. Split into final data + metadata — [`column_selection.select_final_columns`](src/column_selection.py)

Reads the Step 7 dataframe and splits its columns into two outputs, driven by the **`ColumnActions.xlsx`** tag file (columns `Column name` / `Action`):

- columns tagged **`data`** → the final ML data file, saved as **Parquet** in `Step8_FinalData/`
- columns tagged **`metadata`** → a **CSV** sidecar in `Step8_Metadata/` (same rows, so it aligns row-for-row with the data file)
- columns tagged **`-`** (or anything else) → dropped

The tag→column mapping is loaded once via [`column_selection.load_column_tags`](src/column_selection.py); change `ColumnActions.xlsx` to change what lands in each output — no code change needed. (Step 8 is the last step; the old Step 9 "key columns" branch has been removed.)

### Post-pipeline QC — [`post_quality_check.run_post_quality_checks`](src/post_quality_check.py)

After all per-target Parquet files have been written, a lightweight QC pass runs once per input CSV against the concatenated **Step 7** output (the full file — it still carries every column, including ones the final data file drops). Its job is to catch regressions in pipeline-produced columns (label values, score ranges, fingerprint lengths) — not to re-validate the raw input. **Best-effort**: only runs when `--end-at >= 7`; if Step 7 hasn't produced output it writes a one-line "skipped" log and returns.

Outputs (next to the input-QC logs):

- `PostQClog_<YYYYMMDD>_<csv_basename>.log`
- `PostQClog_<YYYYMMDD>_<csv_basename>.xlsx`

See [POST_QC.md](POST_QC.md) for the full check list (23 checks across 5 sections), tuning constants, and design notes.

## Running a Subset of Steps

Use `--start-from N` and `--end-at N` to resume from or stop at a specific step. Earlier steps that have already been saved are loaded from disk; later steps are skipped.

```powershell
# Re-run only fingerprint extraction onward (steps 1-6 are loaded from disk)
python src/Main.py --input-file run.csv --start-from 7

# Run only step 8 (re-derive the final data + metadata) using Step 7's saved Parquet
python src/Main.py --input-file run.csv --start-from 8 --end-at 8

# Run just the early cleaning (steps 1-5), stop before label generation
python src/Main.py --input-file run.csv --end-at 5
```

`--start-from` accepts 0–8 (0 = Quality Check), `--end-at` accepts 0–8. Defaults: `--start-from 0 --end-at 8` (run everything, including QC). Step 0 (QC) runs only when `--start-from 0`.

## Output Layout

For each input CSV, the pipeline creates one `ProcessedData_<csv_basename>/` folder inside `--output-dir` (which defaults to the input file's folder, or the `--input-dir` folder). On a `gs://` input/output this whole tree is written into the bucket. Each step's output lives in its own folder so any step can be re-run from saved checkpoints:

Files in **Steps 1–6** are named by `TARGET_ID`; **Steps 7–8** are named by `UNIQUE_PROTEIN_ID` (`|` in complex IDs is sanitized to `+`; falls back to `TARGET_ID` if the column is absent).

```
ProcessedData_<csv_basename>/
├── QCaircheck<YYYYMMDD>_<csv_basename>.log    # step 0 plain-text log (date the check was run)
├── QCaircheck<YYYYMMDD>_<csv_basename>.xlsx   # step 0 same data in Excel (color-coded)
├── proteins_in_batch.csv        # per-batch protein summary (one row per TARGET_ID)
├── PostQClog_<YYYYMMDD>_<csv_basename>.log    # post-pipeline QC plain-text log
├── PostQClog_<YYYYMMDD>_<csv_basename>.xlsx   # post-pipeline QC Excel (color-coded)
├── Step1_Separated/             # step 1 — split by target           (CSV, <target>.csv)
├── Step2_WithScores/            # step 2 — score columns added       (CSV)
├── Step3_AnomalyFiltered/       # step 3 — anomalies resolved        (CSV)
├── Step4_IsomerHandled/         # step 4 — isomers split             (CSV)
├── Step5_WithNegatives/         # step 5 — masterlist negatives      (CSV)
├── Step6_MLReady/               # step 6 — labels added              (CSV)
├── Step7_WithFingerprints/      # step 7 — full file, FPs + LABEL    (Parquet, <unique_protein_id>.parquet)
├── Step8_FinalData/             # step 8 — 'data'-tagged columns     (Parquet, <unique_protein_id>.parquet)
└── Step8_Metadata/              # step 8 — 'metadata'-tagged columns (CSV, <unique_protein_id>.csv)
```

Step 0 may also write supplementary report CSVs (`FullyDuplicate_rows_report.csv`, `invalid_smiles_report.csv`, `chiral_selectivity_not_allowed_report.csv`, etc.) alongside the QC log files when checks find issues — see [QUALITY_CHECKS.md](QUALITY_CHECKS.md) for the full list.
