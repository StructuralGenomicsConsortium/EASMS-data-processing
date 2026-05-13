# Data Processing Pipeline

What happens to a raw ASMS results CSV once Quality Check passes. Steps 1–9 below run sequentially per input file; each step's output is persisted on disk so any step can be resumed from a saved checkpoint.

For QC details (Step 0), see [QUALITY_CHECKS.md](QUALITY_CHECKS.md). For run instructions and the `--start-from` / `--end-at` flags, see [Running a Subset of Steps](#running-a-subset-of-steps) below or [USAGE.md](USAGE.md).

## Pipeline Steps

The entry point is [src/Main.py](src/Main.py). It iterates over every CSV in `RawData/` and runs Step 0 (QC) first, then the steps below for files that pass QC.

### 1. Split by target — [`separate_protein_files.split_protein_data`](src/separate_protein_files.py)

Groups rows in the raw CSV by `TARGET_ID` and writes one CSV per target into `ProcessedData_<csv_basename>/Separated_Files/`. Requires `PROTEIN_NUMBER`, `ASMS_BATCH_NAME`, and `TARGET_ID` columns.

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

Looks up the master-list file for the current raw CSV via `MasterList_Information.xlsx`, loads it, and adds compounds from the master list as negative samples. Each master list must contain a `SMILES` column.

### 6. Generate ML labels — [`produce_ml_labels.generate_ml_labels`](src/produce_ml_labels.py)

Assigns an `AIRCHECK_LABEL` integer based on `EASMS_ENRICHMENT`, `PVALUE`, `ISOMERS`, and `HAD_DUPLICATE_INTENSITY` (values range from −2 to 4 — see the module's docstring for the exact rules). This is the last CSV-format step.

### 7. Extract fingerprints + rename + binary label — [`fingerprint_extraction.extract_fingerprints`](src/fingerprint_extraction.py)

Three transformations applied together:

- **Fingerprints / descriptors** (via [src/fingerprints.py](src/fingerprints.py), [src/utils.py](src/utils.py)): `MW`, `ALOGP`, and `ECFP4`, `ECFP6`, `FCFP4`, `FCFP6`, `MACCS`, `RDK`, `AVALON`, `TOPTOR`, `ATOMPAIR`.
- **Column renames** (inline in `Main.py`): `TARGET_VALUE` → `TARGET_INTENSITY_VALUE`, `MEAN_NONTARGET_VALUES` → `NONTARGET_INTENSITY_VALUE`.
- **Binary label**: `LABEL = 1 if BINARY_LABEL == "Y" else 0`.

This is the first Parquet-format step (CSV is dropped because the wide fingerprint columns make it slow and huge).

### 8. Select full column set — [`column_selection.select_final_columns`](src/column_selection.py)

Reads the Step 7 output and selects `DesiredColumns` (46 cols: all metadata + scores + labels + fingerprints). Saved as Parquet.

### 9. Select key column set — [`column_selection.select_final_columns`](src/column_selection.py)

Reads the **same Step 7 output** (parallel branch — *not* downstream of Step 8) and selects `DesiredColumns2` (19 cols: IDs, target info, scores, label, MW, ALOGP, fingerprints). Saved as Parquet.

## Running a Subset of Steps

Use `--start-from N` and `--end-at N` to resume from or stop at a specific step. Earlier steps that have already been saved are loaded from disk; later steps are skipped.

```powershell
# Re-run only fingerprint extraction onward (steps 1-6 are loaded from disk)
python src/Main.py --start-from 7

# Run only step 9 (re-derive the key column subset) using Step 7's saved Parquet
python src/Main.py --start-from 9 --end-at 9

# Run just the early cleaning (steps 1-5), stop before label generation
python src/Main.py --end-at 5
```

`--start-from` accepts 0–9 (0 = Quality Check), `--end-at` accepts 0–9. Defaults: `--start-from 0 --end-at 9` (run everything, including QC). Step 0 (QC) runs only when `--start-from 0`.

## Output Layout

For each input CSV, the pipeline creates one `ProcessedData_<csv_basename>/` folder inside `--output-dir` (which defaults to `--input-dir`). Each step's output lives in its own folder so any step can be re-run from saved checkpoints:

```
ProcessedData_<csv_basename>/
├── QCaircheck<YYYYMMDD>_<csv_basename>.log    # step 0 plain-text log (date the check was run)
├── QCaircheck<YYYYMMDD>_<csv_basename>.xlsx   # step 0 same data in Excel (color-coded)
├── Step1_Separated/             # step 1 — split by target           (CSV)
│   └── <target>.csv
├── Step2_WithScores/            # step 2 — score columns added       (CSV)
├── Step3_AnomalyFiltered/       # step 3 — anomalies resolved        (CSV)
├── Step4_IsomerHandled/         # step 4 — isomers split             (CSV)
├── Step5_WithNegatives/         # step 5 — masterlist negatives      (CSV)
├── Step6_MLReady/               # step 6 — labels added              (CSV)
├── Step7_WithFingerprints/      # step 7 — FPs + rename + LABEL      (Parquet)
│   └── <target>.parquet
├── Step8_FullColumns/           # step 8 — full column subset        (Parquet)
└── Step9_KeyColumns/            # step 9 — slim column subset        (Parquet)
```

Step 0 may also write supplementary report CSVs (`FullyDuplicate_rows_report.csv`, `invalid_smiles_report.csv`, etc.) alongside the QC log files when checks find issues — see [QUALITY_CHECKS.md](QUALITY_CHECKS.md) for the full list.
