# ASMS Data Curation Pipeline

This repository contains a Python-based data curation pipeline for processing Affinity Selection Mass Spectrometry (ASMS) datasets. The pipeline prepares data for machine learning by performing cleaning, labeling, and fingerprint extraction.

## Main Features

- Splits protein-specific data into separate files
- Detects and filters out anomalous entries
- Handles isomer corrections
- Adds negative samples from a master list
- Generates binary labels for machine learning
- Extracts chemical fingerprints (e.g., ECFP4, FCFP6, MACCS)
- Saves curated data in both CSV and Parquet formats

## Pipeline Steps

The entry point is [src/Main.py](src/Main.py). It iterates over every CSV in `RawData/` and runs the steps below per file. Each step is implemented in its own module so it can be edited or reused on its own.

### 1. Split by target — [`separate_protein_files.split_protein_data`](src/separate_protein_files.py)

Groups rows in the raw CSV by `TARGET_ID` and writes one CSV per target into `ProcessedData_<csv_basename>/Separated_Files/`. Requires `PROTEIN_NUMBER`, `ASMS_BATCH_NUM`, and `TARGET_ID` columns.

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

Assigns an `AIRCHECK_LABEL` integer based on `EASMS_ENRICHMENT`, `PVALUE`, `ISOMERS`, and `HAD_DUPLICATE_INTENSITY` (values range from −2 to 4 — see the module's docstring for the exact rules). The DataFrame at this point is saved to `MLReady/<target>.csv` and `.parquet`.

### 7. Extract chemical fingerprints — [`fingerprint_extraction.extract_fingerprints`](src/fingerprint_extraction.py)

Uses RDKit (via [src/fingerprints.py](src/fingerprints.py) and [src/utils.py](src/utils.py)) to add per-compound molecular descriptors and fingerprints:
- Descriptors: `MW`, `ALOGP`
- Fingerprints: `ECFP4`, `ECFP6`, `FCFP4`, `FCFP6`, `MACCS`, `RDK`, `AVALON`, `TOPTOR`, `ATOMPAIR`

### 8. Rename / derive columns — inline in [src/Main.py](src/Main.py)

- `TARGET_VALUE` → `TARGET_INTENSITY_VALUE`
- `MEAN_NONTARGET_VALUES` → `NONTARGET_INTENSITY_VALUE`
- Adds `LABEL = 1 if BINARY_LABEL == "Y" else 0`

### 9. Select final columns and save — [`column_selection.select_final_columns`](src/column_selection.py)

Applied twice with two different column lists defined in `Main.py`:

- **`DesiredColumns`** (46 cols — all metadata + scores + labels + fingerprints) → `MLReady_FullColumns/<target>.csv` and `.parquet`
- **`DesiredColumns2`** (19 cols — slim model-input set: IDs, target info, scores, label, MW, ALOGP, fingerprints) → `MLReady_KeyColumns/<target>.csv` and `.parquet`

## Output Layout

For each input CSV, the pipeline creates one `ProcessedData_<csv_basename>/` folder at the dataset root:

```
ProcessedData_<csv_basename>/
├── Separated_Files/         # output of step 1
│   └── <target>.csv
├── MLReady/                 # output of step 6 (labeled, no fingerprints)
│   ├── <target>.csv
│   └── <target>.parquet
├── MLReady_FullColumns/     # output of step 9 (full metadata + fingerprints)
│   ├── <target>.csv
│   └── <target>.parquet
└── MLReady_KeyColumns/      # output of step 9 (slim model-input + fingerprints)
    ├── <target>.csv
    └── <target>.parquet
```

## Data Inputs

The pipeline expects two input folders at the dataset root (the path given via `--path`, or the parent of the current working directory if `--path` is omitted). For instructions on how to run the pipeline, see [USAGE.md](USAGE.md).

### `RawData/`

One or more **ASMS results CSV files**, for example `ASMS_results_2_all.csv`. Each row is a compound–protein measurement with target/non-target intensities, replicates, pool info, and protein metadata. Every CSV in this folder is processed.

### `MasterLists/`

Excel files describing the compound libraries used in the screen. This folder must contain:

- **`MasterList_Information.xlsx`** (required). Maps each raw-data CSV to its corresponding master list. Must have two columns:
  - `FileName` — the filename of a CSV in `RawData/` (e.g. `ASMS_results_2_all.csv`)
  - `MaterListName` — the base name (no extension) of the matching master list `.xlsx` file in `MasterLists/`

- **One `.xlsx` per master list referenced above** (e.g. `Chemdiv+Chiral6k_15k.xlsx`, `Chemdiv_9k.xlsx`). Each must contain at least a `SMILES` column; it's used to draw negative samples for the model.

## Sample Data

For reference, the repo includes two small placeholder folders:

- [RawData_sample/](RawData_sample/)
- [MasterLists_sample/](MasterLists_sample/)

These show the expected file layout and naming. **They are not picked up by the pipeline automatically** — `Main.py` only reads from `RawData/` and `MasterLists/`. To use them, either:

1. Rename the folders by dropping the `_sample` suffix:
   ```powershell
   Rename-Item RawData_sample RawData
   Rename-Item MasterLists_sample MasterLists
   ```
2. Or copy/move the sample files into your own `RawData/` and `MasterLists/` folders.

Your real `RawData/`, `MasterLists/`, and the generated `ProcessedData/` are all gitignored — only the `_sample` versions are tracked in this repo.
