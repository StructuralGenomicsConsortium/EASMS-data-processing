# ASMS Data Curation Pipeline

This repository contains a Python-based data curation pipeline for processing Affinity Selection Mass Spectrometry (ASMS) datasets. The pipeline prepares data for machine learning by performing quality checks, cleaning, labeling, and fingerprint extraction.

## Main Features

- 100+ pre-processing quality checks (file format, filename format, row content, per-column rules) with plain-text + Excel logs
- Runs on a **single file or a folder** of files, with paths on the **local disk or Google Cloud Storage** (`gs://`) — auto-detected from the path
- Splits protein-specific data into separate files
- Detects and filters out anomalous entries
- Handles isomer corrections
- Adds negative samples from a master list
- Generates binary labels for machine learning
- Extracts chemical fingerprints (e.g., ECFP4, FCFP6, MACCS)
- Saves curated data in both CSV and Parquet formats

## Documentation

- **[USAGE.md](USAGE.md)** — environment setup, run commands, `--start-from` / `--end-at` flags, and **running on Google Cloud Storage (GCP)**.
- **[QUALITY_CHECKS.md](QUALITY_CHECKS.md)** — Step 0 (input QC): every check, severity, and report file produced.
- **[PIPELINE.md](PIPELINE.md)** — Steps 1–8 (data processing): what each step does, output layout, resuming.
- **[POST_QC.md](POST_QC.md)** — Post-pipeline QC: 23 checks that run after Step 7 to catch regressions in the pipeline's own output.

## Requirements

The pipeline needs two kinds of input: the **raw data** to process, and a few **config/reference files**. They are specified independently — you are no longer required to lay everything out inside one folder.

1. **Python environment with dependencies installed** — set up `.venv` and install `requirements.txt` (includes `fsspec`/`gcsfs` for GCS support). Step-by-step instructions in [USAGE.md §1](USAGE.md).

2. **Raw data** — one ASMS results CSV, or a folder of them. Each file is named in the convention `asms_<provider>_<batch>_<library>_<YYYYMMDD>.csv`. Point at it with **`--input-file`** (a single CSV) or **`--input-dir`** (a folder; every `*.csv` in it is processed). Paths may be local or `gs://`.

3. **`MasterLists/` folder** — one file per compound library, named `<LIBRARY_NAME>.xlsx` **or** `<LIBRARY_NAME>.csv` (each must contain at least a `SMILES` column, ideally `formula` too). The library for a given raw file is resolved directly from its `LIBRARY_NAME` column — **no mapping file needed**. Default location: this repo's `MasterLists/`; override with `--masterlists-dir`.

4. **`Providers.csv`** — valid provider acronyms and data-generator names (columns `acronym`, `name`, `data_generator_name`). The real file is gitignored (private company info); copy [Providers_sample.csv](Providers_sample.csv) to `Providers.csv` and fill in real values. Default location: this repo root; override with `--providers-csv`.

5. **`RawDataColumns.csv`** — canonical column-name reference. Row 1 lists every column name a raw CSV must contain; row 2 holds data types (informational only). QC fails a file when its columns don't match this list. Default location: this repo root; override with `--meta-csv`.

6. **`ColumnActions.xlsx`** — drives the Step 8 column split. Two columns: `Column name` and `Action`. Columns tagged **`data`** go to the final ML data file (`Step8_FinalData/`, Parquet); columns tagged **`metadata`** go to a metadata CSV (`Step8_Metadata/`); columns tagged **`-`** are dropped. Edit this file to change what lands where — no code change needed. Default location: this repo root; override with `--column-actions`.

By default, the **config/reference files (3–6) are read from this repo** — so you only need to point `--input-file`/`--input-dir` at the data. Each can be overridden to a shared folder or a `gs://` path. See [USAGE.md](USAGE.md) for the full flag list and examples.

## Data Inputs (file formats)

### Raw data (`--input-file` / `--input-dir`)

**ASMS results CSV files.** Each row is a compound–protein measurement with target/non-target intensities, replicates, pool info, and protein metadata. Pass a single file with `--input-file`, or a folder of CSVs with `--input-dir` (every `*.csv` in it is processed — no special subfolder name required). The required column names are defined by `RawDataColumns.csv` (see below).

### `MasterLists/`

One file per compound library used in the screen, named after the library:

- **`<LIBRARY_NAME>.xlsx` or `<LIBRARY_NAME>.csv`** (e.g. `Chemdiv9k.csv`). Each must contain at least a `SMILES` column and a `formula` column; optional `COMPOUND_ID` / `SGC ID for Component` and `SGC ID for Pool` columns are copied onto negative samples when present. Used to draw negative samples and to validate input SMILES / formulas.

The library for a given raw file is found by reading its **`LIBRARY_NAME`** column and loading `<LIBRARY_NAME>.xlsx`/`.csv` from `MasterLists/` (`.xlsx` wins if both exist). There is **no `MasterList_Information.xlsx` mapping file** — just make sure the matching library file exists.

### `Providers.csv`

Three columns:

```csv
acronym,name,data_generator_name
acmecorp,Acme Corp Research Labs,ASMS_ACME_CORP
fakelab,FakeLab Pharmaceuticals Inc,ASMS_FAKELAB
genericrx,GenericRx Therapeutics,ASMS_GENERICRX
```

- `acronym` — the `<provider>` segment of raw CSV filenames and the prefix of `ASMS_BATCH_NAME` values.
- `data_generator_name` — the exact value the `DATA_GENERATOR_NAME` column must contain.

This file is gitignored; [Providers_sample.csv](Providers_sample.csv) has placeholder values. Copy it to `Providers.csv` and replace with real entries.

### `RawDataColumns.csv`

The **canonical column-name reference** for raw ASMS results files. The QC step (Check 7) reads it and compares the columns of each raw CSV against this list — files with missing or extra columns fail QC and are skipped.

Format:

- **Row 1 (header)** — the canonical column names every raw CSV is expected to have (e.g. `COMPOUND_ID, SMILES, ASMS_BATCH_NAME, COMPOUND_FORMULA, POOL_NAME, ...`).
- **Row 2** — the data type per column (`VARCHAR`, `INT`, `FLOAT`, `BOOL`). Only row 1 is used by the checker; row 2 is informational.

Only column **names** are compared (not types and not order). Whitespace around names is stripped and accidental duplicate columns are collapsed, so a stray trailing space won't cause a false failure.

To change which columns are required, edit `RawDataColumns.csv` directly — no code change needed.

## Sample Data

For reference, the repo includes two small placeholder folders that show the expected file layout and naming:

- [RawData_sample/](RawData_sample/) — example raw CSV(s); point `--input-dir` at this folder (or `--input-file` at a file inside it) to try the pipeline.
- [MasterLists_sample/](MasterLists_sample/) — example library file(s); pass `--masterlists-dir MasterLists_sample` to use them.

Your real raw-data and `MasterLists/` folders and the generated `ProcessedData_*/` folders are gitignored — only the `_sample` versions are tracked in this repo.

## What's new

Recent changes to the pipeline:

- **Local *and* Google Cloud Storage** — pass local paths or `gs://` URLs anywhere; the code auto-detects. See the GCP section in [USAGE.md](USAGE.md).
- **File or folder input** — `--input-file` for one CSV, `--input-dir` for a folder. The old required `RawData/` subfolder is gone; output `ProcessedData_<name>/` is created next to the input (or at `--output-dir`).
- **Config as input** — `MasterLists/`, `Providers.csv`, `RawDataColumns.csv`, `ColumnActions.xlsx` default to this repo and are each overridable (`--masterlists-dir`, `--providers-csv`, `--meta-csv`, `--column-actions`), including to `gs://`.
- **Masterlist resolution simplified** — the library is read from the `LIBRARY_NAME` column; `MasterList_Information.xlsx` is no longer used. Master lists may be `.xlsx` **or** `.csv`.
- **Final-step redesign** — the old Step 8/9 `DesiredColumns` selections are replaced by a single Step 8 driven by `ColumnActions.xlsx`: `data`-tagged columns → `Step8_FinalData/` (Parquet), `metadata`-tagged columns → `Step8_Metadata/` (CSV), `-` dropped. Step 7 remains the full ML-ready file; post-pipeline QC now validates Step 7. Steps 7–8 output files are named by `UNIQUE_PROTEIN_ID` (Steps 1–6 by `TARGET_ID`). `--start-from`/`--end-at` now span 0–8.
- **Per-batch protein summary** — each file gets a `proteins_in_batch.csv` (one row per `TARGET_ID`: `PROTEIN_NUMBER`, `UNIPROT_ID`, `PROTEIN_SEQ`, `TARGET_ID`).
- **Column updates** — units added to `INCUBATION_VOLUME (uL)`, `PROTEIN_CONC (uM)`, `COMPOUND_CONC (uM)`, `RT (min)`; `MS_REPRODUCIBILITY` spelling corrected; new `PROTEIN_NAME` (VARCHAR) column with checks; `PROTEIN_ID` renamed to `UNIPROT_ID`; new `UNIQUE_PROTEIN_ID` (VARCHAR, format `<names>_<keys>` with `|`-separated lists for protein complexes) column with checks.
- **TARGET_ID complexes** — `TARGET_ID` now also supports complexes: `<names>_<keys>_<ranges>` as `|`-separated lists of equal length (e.g. `EXOSC7|EXOSC6|EXOSC8_Q15024|Q5RKV6|Q96B26_1_291|1_272|1_276`).
- **CHIRAL_SELECTIVITY report** — failing rows are written to `chiral_selectivity_not_allowed_report.csv` for investigation (the QC log keeps just the count).
- **Faster fingerprints** — Step 7 extraction parses each SMILES once, batches the generators, and dedupes unique SMILES (drop-in; the original provider version is archived as `fingerprint_extraction_original.py`).
