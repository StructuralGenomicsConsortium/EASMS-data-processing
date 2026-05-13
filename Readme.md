# ASMS Data Curation Pipeline

This repository contains a Python-based data curation pipeline for processing Affinity Selection Mass Spectrometry (ASMS) datasets. The pipeline prepares data for machine learning by performing quality checks, cleaning, labeling, and fingerprint extraction.

## Main Features

- 100+ pre-processing quality checks (file format, filename format, row content, per-column rules) with plain-text + Excel logs
- Splits protein-specific data into separate files
- Detects and filters out anomalous entries
- Handles isomer corrections
- Adds negative samples from a master list
- Generates binary labels for machine learning
- Extracts chemical fingerprints (e.g., ECFP4, FCFP6, MACCS)
- Saves curated data in both CSV and Parquet formats

## Documentation

- **[USAGE.md](USAGE.md)** — environment setup, run commands, `--start-from` / `--end-at` flags.
- **[QUALITY_CHECKS.md](QUALITY_CHECKS.md)** — Step 0 (QC): every check, severity, and report file produced.
- **[PIPELINE.md](PIPELINE.md)** — Steps 1–9 (data processing): what each step does, output layout, resuming.

## Requirements

Before running the pipeline, make sure all of the following are in place. Items 2–5 live inside `--input-dir` (defaults to the current working directory).

1. **Python environment with dependencies installed** — set up `.venv` and install `requirements.txt`. Step-by-step instructions in [USAGE.md §1](USAGE.md).

2. **`RawData/` folder** — one or more ASMS results CSV files, each named in the convention `asms_<provider>_<batch>_<library>_<YYYYMMDD>.csv`. Every CSV in this folder is processed.

3. **`MasterLists/` folder** containing:
   - `MasterList_Information.xlsx` — required mapping file (columns `FileName`, `MaterListName`) that links each raw CSV to its master list.
   - One `.xlsx` per registered compound library referenced above (each must contain at least a `SMILES` column).

4. **`Providers.csv`** — list of valid provider acronyms and data-generator names. The real file is gitignored (private company info); copy [Providers_sample.csv](Providers_sample.csv) to `Providers.csv` and replace the entries with real values. Required columns: `acronym`, `name`, `data_generator_name`.

5. **`ASMS Meta Data.csv`** — canonical column-name reference. The first row lists every column name a raw ASMS CSV must contain; the second row holds data types (informational only). QC fails a file when its columns don't exactly match this list.

Expected layout:

```
<input-dir>/
├── RawData/
│   ├── asms_<provider>_<batch>_<library>_<YYYYMMDD>.csv
│   └── ...
├── MasterLists/
│   ├── MasterList_Information.xlsx
│   ├── <library1>.xlsx
│   └── <library2>.xlsx
├── Providers.csv
└── ASMS Meta Data.csv
```

## Data Inputs (file formats)

### `RawData/`

One or more **ASMS results CSV files**. Each row is a compound–protein measurement with target/non-target intensities, replicates, pool info, and protein metadata. Every CSV in this folder is processed. The required column names are defined by `ASMS Meta Data.csv` (see below).

### `MasterLists/`

Excel files describing the compound libraries used in the screen. This folder must contain:

- **`MasterList_Information.xlsx`** (required). Maps each raw-data CSV to its corresponding master list. Must have two columns:
  - `FileName` — the filename of a CSV in `RawData/` (e.g. `asms_acmecorp_01_Chemdiv9k_20260512.csv`)
  - `MaterListName` — the base name (no extension) of the matching master list `.xlsx` file in `MasterLists/`

- **One `.xlsx` per master list referenced above** (e.g. `Chemdiv9k.xlsx`). Each must contain at least a `SMILES` column and a `formula` column. Used to draw negative samples for the model and to validate input SMILES / formulas.

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

### `ASMS Meta Data.csv`

The **canonical column-name reference** for raw ASMS results files. The QC step (Check 7) reads it and compares the columns of each raw CSV against this list — files with missing or extra columns fail QC and are skipped.

Format:

- **Row 1 (header)** — the canonical column names every raw CSV is expected to have (e.g. `COMPOUND_ID, SMILES, ASMS_BATCH_NAME, COMPOUND_FORMULA, POOL_NAME, ...`).
- **Row 2** — the data type per column (`VARCHAR`, `INT`, `FLOAT`, `BOOL`). Only row 1 is used by the checker; row 2 is informational.

Only column **names** are compared (not types and not order). Whitespace around names is stripped and accidental duplicate columns are collapsed, so a stray trailing space won't cause a false failure.

To change which columns are required, edit `ASMS Meta Data.csv` directly — no code change needed.

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

Your real `RawData/`, `MasterLists/`, and the generated `ProcessedData_*/` folders are all gitignored — only the `_sample` versions are tracked in this repo.
