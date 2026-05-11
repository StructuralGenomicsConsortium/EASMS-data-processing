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
