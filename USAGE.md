# Usage

## 1. First-time setup

Create a virtual environment and install dependencies:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

Other shells:

```bash
# Git Bash
python -m venv .venv
source .venv/Scripts/activate
pip install -r requirements.txt
```

```cmd
:: Windows cmd.exe
python -m venv .venv
.venv\Scripts\activate.bat
pip install -r requirements.txt
```

## 2. Activate the environment (every new terminal)

```powershell
.\.venv\Scripts\Activate.ps1
```

## 3. Run the pipeline

You give the pipeline **one** raw-data input — a single file or a folder:

| Flag | Meaning |
|---|---|
| `--input-file` | A single raw ASMS CSV. Output `ProcessedData_<name>/` is created in the file's own folder. |
| `--input-dir` | A folder of raw CSVs (every `*.csv` in it is processed; no `RawData/` subfolder needed). |
| `--output-dir` | Where `ProcessedData_*/` is written. Optional — defaults to the input file's folder / the input dir. |

The **config/reference files default to this repo** and are each overridable:

| Flag | Default |
|---|---|
| `--masterlists-dir` | `<repo>/MasterLists` |
| `--providers-csv` | `<repo>/Providers.csv` |
| `--meta-csv` | `<repo>/RawDataColumns.csv` |
| `--column-actions` | `<repo>/ColumnActions.xlsx` (drives the Step 8 data/metadata column split) |

```powershell
# One file — output lands next to it, config read from the repo
python src/Main.py --input-file "D:\my\dataset\asms_acme_01_lib_20260101.csv"

# A folder of CSVs
python src/Main.py --input-dir "D:\my\dataset\raw"

# Separate output location
python src/Main.py --input-file "D:\my\asms_run.csv" --output-dir "D:\my\results"

# Pull config from a shared folder instead of the repo
python src/Main.py --input-file "D:\my\asms_run.csv" --masterlists-dir "D:\shared\MasterLists"
```

Help text:

```powershell
python src/Main.py --help
```

## 4. Run only a subset of steps

Use `--start-from N` and `--end-at N` to control which steps execute. Step numbers are 1–8 (see [PIPELINE.md](PIPELINE.md) for what each step does). Skipped earlier steps are loaded from their saved output on disk. (Examples use `--input-file`; `--input-dir` works the same way, and both accept `gs://` paths.)

```powershell
# Run only the Quality Check (step 0) and stop
python src/Main.py --input-file run.csv --end-at 0

# Run only steps 1 and 2
python src/Main.py --input-file run.csv --end-at 2

# Run only step 1
python src/Main.py --input-file run.csv --end-at 1

# Re-run from fingerprint extraction onward (steps 1-6 are loaded from disk)
python src/Main.py --input-file run.csv --start-from 7

# Run exactly one step (e.g. step 5)
python src/Main.py --input-file run.csv --start-from 5 --end-at 5
```

Defaults: `--start-from 0 --end-at 8` (run everything, including QC). Quality checks (step 0) run only when `--start-from 0`.

---

# ☁️ Running on Google Cloud Storage (GCP)

The pipeline can read **and** write directly from a GCP bucket — no manual download/upload. Any path you pass (`--input-file`, `--input-dir`, `--output-dir`, or the config overrides) can be a `gs://` URL. The code auto-detects: a normal path uses the local disk, a `gs://...` path uses the bucket. Everything else works exactly the same.

## 1. Install the cloud dependencies

GCS support needs **`fsspec`** and **`gcsfs`**. They are in `requirements.txt`, so a normal install covers them:

```powershell
pip install -r requirements.txt
```

(or install just these two: `pip install fsspec gcsfs`).

## 2. Authenticate (once per machine)

The code uses your **Application Default Credentials** — there are no keys in the code. Log in once:

```powershell
gcloud auth application-default login
```

This opens a browser; after you sign in, `gcsfs` picks up the credentials automatically. You also need read/write access to the bucket. If you skip this step you'll get a clear credentials error (not a hang).

## 3. Run against the bucket

```powershell
# One file in the bucket — output ProcessedData_<name>/ is written next to it
python src/Main.py --input-file gs://my-bucket/asms/asms_acme_01_lib_20260101.csv

# A folder of CSVs in the bucket
python src/Main.py --input-dir gs://my-bucket/asms/raw

# Read from the bucket, write results somewhere else
python src/Main.py --input-file gs://my-bucket/asms/run.csv --output-dir gs://my-bucket/asms/results
```

By default the config/reference files are still read from the **local repo**, so you don't have to copy them into the bucket. To pull them from the bucket instead, override the paths:

```powershell
python src/Main.py --input-file gs://my-bucket/asms/run.csv `
  --masterlists-dir gs://my-bucket/config/MasterLists `
  --providers-csv   gs://my-bucket/config/Providers.csv `
  --meta-csv        "gs://my-bucket/config/RawDataColumns.csv"
```

Mixing is fine — e.g. read from `gs://`, write locally, or vice-versa. Each path is detected independently.

## Notes & gotchas

- **`--input-file` vs `--input-dir`:** use `--input-file` for a single object and `--input-dir` for a folder. Passing a *file* to `--input-dir` is rejected with a clear error.
- **All step outputs go to the bucket** — Step 1–9 CSV/Parquet files, QC and post-QC `.log`/`.xlsx`, and any report CSVs all land under `ProcessedData_<name>/` in `--output-dir`.
- **GCS has no real folders** — a "folder" is just a prefix; it appears only once an object is written under it.

## Verify the GCS results

```powershell
gcloud storage ls "gs://my-bucket/asms/ProcessedData_<name>/"
```
