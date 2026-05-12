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

The pipeline takes two paths:

- `--input-dir` — folder containing `RawData/`, `MasterLists/`, and `Providers.csv`. Defaults to current directory.
- `--output-dir` — folder where `ProcessedData_<csv>/` is created. Defaults to `--input-dir`.

```powershell
# Run from the repo root with everything in place
python src/Main.py

# Specify input only (outputs land inside the same folder)
python src/Main.py --input-dir "D:\my\dataset"

# Separate input and output
python src/Main.py --input-dir "D:\my\dataset" --output-dir "D:\my\results"
```

Help text:

```powershell
python src/Main.py --help
```

## 4. Run only a subset of steps

Use `--start-from N` and `--end-at N` to control which steps execute. Step numbers are 1–9 (see [Readme.md](Readme.md#pipeline-steps) for what each step does). Skipped earlier steps are loaded from their saved output on disk.

```powershell
# Run only the Quality Check (step 0) and stop
python src/Main.py --end-at 0

# Run only steps 1 and 2
python src/Main.py --end-at 2

# Run only step 1
python src/Main.py --end-at 1

# Re-run from fingerprint extraction onward (steps 1-6 are loaded from disk)
python src/Main.py --start-from 7

# Run exactly one step (e.g. step 5)
python src/Main.py --start-from 5 --end-at 5
```

Defaults: `--start-from 0 --end-at 9` (run everything, including QC). Quality checks (step 0) run only when `--start-from 0`.
