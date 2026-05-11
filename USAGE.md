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

Pass your dataset root directory (the folder containing `RawData/` and `MasterLists/`):

```powershell
python src/Main.py --path "D:\0000-UHN\03-DataAndCodes\Data\ASMS\EASMS_05Feb2026_batch2all"
```

If you omit `--path`, the pipeline defaults to the parent of the current working directory:

```powershell
python src/Main.py
```

Help text:

```powershell
python src/Main.py --help
```
