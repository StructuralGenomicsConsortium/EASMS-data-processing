"""
Pre-processing quality checks for raw input files.

Each check is a small function that returns (passed: bool, message: str).
`run_quality_checks` orchestrates them, writes a log file, and returns
True if every check passed.
"""

import codecs
import os
from datetime import datetime

import pandas as pd


# File size limits
MAX_FILE_SIZE_GB = 10
MAX_FILE_SIZE_BYTES = MAX_FILE_SIZE_GB * 1024 ** 3


# ---------- Individual checks ----------

def check_file_opens(file_path):
    """File can be opened for reading."""
    try:
        with open(file_path, "rb") as f:
            f.read(1)
    except OSError as e:
        return False, f"could not open file: {e}"
    return True, "file opens for reading"


def check_is_csv(file_path):
    """File has a .csv extension."""
    ext = os.path.splitext(file_path)[1].lower()
    if ext == ".csv":
        return True, "extension is '.csv'"
    return False, f"expected '.csv', got '{ext or '(no extension)'}'"


def check_file_not_empty(file_path):
    """File size is greater than zero bytes."""
    try:
        size = os.path.getsize(file_path)
    except OSError as e:
        return False, f"could not stat file: {e}"
    if size == 0:
        return False, "file is empty (0 bytes)"
    return True, f"file size is {size:,} bytes ({size / (1024 ** 2):.2f} MB)"


def check_file_size_under_limit(file_path):
    """File size is below the configured limit (default 10 GB)."""
    try:
        size = os.path.getsize(file_path)
    except OSError as e:
        return False, f"could not stat file: {e}"
    if size > MAX_FILE_SIZE_BYTES:
        return False, (
            f"file size {size / (1024 ** 3):.2f} GB exceeds the "
            f"{MAX_FILE_SIZE_GB} GB limit"
        )
    return True, (
        f"file size {size / (1024 ** 3):.4f} GB is within the "
        f"{MAX_FILE_SIZE_GB} GB limit"
    )


def check_encoding_is_utf8(file_path, chunk_size=1024 * 1024):
    """File contents decode cleanly as UTF-8 (reads the whole file in chunks)."""
    decoder = codecs.getincrementaldecoder("utf-8")()
    try:
        with open(file_path, "rb") as f:
            while True:
                chunk = f.read(chunk_size)
                if not chunk:
                    decoder.decode(b"", final=True)
                    break
                decoder.decode(chunk, final=False)
    except UnicodeDecodeError as e:
        return False, f"not UTF-8: {e}"
    except OSError as e:
        return False, f"read error: {e}"
    return True, "decodes as UTF-8"


def check_csv_parseable(file_path):
    """File parses as CSV via pandas; also report rows, columns, and column names."""
    try:
        df = pd.read_csv(file_path, nrows=5)
    except Exception as e:
        return False, f"pandas could not parse as CSV: {e}"
    n_cols = len(df.columns)
    if n_cols < 1:
        return False, "parsed but has zero columns"

    # Count total data rows (line count minus header). Streamed so it stays
    # cheap even for multi-GB files.
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            line_count = sum(1 for _ in f)
        rows_msg = f"{max(line_count - 1, 0):,} rows"
    except OSError:
        rows_msg = "row count unavailable"

    col_list = ", ".join(df.columns)
    return True, (
        f"parsed as CSV with {n_cols} columns and {rows_msg}. "
        f"Columns: {col_list}"
    )


# Register all checks here. Order is preserved in the log; cheap and
# foundational checks come first so failures surface early.
CHECKS = [
    ("File opens without errors",                check_file_opens),
    ("File is a CSV (extension)",                check_is_csv),
    ("File is not empty",                        check_file_not_empty),
    (f"File size is under {MAX_FILE_SIZE_GB} GB", check_file_size_under_limit),
    ("File encoding is UTF-8",                   check_encoding_is_utf8),
    ("File is a CSV (parseable content)",        check_csv_parseable),
]


# ---------- Orchestrator ----------

def run_quality_checks(file_path, log_dir):
    """
    Run every check in CHECKS against `file_path` and write a log file.

    Args:
        file_path (str): Path to the raw input file being validated.
        log_dir (str):   Directory where the log file should be saved.
                         Created if it does not exist.

    Returns:
        bool: True if every check passed, False otherwise.
    """
    os.makedirs(log_dir, exist_ok=True)
    file_name = os.path.basename(file_path)
    base_name = os.path.splitext(file_name)[0]
    log_path = os.path.join(log_dir, f"QCLog-{base_name}.log")

    all_passed = True
    with open(log_path, "w", encoding="utf-8") as log:
        log.write("Quality Check Log\n")
        log.write(f"File:      {file_name}\n")
        log.write(f"Generated: {datetime.now().isoformat(timespec='seconds')}\n")
        log.write("=" * 60 + "\n\n")

        for i, (description, check_fn) in enumerate(CHECKS, start=1):
            try:
                passed, message = check_fn(file_path)
            except Exception as e:
                passed, message = False, f"check raised an exception: {e}"
            status = "PASS" if passed else "FAIL"
            log.write(f"Check {i}: {description}\n")
            log.write(f"  Result : {status}\n")
            log.write(f"  Detail : {message}\n\n")
            if not passed:
                all_passed = False

        log.write("=" * 60 + "\n")
        log.write(f"Overall : {'PASS' if all_passed else 'FAIL'}\n")

    return all_passed
