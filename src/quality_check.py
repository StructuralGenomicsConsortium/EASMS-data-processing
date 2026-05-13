"""
Pre-processing quality checks for raw input files.

Each check is a small function with signature:
    check(file_path, **context) -> (passed: bool, message: str)

`run_quality_checks` orchestrates them, writes a log file grouped by section,
and returns True if every check passed.
"""

import codecs
import os
import re
from datetime import datetime

import pandas as pd


# ---------- Configuration ----------

MAX_FILE_SIZE_GB = 10
MAX_FILE_SIZE_BYTES = MAX_FILE_SIZE_GB * 1024 ** 3

MIN_BATCH_NUMBER = 0
MAX_BATCH_NUMBER = 10000

# Filename format: asms_<provider>_<batch>_<library>_<date>.csv
# Library names may contain underscores; date (YYYYMMDD) anchors the tail.
FILENAME_RE = re.compile(
    r"^asms_(?P<provider>[a-z]+)_(?P<batch>\d{1,5})_(?P<library>.+)_(?P<date>\d{8})\.csv$",
    re.IGNORECASE,
)

# Allowed characters in the filename (alphanumeric, underscore, period, hyphen).
FILENAME_ALLOWED_RE = re.compile(r"^[A-Za-z0-9_.\-]+$")


# ---------- Helpers ----------

def _parse_filename(file_path):
    """Returns the regex match object or None."""
    return FILENAME_RE.match(os.path.basename(file_path))


def _load_providers(providers_csv_path):
    """Load valid provider acronyms from a CSV with an `acronym` column.

    Returns a list of lowercase strings, or None if the file is missing.
    """
    if not providers_csv_path or not os.path.exists(providers_csv_path):
        return None
    try:
        df = pd.read_csv(providers_csv_path)
    except Exception:
        return None
    if "acronym" not in df.columns:
        return []
    return [str(a).strip().lower() for a in df["acronym"].dropna()]


def _list_libraries(masterlist_dir):
    """Return library names (filename stems) from MasterLists/, excluding the
    MasterList_Information mapping file. None if directory is missing.
    """
    if not masterlist_dir or not os.path.isdir(masterlist_dir):
        return None
    libs = []
    for name in os.listdir(masterlist_dir):
        if name == "MasterList_Information.xlsx":
            continue
        stem, ext = os.path.splitext(name)
        if ext.lower() in (".xlsx", ".xls", ".csv"):
            libs.append(stem)
    return libs


def _load_meta_columns(meta_csv_path):
    """Load the valid column names from the ASMS Meta Data reference CSV.

    The reference file's header row lists the canonical column names; the
    second row holds data types. Whitespace is stripped and duplicates
    (e.g. an accidental trailing-space variant) are collapsed.

    Returns a list of strings, or None if the file is missing/unreadable.
    """
    if not meta_csv_path or not os.path.exists(meta_csv_path):
        return None
    try:
        df = pd.read_csv(meta_csv_path, nrows=0)
    except Exception:
        return None
    seen = []
    for col in df.columns:
        name = str(col).strip()
        if name and name not in seen:
            seen.append(name)
    return seen


# ---------- File-format checks ----------

def check_file_opens(file_path, **_):
    """File can be opened for reading."""
    try:
        with open(file_path, "rb") as f:
            f.read(1)
    except OSError as e:
        return False, f"could not open file: {e}"
    return True, "file opens for reading"


def check_is_csv(file_path, **_):
    """File has a .csv extension."""
    ext = os.path.splitext(file_path)[1].lower()
    if ext == ".csv":
        return True, "extension is '.csv'"
    return False, f"expected '.csv', got '{ext or '(no extension)'}'"


def check_file_not_empty(file_path, **_):
    """File size is greater than zero bytes."""
    try:
        size = os.path.getsize(file_path)
    except OSError as e:
        return False, f"could not stat file: {e}"
    if size == 0:
        return False, "file is empty (0 bytes)"
    return True, f"file size is {size:,} bytes ({size / (1024 ** 2):.2f} MB)"


def check_file_size_under_limit(file_path, **_):
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


def check_encoding_is_utf8(file_path, chunk_size=1024 * 1024, **_):
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


def check_csv_parseable(file_path, **_):
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


def check_columns_match_metadata(file_path, meta_columns=None, **_):
    """File's column names match the reference list in ASMS Meta Data.csv."""
    if meta_columns is None:
        return False, "metadata reference unavailable (ASMS Meta Data.csv not found)"
    if not meta_columns:
        return False, "metadata reference has no columns"

    try:
        df = pd.read_csv(file_path, nrows=0)
    except Exception as e:
        return False, f"could not read file columns: {e}"

    file_cols = [str(c).strip() for c in df.columns]
    valid = set(meta_columns)
    file_set = set(file_cols)

    missing = [c for c in meta_columns if c not in file_set]
    extra   = [c for c in file_cols    if c not in valid]

    if not missing and not extra:
        return True, "columns match ASMS Meta Data.csv reference"

    parts = ["columns do not match ASMS Meta Data.csv"]
    if missing:
        parts.append(f"missing required columns ({len(missing)}): {missing}")
    if extra:
        parts.append(f"extra columns not in reference ({len(extra)}): {extra}")
    return False, "; ".join(parts)


# ---------- Filename-format checks ----------

def check_filename_no_special_chars(file_path, **_):
    """Filename has only alphanumerics, underscores, periods, and hyphens."""
    name = os.path.basename(file_path)
    if FILENAME_ALLOWED_RE.match(name):
        return True, "no special characters or spaces in filename"
    bad = sorted(set(c for c in name if not re.match(r"[A-Za-z0-9_.\-]", c)))
    return False, f"filename contains disallowed character(s): {bad}"


def check_filename_starts_with_asms(file_path, **_):
    """Filename begins with 'asms_'."""
    name = os.path.basename(file_path)
    if name.lower().startswith("asms_"):
        return True, "filename starts with 'asms_'"
    return False, f"filename should start with 'asms_', got '{name[:10]}...'"


def check_filename_overall_format(file_path, **_):
    """Filename matches asms_<provider>_<batch>_<library>_<date>.csv."""
    match = _parse_filename(file_path)
    if match:
        return True, (
            f"parsed: provider='{match.group('provider')}', "
            f"batch='{match.group('batch')}', "
            f"library='{match.group('library')}', "
            f"date='{match.group('date')}'"
        )
    return False, (
        "filename does not match 'asms_<provider>_<batchN>_<library>_<YYYYMMDD>.csv'"
    )


def check_provider_acronym(file_path, providers=None, **_):
    """Provider acronym in the filename is in the registered list."""
    match = _parse_filename(file_path)
    if not match:
        return False, "filename did not parse; cannot extract provider"
    provider = match.group("provider").lower()

    if providers is None:
        return False, "providers list unavailable (Providers.csv not found)"
    if not providers:
        return False, "providers list is empty"
    if provider in providers:
        return True, f"provider '{provider}' is registered"
    return False, f"provider '{provider}' not in registered list: {providers}"


def check_batch_number_range(file_path, **_):
    """Batch number in the filename is between MIN_BATCH_NUMBER and MAX_BATCH_NUMBER."""
    match = _parse_filename(file_path)
    if not match:
        return False, "filename did not parse; cannot extract batch number"
    batch_str = match.group("batch")
    try:
        batch_int = int(batch_str)
    except ValueError:
        return False, f"batch number '{batch_str}' is not an integer"
    if MIN_BATCH_NUMBER <= batch_int <= MAX_BATCH_NUMBER:
        return True, f"batch number {batch_int} (from '{batch_str}') is in range [{MIN_BATCH_NUMBER}, {MAX_BATCH_NUMBER}]"
    return False, (
        f"batch number {batch_int} is outside [{MIN_BATCH_NUMBER}, {MAX_BATCH_NUMBER}]"
    )


def check_library_name(file_path, libraries=None, **_):
    """Library name in the filename matches a file in MasterLists/."""
    match = _parse_filename(file_path)
    if not match:
        return False, "filename did not parse; cannot extract library name"
    library = match.group("library")

    if libraries is None:
        return False, "libraries list unavailable (MasterLists/ not found)"
    if not libraries:
        return False, "no registered libraries (MasterLists/ is empty)"
    if library in libraries:
        return True, f"library '{library}' is registered"
    return False, f"library '{library}' not in registered list: {libraries}"


def check_date_valid_and_not_future(file_path, **_):
    """Date in the filename is a valid YYYYMMDD and not in the future."""
    match = _parse_filename(file_path)
    if not match:
        return False, "filename did not parse; cannot extract date"
    date_str = match.group("date")
    try:
        date_obj = datetime.strptime(date_str, "%Y%m%d").date()
    except ValueError as e:
        return False, f"date '{date_str}' is not a valid YYYYMMDD: {e}"
    today = datetime.now().date()
    if date_obj > today:
        return False, f"date {date_obj.isoformat()} is in the future (today is {today.isoformat()})"
    return True, f"date {date_obj.isoformat()} is valid and not in the future"


# ---------- Row content checks ----------

def check_no_duplicate_rows(file_path, output_dir=None, **_):
    """No two rows in the file are identical across all columns.

    Emits a WARN (not FAIL) when duplicates exist, because Step 3
    (anomaly_selection) will drop them anyway. Also writes
    `duplicate_rows_report.csv` into `output_dir` listing every row that
    is part of a duplicate group (all copies, not just the dropped ones).
    """
    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        return False, f"could not read file: {e}"

    dup_mask_first = df.duplicated(keep="first")
    n_dups = int(dup_mask_first.sum())
    if n_dups == 0:
        return True, f"no fully duplicate rows (checked {len(df):,} rows)"

    # Write a report of every row in any duplicate group (keep=False shows all copies)
    report_msg = ""
    if output_dir:
        try:
            os.makedirs(output_dir, exist_ok=True)
            report_df = df[df.duplicated(keep=False)].copy()
            report_df.insert(0, "FileLine", report_df.index + 2)  # +2: 1-index + header
            report_path = os.path.join(output_dir, "duplicate_rows_report.csv")
            report_df.to_csv(report_path, index=False)
            report_msg = f"; report saved to {report_path}"
        except Exception as e:
            report_msg = f"; failed to write report: {e}"

    # 1-indexed line numbers including the header row (so row 0 in df is line 2)
    dup_line_nums = [int(i) + 2 for i in dup_mask_first[dup_mask_first].index[:5]]
    suffix = f" (and {n_dups - 5} more)" if n_dups > 5 else ""
    return True, (
        f"found {n_dups:,} fully duplicate row(s) "
        f"(extras will be removed in Step 3). "
        f"First duplicates at file line(s): {dup_line_nums}{suffix}{report_msg}"
    ), "WARN"


# ---------- Section registry ----------

SECTIONS = [
    ("File Format Checks", [
        ("File opens without errors",                check_file_opens),
        ("File is a CSV (extension)",                check_is_csv),
        ("File is not empty",                        check_file_not_empty),
        (f"File size is under {MAX_FILE_SIZE_GB} GB", check_file_size_under_limit),
        ("File encoding is UTF-8",                   check_encoding_is_utf8),
        ("File is a CSV (parseable content)",        check_csv_parseable),
        ("Columns match ASMS Meta Data.csv reference", check_columns_match_metadata),
    ]),
    ("Filename Format Checks", [
        ("Filename has no special characters or spaces", check_filename_no_special_chars),
        ("Filename starts with 'asms_'",                 check_filename_starts_with_asms),
        ("Filename matches overall format",              check_filename_overall_format),
        ("Provider acronym is registered",               check_provider_acronym),
        ("Batch number is in valid range",               check_batch_number_range),
        ("Library name is registered",                   check_library_name),
        ("Date is valid YYYYMMDD and not in the future", check_date_valid_and_not_future),
    ]),
    ("Row Content Checks", [
        ("No fully duplicate rows", check_no_duplicate_rows),
    ]),
]


# ---------- Orchestrator ----------

SEPARATOR = "=" * 60


def run_quality_checks(file_path, log_dir, providers_csv=None, masterlist_dir=None,
                       meta_csv=None):
    """
    Run every check in SECTIONS against `file_path` and write a sectioned log.

    Args:
        file_path (str):       Path to the raw input file being validated.
        log_dir (str):         Directory where the log file is written.
        providers_csv (str):   Path to Providers.csv (acronym list). Optional.
        masterlist_dir (str):  Path to the MasterLists/ folder. Optional.
        meta_csv (str):        Path to ASMS Meta Data.csv (valid columns). Optional.

    Returns:
        tuple[bool, list[tuple[str, str]]]:
            (all_passed, failed_checks). `failed_checks` is a list of
            (description, message) pairs for every check that failed; empty
            when all_passed is True.
    """
    os.makedirs(log_dir, exist_ok=True)
    file_name = os.path.basename(file_path)
    base_name = os.path.splitext(file_name)[0]
    today = datetime.now().strftime("%Y%m%d")
    log_path = os.path.join(log_dir, f"QCaircheck{today}_{base_name}.log")

    # Pre-load context that several checks need. `output_dir` lets checks
    # write supplementary files (e.g. duplicate_rows_report.csv) next to the log.
    context = {
        "providers":    _load_providers(providers_csv),
        "libraries":    _list_libraries(masterlist_dir),
        "meta_columns": _load_meta_columns(meta_csv),
        "output_dir":   log_dir,
    }

    all_passed = True
    failed_checks = []   # collected for the caller to print a nice console message
    rows = []            # accumulated for the Excel report
    check_idx = 0
    generated_at = datetime.now().isoformat(timespec="seconds")
    with open(log_path, "w", encoding="utf-8") as log:
        log.write("Quality Check Log\n")
        log.write(f"File:      {file_name}\n")
        log.write(f"Generated: {generated_at}\n")
        log.write(SEPARATOR + "\n")

        for section_name, checks in SECTIONS:
            log.write("\n")
            log.write(f"{section_name}\n")
            log.write(SEPARATOR + "\n\n")

            for description, check_fn in checks:
                check_idx += 1
                try:
                    result = check_fn(file_path, **context)
                    # Checks can return 2-tuple (passed, msg) or 3-tuple
                    # (passed, msg, status) where status is "PASS", "FAIL", or "WARN".
                    if isinstance(result, tuple) and len(result) == 3:
                        passed, message, status = result
                    else:
                        passed, message = result
                        status = "PASS" if passed else "FAIL"
                except Exception as e:
                    passed, message, status = False, f"check raised an exception: {e}", "FAIL"
                log.write(f"Check {check_idx}: {description}\n")
                log.write(f"  Result : {status}\n")
                log.write(f"  Detail : {message}\n\n")
                rows.append({
                    "Section":  section_name,
                    "Check #":  check_idx,
                    "Criteria": description,
                    "Status":   status,
                    "Detail":   message,
                })
                if status == "FAIL":
                    all_passed = False
                    failed_checks.append((description, message))

        log.write(SEPARATOR + "\n")
        log.write(f"Overall : {'PASS' if all_passed else 'FAIL'}\n")

    # Write the Excel companion next to the .log
    excel_path = os.path.join(log_dir, f"QCaircheck{today}_{base_name}.xlsx")
    try:
        _write_excel_report(
            rows=rows,
            excel_path=excel_path,
            file_name=file_name,
            generated_at=generated_at,
            overall_status="PASS" if all_passed else "FAIL",
        )
    except Exception as e:
        # Excel write failure should not block QC; record in log
        with open(log_path, "a", encoding="utf-8") as log:
            log.write(f"\n(Note: failed to write Excel report: {e})\n")

    return all_passed, failed_checks


def _write_excel_report(rows, excel_path, file_name, generated_at, overall_status):
    """Write a color-coded Excel version of the QC log.

    Layout (one sheet, "QC Results"):
      Row 1: File:      <file_name>
      Row 2: Generated: <iso timestamp>
      Row 3: Overall:   <PASS/FAIL>
      Row 4: (blank)
      Row 5: Header (Section / Check # / Criteria / Status / Detail)
      Row 6+: Data, with each row's background color tied to its Status.
    """
    from openpyxl import Workbook
    from openpyxl.styles import PatternFill, Font, Alignment

    wb = Workbook()
    ws = wb.active
    ws.title = "QC Results"

    # Metadata block
    bold = Font(bold=True)
    ws["A1"] = "File:";      ws["B1"] = file_name
    ws["A2"] = "Generated:"; ws["B2"] = generated_at
    ws["A3"] = "Overall:";   ws["B3"] = overall_status
    for addr in ("A1", "A2", "A3"):
        ws[addr].font = bold

    # Header row
    header_row = 5
    headers = ["Section", "Check #", "Criteria", "Status", "Detail"]
    for col_idx, name in enumerate(headers, start=1):
        cell = ws.cell(row=header_row, column=col_idx, value=name)
        cell.font = bold

    # Data rows with color coding by status
    status_fills = {
        "PASS": PatternFill(start_color="C6EFCE", end_color="C6EFCE", fill_type="solid"),
        "FAIL": PatternFill(start_color="FFC7CE", end_color="FFC7CE", fill_type="solid"),
        "WARN": PatternFill(start_color="FFEB9C", end_color="FFEB9C", fill_type="solid"),
    }
    wrap = Alignment(wrap_text=True, vertical="top")
    for offset, row in enumerate(rows, start=1):
        excel_row = header_row + offset
        values = [row["Section"], row["Check #"], row["Criteria"], row["Status"], row["Detail"]]
        for col_idx, val in enumerate(values, start=1):
            cell = ws.cell(row=excel_row, column=col_idx, value=val)
            cell.fill = status_fills.get(row["Status"], status_fills["PASS"])
            cell.alignment = wrap

    # Column widths
    widths = {"A": 26, "B": 9, "C": 50, "D": 8, "E": 100}
    for col, w in widths.items():
        ws.column_dimensions[col].width = w

    # Freeze the header row so it stays visible when scrolling
    ws.freeze_panes = f"A{header_row + 1}"

    wb.save(excel_path)
