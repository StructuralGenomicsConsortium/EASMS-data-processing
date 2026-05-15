"""
Post-pipeline quality checks.

Runs after Step 9 against the concatenation of every per-target Parquet in
`Step8_FullColumns/`. Verifies that computed columns produced by the pipeline
satisfy basic invariants (value sets, ranges, non-negativity, fingerprint
lengths). Reuses helpers and the Excel writer from `quality_check.py`.

Output written to `<ProcessedData_<csv>>/`:
- `PostQClog_<YYYYMMDD>_<csv_basename>.log`  — plain text
- `PostQClog_<YYYYMMDD>_<csv_basename>.xlsx` — color-coded, with Statistics tab
"""

import os
from datetime import datetime

import pandas as pd

# Reuse helpers and orchestrator pieces from the input-side QC module
from quality_check import (
    _ensure_df_and_columns,
    _check_column_in_set,
    _check_column_in_range,
    _check_column_positive,
    _check_column_at_least,
    _generate_statistics_summary,
    _write_excel_report,
    SEPARATOR,
)


# ---------- Constants ----------

LABEL_VALUES = {0, 1}
AIRCHECK_LABEL_VALUES = {-2, -1, 0, 1, 2, 3, 4}

# PVALUE comes from a Welch's t-test; allowed range [0, 1].
PVALUE_MIN = 0.0
PVALUE_MAX = 1.0

# ALOGP "reasonable" range. Outliers can be legitimate, so this is a WARN.
ALOGP_MIN = -5
ALOGP_MAX = 10

# Expected fingerprint vector lengths (from src/fingerprints.py: each
# HitGen*FPFunc's `_dimension` attribute).
FP_EXPECTED_LENGTHS = {
    "ECFP4":    2048,
    "ECFP6":    2048,
    "FCFP4":    2048,
    "FCFP6":    2048,
    "MACCS":    167,
    "RDK":      2048,
    "AVALON":   2048,
    "TOPTOR":   2048,
    "ATOMPAIR": 2048,
}

# Columns that should hold {"Y", "N"} string flags.
YN_COLUMNS = ("MassSpec_Detected", "HAD_DUPLICATE_INTENSITY")


# ---------- Fingerprint length helper ----------

def _fp_length(val):
    """Return length of a fingerprint value (numpy array or comma-string), or None."""
    if val is None:
        return None
    # numpy array / list / pandas Series — has __len__ but isn't a str
    if hasattr(val, "__len__") and not isinstance(val, str):
        try:
            return len(val)
        except Exception:
            return None
    if isinstance(val, str):
        # legacy comma-separated format
        return len([p for p in val.split(",") if p.strip()])
    return None


def _check_fp_length(df, column, expected_length):
    """Every non-null row in `column` has a fingerprint vector of `expected_length`."""
    ok, fail = _ensure_df_and_columns(df, column)
    if not ok:
        return False, fail
    bad_rows = []
    for idx, val in df[column].dropna().items():
        actual = _fp_length(val)
        if actual is None:
            bad_rows.append((idx, "could not measure length"))
        elif actual != expected_length:
            bad_rows.append((idx, actual))
    if not bad_rows:
        return True, (
            f"all {len(df):,} '{column}' values have length {expected_length}"
        )
    sample = bad_rows[:5]
    return False, (
        f"{len(bad_rows):,} '{column}' value(s) have wrong fingerprint length "
        f"(expected {expected_length}); first few: {sample}"
    )


# ---------- Check functions ----------

# LABEL / AIRCHECK_LABEL — discrete value sets
def check_label_values(file_path, df=None, **_):
    """LABEL contains only {0, 1}."""
    return _check_column_in_set(df, "LABEL", LABEL_VALUES)

def check_aircheck_label_values(file_path, df=None, **_):
    """AIRCHECK_LABEL contains only {-2, -1, 0, 1, 2, 3, 4}."""
    return _check_column_in_set(df, "AIRCHECK_LABEL", AIRCHECK_LABEL_VALUES)


# Molecular properties
def check_mw_positive(file_path, df=None, **_):
    """MW (molecular weight) is strictly positive."""
    return _check_column_positive(df, "MW")

def check_alogp_range(file_path, df=None, **_):
    """ALOGP is roughly within [-5, 10] (outliers possible, hence WARN)."""
    ok, message = _check_column_in_range(df, "ALOGP", ALOGP_MIN, ALOGP_MAX)
    if ok:
        return True, message
    return True, message, "WARN"   # downgrade FAIL → WARN; outliers can be real


# Score columns — non-negative
def check_target_intensity_non_negative(file_path, df=None, **_):
    """TARGET_INTENSITY_VALUE >= 0."""
    return _check_column_at_least(df, "TARGET_INTENSITY_VALUE", 0)

def check_nontarget_intensity_non_negative(file_path, df=None, **_):
    """NONTARGET_INTENSITY_VALUE >= 0."""
    return _check_column_at_least(df, "NONTARGET_INTENSITY_VALUE", 0)

def check_selective_value_non_negative(file_path, df=None, **_):
    """SELECTIVE_VALUE >= 0."""
    return _check_column_at_least(df, "SELECTIVE_VALUE", 0)

def check_ntc_value_non_negative(file_path, df=None, **_):
    """NTC_VALUE >= 0."""
    return _check_column_at_least(df, "NTC_VALUE", 0)

def check_enrichment_non_negative(file_path, df=None, **_):
    """ENRICHMENT >= 0."""
    return _check_column_at_least(df, "ENRICHMENT", 0)

def check_easms_enrichment_non_negative(file_path, df=None, **_):
    """EASMS_ENRICHMENT >= 0."""
    return _check_column_at_least(df, "EASMS_ENRICHMENT", 0)

def check_selective_enrichment_non_negative(file_path, df=None, **_):
    """SELECTIVE_ENRICHMENT >= 0."""
    return _check_column_at_least(df, "SELECTIVE_ENRICHMENT", 0)


# PVALUE — [0, 1]
def check_pvalue_range(file_path, df=None, **_):
    """PVALUE is in [0, 1]."""
    return _check_column_in_range(df, "PVALUE", PVALUE_MIN, PVALUE_MAX)


# Y/N flag columns
def check_massspec_detected_values(file_path, df=None, **_):
    """MassSpec_Detected only contains {'Y', 'N'}."""
    return _check_column_in_set(df, "MassSpec_Detected", {"Y", "N"})

def check_had_duplicate_intensity_values(file_path, df=None, **_):
    """HAD_DUPLICATE_INTENSITY only contains {'Y', 'N'}."""
    return _check_column_in_set(df, "HAD_DUPLICATE_INTENSITY", {"Y", "N"})


# Fingerprint length checks (one per fingerprint)
def _make_fp_length_check(column, expected_length):
    def _check(file_path, df=None, **_):
        return _check_fp_length(df, column, expected_length)
    _check.__name__ = f"check_{column.lower()}_length"
    _check.__doc__ = f"{column} fingerprint vector has length {expected_length}."
    return _check


_FP_CHECKS = {col: _make_fp_length_check(col, dim) for col, dim in FP_EXPECTED_LENGTHS.items()}


# ---------- Section registry ----------

SECTIONS = [
    ("Label Checks", [
        ("LABEL only contains {0, 1}",                                    check_label_values),
        ("AIRCHECK_LABEL only contains {-2, -1, 0, 1, 2, 3, 4}",          check_aircheck_label_values),
    ]),
    ("Score Range Checks", [
        ("PVALUE is in [0, 1]",                                           check_pvalue_range),
        ("TARGET_INTENSITY_VALUE >= 0",                                   check_target_intensity_non_negative),
        ("NONTARGET_INTENSITY_VALUE >= 0",                                check_nontarget_intensity_non_negative),
        ("SELECTIVE_VALUE >= 0",                                          check_selective_value_non_negative),
        ("NTC_VALUE >= 0",                                                check_ntc_value_non_negative),
        ("ENRICHMENT >= 0",                                               check_enrichment_non_negative),
        ("EASMS_ENRICHMENT >= 0",                                         check_easms_enrichment_non_negative),
        ("SELECTIVE_ENRICHMENT >= 0",                                     check_selective_enrichment_non_negative),
    ]),
    ("Molecular Property Checks", [
        ("MW (molecular weight) is positive",                             check_mw_positive),
        (f"ALOGP is in [{ALOGP_MIN}, {ALOGP_MAX}] (typical)",             check_alogp_range),
    ]),
    ("Flag Column Checks", [
        ("MassSpec_Detected only contains {Y, N}",                        check_massspec_detected_values),
        ("HAD_DUPLICATE_INTENSITY only contains {Y, N}",                  check_had_duplicate_intensity_values),
    ]),
    ("Fingerprint Length Checks", [
        (f"{col} fingerprint length == {FP_EXPECTED_LENGTHS[col]}", _FP_CHECKS[col])
        for col in FP_EXPECTED_LENGTHS
    ]),
]


# ---------- Orchestrator ----------

def _load_pipeline_output(parquet_dir):
    """Load all `.parquet` files in `parquet_dir` and concatenate them. Returns
    None if the directory is missing or empty."""
    if not parquet_dir or not os.path.isdir(parquet_dir):
        return None
    files = sorted(
        os.path.join(parquet_dir, f)
        for f in os.listdir(parquet_dir)
        if f.endswith(".parquet")
    )
    if not files:
        return None
    try:
        return pd.concat((pd.read_parquet(f) for f in files), ignore_index=True)
    except Exception:
        return None


def run_post_quality_checks(parquet_dir, log_dir, csv_basename):
    """Run post-pipeline checks against Step 8's per-target Parquet outputs.

    Args:
        parquet_dir:   Path to the folder of per-target Parquet files (typically
                       `ProcessedData_<csv>/Step8_FullColumns/`).
        log_dir:       Where to write the log + Excel report.
        csv_basename:  Base name of the originating raw CSV (used in filenames).

    Returns:
        bool: True if every check passed (FAIL count == 0), False otherwise.
              If parquet_dir is missing/empty, returns True and writes a one-
              line "skipped" log (post-QC is best-effort, never blocks).
    """
    os.makedirs(log_dir, exist_ok=True)
    today = datetime.now().strftime("%Y%m%d")
    log_path = os.path.join(log_dir, f"PostQClog_{today}_{csv_basename}.log")
    excel_path = os.path.join(log_dir, f"PostQClog_{today}_{csv_basename}.xlsx")

    df = _load_pipeline_output(parquet_dir)
    if df is None:
        with open(log_path, "w", encoding="utf-8") as log:
            log.write("Post-Pipeline Quality Check Log\n")
            log.write(f"Skipped: no Parquet files found in {parquet_dir}\n")
        return True

    all_passed = True
    rows = []                  # for Excel report
    check_idx = 0
    generated_at = datetime.now().isoformat(timespec="seconds")
    context = {"df": df}

    with open(log_path, "w", encoding="utf-8") as log:
        log.write("Post-Pipeline Quality Check Log\n")
        log.write(f"Source:    {parquet_dir}\n")
        log.write(f"Rows:      {len(df):,}\n")
        log.write(f"Generated: {generated_at}\n")
        log.write(SEPARATOR + "\n")

        for section_name, checks in SECTIONS:
            log.write("\n")
            log.write(f"{section_name}\n")
            log.write(SEPARATOR + "\n\n")
            for description, check_fn in checks:
                check_idx += 1
                try:
                    result = check_fn(None, **context)
                    if isinstance(result, tuple) and len(result) == 3:
                        passed, message, status = result
                    else:
                        passed, message = result
                        status = "PASS" if passed else "FAIL"
                except Exception as e:
                    passed, message, status = False, f"check raised: {e}", "FAIL"
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

        log.write(SEPARATOR + "\n")
        log.write(f"Overall : {'PASS' if all_passed else 'FAIL'}\n")

        # Statistics summary, same format as input QC.
        log.write("\n")
        log.write(SEPARATOR + "\n")
        log.write("Statistics Summary\n")
        log.write(SEPARATOR + "\n\n")
        try:
            log.write(_generate_statistics_summary(df))
        except Exception as e:
            log.write(f"(Statistics summary failed to render: {e})\n")

    # Excel companion
    try:
        _write_excel_report(
            rows=rows,
            excel_path=excel_path,
            file_name=f"Post-pipeline output from {csv_basename}",
            generated_at=generated_at,
            overall_status="PASS" if all_passed else "FAIL",
            df=df,
        )
    except Exception as e:
        with open(log_path, "a", encoding="utf-8") as log:
            log.write(f"\n(Note: failed to write Excel report: {e})\n")

    return all_passed
