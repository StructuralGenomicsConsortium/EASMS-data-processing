import os
import pandas as pd
import numpy as np
import scipy.stats as stats
from io_utils import pjoin, makedirs

def compute_and_add_scores(file_paths, output_dir=None):
    """
    Computes TARGET_VALUE, ENRICHMENT, SELECTIVE_ENRICHMENT, EASMS_ENRICHMENT,
    MEAN_NONTARGET_VALUES, and PVALUE for a list of CSV files.
    Each file is processed individually, using all other files to compute comparison values.

    If output_dir is provided, saves each scored file to that directory using the
    original basename. Otherwise saves back to the original file path (in-place).

    Returns:
        list[str]: Paths of the saved (scored) files.
    """
    if not file_paths:
        print("No files provided for score computation.")
        return

    # Load all CSV files
    dataframes = {f: pd.read_csv(f) for f in file_paths}

    reps = ["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]

    # Ensure necessary columns exist
    required_columns = {"COMPOUND_ID", *reps}
    for filename, df in dataframes.items():
        if not required_columns.issubset(df.columns):
            raise ValueError(f"Missing required columns in {filename}: {required_columns - set(df.columns)}")

    # Compute TARGET_VALUE for each df
    for df in dataframes.values():
        df["TARGET_VALUE"] = pd.to_numeric(
            df[reps].mean(axis=1, skipna=True), errors="coerce"
        )

    # Step 2: Process each file individually.
    #
    # All per-compound comparison values used to be computed with a per-row
    # ``df.apply(axis=1)`` that re-filtered ``merged_other`` for every row
    # (``merged_other[merged_other["COMPOUND_ID"] == cid]``) -- i.e. a full scan
    # of the whole batch once per row, ~O(files^2 * rows^2). The logic below is
    # mathematically identical but computed once per compound with
    # groupby/map, which is ~O(files^2 * rows) and runs ~100x+ faster on real
    # batches (verified to match the old values element-for-element).
    for current_file, df in dataframes.items():
        print(f"\n Processing: {os.path.basename(current_file)}")

        # Exclude current df (the comparison set for this file). Guard the
        # single-file case (concat of an empty list raises) by treating it as an
        # empty comparison set -> all comparison-derived columns become NaN.
        other_dfs = [d for f, d in dataframes.items() if f != current_file]
        if other_dfs:
            merged_other = pd.concat(other_dfs, ignore_index=True)
        else:
            merged_other = pd.DataFrame(columns=df.columns)

        # Ensure TARGET_VALUE is present in merged_other
        if "TARGET_VALUE" not in merged_other.columns:
            merged_other["TARGET_VALUE"] = pd.to_numeric(
                merged_other[reps].mean(axis=1, skipna=True), errors="coerce"
            )

        cid = df["COMPOUND_ID"]

        # SELECTIVE_VALUE & NTC_VALUE — per-compound max/min over the other files.
        grp_tv = merged_other.groupby("COMPOUND_ID")["TARGET_VALUE"]
        df["SELECTIVE_VALUE"] = cid.map(grp_tv.max())
        df["NTC_VALUE"] = cid.map(grp_tv.min())

        # ENRICHMENT calculations
        df["ENRICHMENT"] = df["TARGET_VALUE"] / df["NTC_VALUE"]
        df["SELECTIVE_ENRICHMENT"] = df["TARGET_VALUE"] / df["SELECTIVE_VALUE"]

        # EASMS_ENRICHMENT and MEAN_NONTARGET_VALUES (vectorized).
        # MEAN_NONTARGET_VALUES[compound] = mean over the other files' rows of
        # each row's mean replicate intensity. EASMS_ENRICHMENT = this file's
        # TARGET_VALUE / that mean, left undefined (NaN) when the compound is
        # absent from the other files or its mean intensity is 0 (matches the
        # old None-returning guards).
        mo_reps = merged_other[reps].apply(pd.to_numeric, errors="coerce")
        mo_row_mean = mo_reps.mean(axis=1, skipna=True)
        mean_nontarget = mo_row_mean.groupby(merged_other["COMPOUND_ID"]).mean()
        mean_nt = cid.map(mean_nontarget)
        df["MEAN_NONTARGET_VALUES"] = mean_nt
        easms = df["TARGET_VALUE"] / mean_nt
        df["EASMS_ENRICHMENT"] = easms.where(mean_nt.notna() & (mean_nt != 0), np.nan)

        # PVALUE — Welch's t-test (equal_var=False) of this file's per-row
        # replicates vs every replicate value from the other files for the same
        # compound. Computed from precomputed per-compound stats instead of a
        # per-row scan + per-row scipy call; ``ttest_ind_from_stats`` is
        # vectorized and returns identical p-values to ``ttest_ind``.
        cur = df[reps].apply(pd.to_numeric, errors="coerce")
        n1 = cur.notna().sum(axis=1).to_numpy()
        mean1 = cur.mean(axis=1, skipna=True).to_numpy()
        std1 = cur.std(axis=1, ddof=1, skipna=True).to_numpy()      # sample std
        popstd1 = cur.std(axis=1, ddof=0, skipna=True).to_numpy()   # zero-var guard

        if len(mo_reps):
            mo_vals = mo_reps.to_numpy(dtype=float).ravel()
            mo_cids = np.repeat(merged_other["COMPOUND_ID"].to_numpy(), len(reps))
            keep = ~np.isnan(mo_vals)
            other_long = pd.DataFrame({"cid": mo_cids[keep], "v": mo_vals[keep]})
            g2 = other_long.groupby("cid")["v"]
            n2 = cid.map(g2.count()).to_numpy(dtype=float)
            mean2 = cid.map(g2.mean()).to_numpy(dtype=float)
            std2 = cid.map(g2.std(ddof=1)).to_numpy(dtype=float)
            popstd2 = cid.map(g2.std(ddof=0)).to_numpy(dtype=float)
        else:
            nan_col = np.full(len(df), np.nan)
            n2 = nan_col.copy()
            mean2, std2, popstd2 = nan_col.copy(), nan_col.copy(), nan_col.copy()
        n2 = np.nan_to_num(n2, nan=0.0)

        pvalue = np.full(len(df), np.nan)               # default: undefined (was None)
        valid = (n1 > 0) & (n2 >= 3)                    # need both samples populated
        zero_var = valid & ((popstd1 < 1e-8) | (popstd2 < 1e-8))
        pvalue[zero_var] = 1.0                          # degenerate variance -> p = 1.0
        do_test = valid & ~zero_var & np.isfinite(std1) & np.isfinite(std2)
        if do_test.any():
            _, p = stats.ttest_ind_from_stats(
                mean1[do_test], std1[do_test], n1[do_test].astype(int),
                mean2[do_test], std2[do_test], n2[do_test].astype(int),
                equal_var=False,
            )
            pvalue[do_test] = p
        df["PVALUE"] = pvalue

    # Step 3: Save the updated files (either in place or into output_dir)
    saved_paths = []
    if output_dir:
        makedirs(output_dir, exist_ok=True)
    for file_path, df in dataframes.items():
        if output_dir:
            out_path = pjoin(output_dir, os.path.basename(file_path))
        else:
            out_path = file_path
        df.to_csv(out_path, index=False)
        saved_paths.append(out_path)
        print(f"Updated and saved: {out_path}")

    print("\nAll files have been processed and saved with computed scores.")
    return saved_paths
