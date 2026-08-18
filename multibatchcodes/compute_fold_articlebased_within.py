# -*- coding: utf-8 -*-
"""Add the article's within-batch enrichment fold as `fold_articlebased_within`.

Formula (equation 1 of the E-ASMS paper):

                          A_POI  x  N
    Enrichment Fold = ---------------------
                            sum(A_i)

    A_POI  peak intensity of the compound on the protein of interest (this row)
    N      number of BACKGROUND proteins from the SAME BATCH (the paper: N = 7)
    A_i    peak intensity of the SAME compound on the i-th background protein
           of that same batch

So the background is strictly WITHIN batch: for a given compound, the other
proteins that were run together with it in the same `ASMS_BATCH_NAME`. The
expression is just  A_POI / mean(background A_i).

Peak intensity A is taken as the mean of POS_INT_REP1/2/3 (nan-safe). Because
the same aggregation is used in the numerator and the denominator, using the
sum of the replicates instead of the mean would give exactly the same fold.

In this data every batch holds exactly 8 proteins and every
(batch, compound) group has exactly 8 rows, so N = 7 everywhere -- identical to
the paper. The code does not assume that: N is computed per group as
(number of proteins with a MEASURED intensity - 1), and groups with N = 0 get
NaN. The 13 rows of MultiBatch_20to36 whose replicates are all missing get NaN
themselves and drop out of the background of their own group, which there means
N = 6 for the other 7 proteins of that (batch, compound) group.

The column is appended IN PLACE to every *.csv of each folder in DIRS (each
folder is computed independently from its own files, so a folder that holds a
subset of the data still gets a correct within-batch background). Files are
written via a temp file + atomic replace, so an interrupted run cannot leave a
half-written CSV.

Run from the repo root:
    python multibatchcodes/compute_fold_articlebased_within.py
    python multibatchcodes/compute_fold_articlebased_within.py --validate
"""

import os
import sys
import glob
import warnings

import numpy as np
import pandas as pd

DIRS = [
    "MultiBatch_20to22",
    "MultiBatch_20to22_modified1",
    "MultiBatch_20to22_modified2",
    "MultiBatch_20to22_modified3",
    "MultiBatch_20to22_modified4",
    "MultiBatch_20to22_modified5",
    "MultiBatch_20to36",
    "MultiBatch_20to36_modified1",
]

REP = ["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]
KEYS = ["ASMS_BATCH_NAME", "COMPOUND_ID"]
NEW_COL = "fold_articlebased_within"


def peak_intensity(df):
    """A = mean of the three replicate peak intensities.

    nan-safe: a row with some replicates missing uses the ones it has; a row
    with all three missing (13 such rows exist in MultiBatch_20to36) gets NaN
    and is then excluded from every background it would take part in.
    """
    r = df[REP].to_numpy(dtype=float)
    with warnings.catch_warnings():
        # all-NaN rows -> "Mean of empty slice"; NaN is the intended answer
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmean(r, axis=1)


def batch_background_stats(paths):
    """Pass 1: per (batch, compound) total intensity S and MEASURED count K.

    K counts only rows with a defined A (`count` skips NaN), so a protein whose
    intensity is missing is not counted as one of the N background proteins --
    it contributes neither to N nor to sum(A_i). Otherwise N would include a
    protein that adds nothing to the denominator and the fold would be inflated.

    Reads only the columns needed for the maths, so the whole folder never has
    to be held in memory at full width.
    """
    parts = []
    for p in paths:
        d = pd.read_csv(p, usecols=KEYS + REP, low_memory=False)
        parts.append(pd.DataFrame({
            "ASMS_BATCH_NAME": d["ASMS_BATCH_NAME"],
            "COMPOUND_ID": d["COMPOUND_ID"],
            "_A": peak_intensity(d),
        }))
    big = pd.concat(parts, ignore_index=True)
    stats = big.groupby(KEYS, sort=False)["_A"].agg(S="sum", K="count").reset_index()
    return stats, big


def fold_from_stats(A, S, K):
    """A_POI * N / sum(A_i)  with  N = K - 1  and  sum(A_i) = S - A_POI."""
    N = K - 1.0
    denom = S - A
    with np.errstate(divide="ignore", invalid="ignore"):
        fold = A * N / denom
    fold = np.where(N <= 0, np.nan, fold)          # no background proteins
    fold = np.where(denom <= 0, np.nan, fold)      # empty/degenerate background
    return fold


def process_dir(src):
    paths = sorted(glob.glob(os.path.join(src, "*.csv")))
    if not paths:
        print(f"[skip] {src}/ has no CSV files")
        return None

    stats, big = batch_background_stats(paths)
    sizes = stats["K"].value_counts().sort_index().to_dict()
    print(f"\n== {src}/  {len(paths)} files | {len(big):,} rows | "
          f"{stats['ASMS_BATCH_NAME'].nunique()} batches | "
          f"{len(stats):,} (batch, compound) groups | group sizes {sizes}")

    n_nan = 0
    n_rows = 0
    for p in paths:
        df = pd.read_csv(p, low_memory=False)
        A = peak_intensity(df)
        m = df[KEYS].merge(stats, on=KEYS, how="left")   # left merge keeps row order
        fold = fold_from_stats(A, m["S"].to_numpy(), m["K"].to_numpy())

        df[NEW_COL] = fold                               # overwrites if re-run
        tmp = p + ".tmp"
        df.to_csv(tmp, index=False)
        os.replace(tmp, p)

        n_rows += len(df)
        n_nan += int(np.isnan(fold).sum())

    print(f"   wrote {NEW_COL} to {len(paths)} files | {n_rows:,} rows | "
          f"{n_nan:,} NaN")
    return stats


def validate(src, n=300, seed=0):
    """Recompute the fold row-by-row, straight from the CSVs, and compare.

    Deliberately naive: for each sampled row it re-reads the compound's rows in
    the same batch from the other protein files and applies equation 1 literally.
    """
    paths = sorted(glob.glob(os.path.join(src, "*.csv")))
    frames = {}
    for p in paths:
        d = pd.read_csv(p, usecols=KEYS + REP + [NEW_COL], low_memory=False)
        d["_A"] = peak_intensity(d)
        frames[p] = d

    rng = np.random.default_rng(seed)
    all_rows = [(p, i) for p, d in frames.items() for i in range(len(d))]
    picks = rng.choice(len(all_rows), size=min(n, len(all_rows)), replace=False)

    max_err = 0.0
    checked = 0
    for j in picks:
        p, i = all_rows[j]
        row = frames[p].iloc[i]
        if np.isnan(row["_A"]):          # unmeasured target -> fold is NaN by design
            continue
        bg = []
        for q, d in frames.items():
            if q == p:
                continue
            hit = d[(d["ASMS_BATCH_NAME"] == row["ASMS_BATCH_NAME"]) &
                    (d["COMPOUND_ID"] == row["COMPOUND_ID"])]
            bg.extend(a for a in hit["_A"].tolist() if not np.isnan(a))
        if not bg:
            continue
        expected = row["_A"] * len(bg) / sum(bg)
        got = row[NEW_COL]
        max_err = max(max_err, abs(expected - got) / max(abs(expected), 1e-12))
        checked += 1
    print(f"[validate] {src}: {checked} rows re-derived from scratch | "
          f"max relative error = {max_err:.2e}")
    return max_err


def main():
    for src in DIRS:
        if not os.path.isdir(src):
            print(f"[skip] {src}/ not found")
            continue
        process_dir(src)

    if "--validate" in sys.argv:
        validate("MultiBatch_20to22")
        validate("MultiBatch_20to36", n=120)


if __name__ == "__main__":
    main()
