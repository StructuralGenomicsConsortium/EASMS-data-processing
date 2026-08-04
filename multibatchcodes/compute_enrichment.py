# -*- coding: utf-8 -*-
"""Compute per-compound enrichment + t-test across ALL proteins.

Input  -> MultiBatch_20to22/  (batches 20, 21, 22)
Output -> MultiBatch_20to22_modified1/

For every protein file in MultiBatch_20to22/, and every row (one compound vs that
protein), we compare the compound's signal on THIS protein (the "target") to its
signal on every OTHER protein (the "non-target" background, a leave-one-file-out
group keyed by COMPOUND_ID):

    target reps      = POS_INT_REP1/2/3 on this protein            (n1 = 3)
    non-target reps  = POS_INT_REP1/2/3 of the SAME COMPOUND_ID
                       pooled over all other protein files          (n2 = 3*(K-1))

    enrichment_shay  = mean(target reps) / mean(non-target reps)
    pvalue_shay      = two-sample Student's t-test (EQUAL variance, i.e. pooled;
                       NOT Welch) between the target reps and non-target reps
    label_shay       = 1 if enrichment_shay >= 5 AND pvalue_shay <= 0.05 else 0

This generalises an earlier calculation that used only 7 other proteins as the
background; here the background is every other protein present for the compound.

Each output file in MultiBatch_20to22_modified1/ is the original
MultiBatch_20to22/ file with these three extra columns appended (all other
columns untouched).

The t-test is vectorised over all ~290k rows using per-compound leave-one-out
sums (sum and sum-of-squares), which is mathematically identical to calling
scipy.stats.ttest_ind(..., equal_var=True) row by row but far faster. The
`--validate` flag checks that equivalence on a random sample.

Run from the repo root:  python multibatchcodes/compute_enrichment.py
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
from scipy import stats

SRC_DIR = "MultiBatch_20to22"
OUT_DIR = "MultiBatch_20to22_modified1"
REP = ["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]
ENR_THRESHOLD = 5.0
P_THRESHOLD = 0.05
N1 = len(REP)  # 3 target replicates, no missing values in this data


def compute_enrichment_and_pvalue(big):
    """Add enrichment_shay / pvalue_shay / label_shay columns to `big` in place.

    `big` must have COMPOUND_ID and the three REP columns. Rows are grouped by
    COMPOUND_ID; for each row the non-target group is that compound's replicates
    from all OTHER rows (i.e. other protein files).
    """
    r = big[REP].to_numpy(dtype=float)
    sA = r.sum(axis=1)                 # per-row sum of the 3 target reps
    qA = (r ** 2).sum(axis=1)          # per-row sum of squares

    big["_sA"] = sA
    big["_qA"] = qA
    grp = big.groupby("COMPOUND_ID")
    S = grp["_sA"].transform("sum").to_numpy()   # sum over ALL reps of this compound
    Q = grp["_qA"].transform("sum").to_numpy()   # sum of squares over ALL reps
    K = grp["_sA"].transform("size").to_numpy()  # number of protein files for compound
    big.drop(columns=["_sA", "_qA"], inplace=True)

    n1 = N1
    n2 = n1 * (K - 1)                  # non-target replicate count (leave-one-out)

    # non-target = everything for this compound minus this row's target reps
    sB = S - sA
    qB = Q - qA

    mean1 = sA / n1                    # target mean intensity
    mean2 = sB / n2                    # non-target (background) mean intensity
    enrichment = mean1 / mean2

    # Student's two-sample t-test (pooled / equal-variance):
    #   s^2 = ((n1-1)s1^2 + (n2-1)s2^2) / (n1+n2-2)
    #   t   = (mean1-mean2) / sqrt(s^2 (1/n1 + 1/n2)),  df = n1+n2-2
    with np.errstate(divide="ignore", invalid="ignore"):
        s1sq = (qA - sA ** 2 / n1) / (n1 - 1)
        s2sq = (qB - sB ** 2 / n2) / (n2 - 1)
        pooled = ((n1 - 1) * s1sq + (n2 - 1) * s2sq) / (n1 + n2 - 2)
        se = np.sqrt(pooled * (1.0 / n1 + 1.0 / n2))
        tstat = (mean1 - mean2) / se
        df = n1 + n2 - 2
        pval = 2.0 * stats.t.sf(np.abs(tstat), df)
    # se == 0 means both groups are constant (e.g. all pinned at the 3000 floor):
    # equal means -> no difference (p = 1); different means -> p ~ 0. This matches
    # scipy's limiting behaviour and keeps p in [0, 1] instead of NaN.
    degenerate = se == 0
    pval = np.where(degenerate, np.where(mean1 == mean2, 1.0, 0.0), pval)

    big["enrichment_shay"] = enrichment
    big["pvalue_shay"] = pval
    big["label_shay"] = (
        (enrichment >= ENR_THRESHOLD) & (pval <= P_THRESHOLD)
    ).astype(int)
    return big


def validate(big, n=400, seed=0):
    """Cross-check the vectorised t-test against scipy on a random sample.

    `big` has a RangeIndex (ignore_index concat), so index == row position. For
    each sampled row we build its non-target group by removing exactly that one
    row's position from its compound's rows (true leave-one-out), then compare
    scipy.ttest_ind(equal_var=True) to the vectorised pvalue_shay.
    """
    import warnings

    rng = np.random.default_rng(seed)
    r = big[REP].to_numpy(dtype=float)
    cid_col = big["COMPOUND_ID"].to_numpy()
    pos_by_compound = {cid: sub.index.to_numpy() for cid, sub in big.groupby("COMPOUND_ID")}

    idx = rng.choice(len(big), size=min(n, len(big)), replace=False)
    max_abs_err = 0.0
    compared = 0
    for i in idx:
        positions = pos_by_compound[cid_col[i]]
        others = positions[positions != i]      # drop exactly this one row
        target = r[i]
        nontarget = r[others].ravel()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, p_scipy = stats.ttest_ind(target, nontarget, equal_var=True)
        if np.isfinite(p_scipy):                 # skip all-constant (scipy -> NaN)
            max_abs_err = max(max_abs_err, abs(p_scipy - big.iat[i, big.columns.get_loc("pvalue_shay")]))
            compared += 1
    print(f"[validate] {compared}/{len(idx)} finite rows | "
          f"max |p_vec - p_scipy| = {max_abs_err:.2e}")
    return max_abs_err


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    paths = sorted(glob.glob(os.path.join(SRC_DIR, "*.csv")))
    print(f"Loading {len(paths)} files from {SRC_DIR}/ ...")
    dfs = {os.path.basename(p): pd.read_csv(p) for p in paths}

    # lightweight combined view (COMPOUND_ID + reps + source tag) for the maths
    parts = []
    for name, df in dfs.items():
        part = df[["COMPOUND_ID"] + REP].copy()
        part["_src"] = name
        parts.append(part)
    big = pd.concat(parts, ignore_index=True)
    print(f"{len(big):,} total rows | {big['COMPOUND_ID'].nunique():,} unique compounds")

    compute_enrichment_and_pvalue(big)

    if "--validate" in sys.argv:
        validate(big)

    # write each original file back out with the 3 new columns appended
    newcols = ["enrichment_shay", "pvalue_shay", "label_shay"]
    total_hits = 0
    for name, df in dfs.items():
        mask = (big["_src"] == name).to_numpy()
        for c in newcols:
            df[c] = big.loc[mask, c].to_numpy()
        hits = int(df["label_shay"].sum())
        total_hits += hits
        df.to_csv(os.path.join(OUT_DIR, name), index=False)
        print(f"  {name:16} {len(df):>6,} rows | {hits:>4} label_shay=1")

    print(f"\nSaved {len(dfs)} files to {OUT_DIR}/  |  "
          f"{total_hits:,} total label_shay=1 "
          f"(enrichment>={ENR_THRESHOLD:g} & p<={P_THRESHOLD:g})")


if __name__ == "__main__":
    main()
