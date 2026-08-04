# -*- coding: utf-8 -*-
"""Median-based enrichment, within- and across-batch -> MultiBatch_20to22_modified2/.

A median-based re-do of the enrichment calculation (the mean-based version lives
in MultiBatch_20to22_modified1/). Median is more robust to the 3000 floor and to a
single noisy replicate/protein.

Batches come from ASMS_BATCH_NAME (3 batches x 8 proteins here):
    within-batch  -> the 7 OTHER proteins in the same batch
    across-batch  -> all 23 OTHER proteins

Per molecule row (one compound on one protein):

    target_median            = median of the row's 3 replicates (POS_INT_REP1/2/3)
    n_detected               = # replicates above the 3000 floor (0..3)
    signal_flag              = most-interesting category present in the 3 reps:
                               "3/2/1 strong"  (>=1e6, incl. very abundant),
                               "3/2/1 moderate" (1e4..1e6), else "none" (weak/floor)

Non-target background = the SAME compound's per-protein medians in the other
proteins (one median per protein -> "median of each protein first"), then:

    nontarget_median_within  = median of the 7 same-batch other-protein medians
    enrichment_within        = target_median / nontarget_median_within
    pvalue_within            = Mann-Whitney U, target's 3 reps vs those 7 medians
    label_within             = 1 if enrichment_within >= 5 and pvalue_within <= 0.05

    nontarget_median_across  = median of the 23 other-protein medians (all batches)
    enrichment_across        = target_median / nontarget_median_across
    pvalue_across            = Mann-Whitney U, target's 3 reps vs those 23 medians
    label_across             = 1 if enrichment_across >= 5 and pvalue_across <= 0.05

Notes
- The p-value uses the SAME background sample as the effect size (the per-protein
  medians), so significance and magnitude are consistent. Within-batch has only 7
  background proteins, so pvalue_within is coarse (min two-sided p ~ 0.017).
- 3000 (floor) values are kept in the medians (a direct median swap of the mean
  version); a compound undetected off-target gives a ~3000 denominator, so
  enrichment then measures rise above background/floor.

Each output file = the original MultiBatch_20to22 file with the columns above
appended (nothing else changed). Run from the repo root:

    python multibatchcodes/compute_enrichment_median.py [--validate]
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.stats import mannwhitneyu

SRC_DIR = "MultiBatch_20to22"
OUT_DIR = "MultiBatch_20to22_modified2"
REP = ["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]
FLOOR = 3000
STRONG_MIN = 1e6      # strong or very abundant
MODERATE_MIN = 1e4    # moderate lower bound
ENR_THRESHOLD = 5.0
P_THRESHOLD = 0.05


def signal_flags(reps):
    """Per-row flag string from the 3 replicate categories (strong outranks moderate)."""
    n_strong = (reps >= STRONG_MIN).sum(axis=1)
    n_mod = ((reps >= MODERATE_MIN) & (reps < STRONG_MIN)).sum(axis=1)
    flag = np.full(len(reps), "none", dtype=object)
    m = n_mod > 0
    flag[m] = [f"{k} moderate" for k in n_mod[m]]
    s = n_strong > 0
    flag[s] = [f"{k} strong" for k in n_strong[s]]     # applied last -> overrides moderate
    return flag


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    paths = sorted(glob.glob(os.path.join(SRC_DIR, "*.csv")))
    print(f"Loading {len(paths)} files from {SRC_DIR}/ ...")
    dfs = {os.path.basename(p): pd.read_csv(p) for p in paths}

    # combined view for the maths
    parts = []
    for name, df in dfs.items():
        part = df[["COMPOUND_ID", "ASMS_BATCH_NAME"] + REP].copy()
        part["_src"] = name
        parts.append(part)
    big = pd.concat(parts, ignore_index=True)

    reps = big[REP].to_numpy(dtype=float)
    target_median = np.median(reps, axis=1)
    n_detected = (reps > FLOOR).sum(axis=1)
    flag = signal_flags(reps)
    comp = big["COMPOUND_ID"].to_numpy()
    batch = big["ASMS_BATCH_NAME"].to_numpy()
    n = len(big)
    print(f"{n:,} rows | {big['COMPOUND_ID'].nunique():,} compounds | "
          f"{len(set(batch))} batches")

    nt_within = np.full(n, np.nan)
    nt_across = np.full(n, np.nan)
    p_within = np.full(n, np.nan)
    p_across = np.full(n, np.nan)

    idx_by_comp = defaultdict(list)
    for i, c in enumerate(comp):
        idx_by_comp[c].append(i)

    def mwu(t, bg):
        # target 3 reps vs background per-protein medians; two-sided
        if len(bg) == 0:
            return np.nan
        if (t == FLOOR).all():          # target never detected -> not significant
            return 1.0
        return mannwhitneyu(t, bg, alternative="two-sided").pvalue

    for c, idxs in idx_by_comp.items():
        idxs = np.array(idxs)
        vals = target_median[idxs]      # per-protein medians for this compound
        bts = batch[idxs]
        for j, i in enumerate(idxs):
            # across-batch: every other protein
            others = np.delete(vals, j)
            if others.size:
                nt_across[i] = np.median(others)
                p_across[i] = mwu(reps[i], others)
            # within-batch: other proteins sharing this batch
            same = (bts == bts[j])
            same[j] = False
            within = vals[same]
            if within.size:
                nt_within[i] = np.median(within)
                p_within[i] = mwu(reps[i], within)

    enr_within = target_median / nt_within
    enr_across = target_median / nt_across
    label_within = ((enr_within >= ENR_THRESHOLD) & (p_within <= P_THRESHOLD)).astype(int)
    label_across = ((enr_across >= ENR_THRESHOLD) & (p_across <= P_THRESHOLD)).astype(int)

    out = {
        "target_median": target_median,
        "n_detected": n_detected,
        "signal_flag": flag,
        "nontarget_median_within": nt_within,
        "enrichment_within": enr_within,
        "pvalue_within": p_within,
        "label_within": label_within,
        "nontarget_median_across": nt_across,
        "enrichment_across": enr_across,
        "pvalue_across": p_across,
        "label_across": label_across,
    }
    for k, v in out.items():
        big[k] = v

    if "--validate" in sys.argv:
        validate(big, reps, target_median, comp, batch)

    newcols = list(out.keys())
    tw = ta = 0
    for name, df in dfs.items():
        mask = (big["_src"] == name).to_numpy()
        for c in newcols:
            df[c] = big.loc[mask, c].to_numpy()
        tw += int(df["label_within"].sum())
        ta += int(df["label_across"].sum())
        df.to_csv(os.path.join(OUT_DIR, name), index=False)
        print(f"  {name:16} {len(df):>6,} rows | within={int(df['label_within'].sum()):>4} "
              f"across={int(df['label_across'].sum()):>4}")

    print(f"\nSaved {len(dfs)} files to {OUT_DIR}/  |  "
          f"label_within=1: {tw:,}  label_across=1: {ta:,}")
    # signal_flag breakdown
    print("\nsignal_flag counts (all rows):")
    print(big["signal_flag"].value_counts().to_string())


def validate(big, reps, target_median, comp, batch, n_check=300, seed=0):
    """Recompute nt medians + MWU p-values for a random sample the slow, direct way."""
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(big), size=min(n_check, len(big)), replace=False)
    max_nt = max_p = 0.0
    for i in idx:
        c, b = comp[i], batch[i]
        oth = np.where((comp == c) & (np.arange(len(big)) != i))[0]
        med_oth = target_median[oth]
        # across
        if med_oth.size:
            max_nt = max(max_nt, abs(np.median(med_oth) - big["nontarget_median_across"].iloc[i]))
            if not (reps[i] == FLOOR).all():
                p = mannwhitneyu(reps[i], med_oth, alternative="two-sided").pvalue
                max_p = max(max_p, abs(p - big["pvalue_across"].iloc[i]))
    print(f"[validate] {len(idx)} rows | max abs-diff nontarget_median_across = {max_nt:.2e} | "
          f"max abs-diff pvalue_across = {max_p:.2e}")


if __name__ == "__main__":
    main()
