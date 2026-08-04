# -*- coding: utf-8 -*-
"""Percentile- (rank-) normalized enrichment -> MultiBatch_20to22_modified4/.

modified4 = every column of MultiBatch_20to22_modified3/ PLUS a percentile-
normalized set (suffix _pctnorm).

Why: median normalization (modified3) only corrected the batch INTENSITY LEVEL,
not the DETECTION-RATE difference (sgcto_22 detects 72% of rows vs 17-36%), so a
residual batch effect remained in the across-batch calls. Percentile (rank)
normalization corrects the WHOLE distribution, so both axes are handled.

Normalization (rank -> common reference distribution, a.k.a. quantile norm):

    1. within each batch, convert every replicate value to its percentile
       (rank / n) among that batch's replicate values,
    2. map that percentile to the intensity at the same percentile of the POOLED
       (all-batch) replicate distribution -> the "reference".

Because the mapping is to reference INTENSITIES (not raw 0-1 percentiles), the
enrichment ratio stays meaningful. The 3000 floor is preserved: the floor block
sits below the pooled ~58%-floor level in every batch, so it maps back to ~3000.

Added columns (computed like the modified2/3 columns but on percentile-normalized
intensities):

    target_median_pctnorm
    nontarget_median_within_pctnorm, enrichment_within_pctnorm, pvalue_within_pctnorm, label_within_pctnorm
    nontarget_median_across_pctnorm, enrichment_across_pctnorm, pvalue_across_pctnorm, label_across_pctnorm

Notes
- Unlike median scaling (linear), percentile normalization is NON-linear, so the
  *_within_pctnorm ENRICHMENT differs slightly from the raw within enrichment.
  Within-batch p-values are unchanged (rank-based, monotonic within a batch).
- Trade-off: rank normalization is scale-free, so "enrichment >= 5" is no longer a
  literal 5x fold in raw intensity; treat it as a ranking/selection score.

Run from the repo root:  python multibatchcodes/compute_enrichment_percentile.py
"""

import os
import glob
import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.stats import mannwhitneyu, rankdata

SRC_DIR = "MultiBatch_20to22_modified3"
OUT_DIR = "MultiBatch_20to22_modified4"
REP = ["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]
FLOOR = 3000
ENR_THRESHOLD = 5.0
P_THRESHOLD = 0.05


def enrichment_and_p(reps, tmed, comp, batch, ndet):
    """Per-protein-median background enrichment + Mann-Whitney p, within & across."""
    n = len(tmed)
    nt_w = np.full(n, np.nan); nt_a = np.full(n, np.nan)
    p_w = np.full(n, np.nan); p_a = np.full(n, np.nan)
    idx_by_comp = defaultdict(list)
    for i, c in enumerate(comp):
        idx_by_comp[c].append(i)
    for c, idxs in idx_by_comp.items():
        idxs = np.array(idxs)
        vals = tmed[idxs]; bts = batch[idxs]
        for j, i in enumerate(idxs):
            others = np.delete(vals, j)
            if others.size:
                nt_a[i] = np.median(others)
                p_a[i] = 1.0 if ndet[i] == 0 else mannwhitneyu(
                    reps[i], others, alternative="two-sided").pvalue
            same = (bts == bts[j]); same[j] = False
            within = vals[same]
            if within.size:
                nt_w[i] = np.median(within)
                p_w[i] = 1.0 if ndet[i] == 0 else mannwhitneyu(
                    reps[i], within, alternative="two-sided").pvalue
    return nt_w, nt_a, p_w, p_a


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    paths = sorted(glob.glob(os.path.join(SRC_DIR, "*.csv")))
    print(f"Loading {len(paths)} files from {SRC_DIR}/ ...")
    dfs = {os.path.basename(p): pd.read_csv(p) for p in paths}
    big = pd.concat([df.assign(_src=name) for name, df in dfs.items()], ignore_index=True)

    reps = big[REP].to_numpy(dtype=float)
    batch = big["ASMS_BATCH_NAME"].to_numpy()
    comp = big["COMPOUND_ID"].to_numpy()
    ndet = big["n_detected"].to_numpy()

    # ---- percentile (rank) normalization to the pooled reference distribution ----
    pooled_sorted = np.sort(reps.ravel())
    grid = np.linspace(0.0, 1.0, len(pooled_sorted))
    reps_n = np.empty_like(reps)
    print("\npercentile-normalizing each batch to the pooled reference ...")
    for b in sorted(set(batch)):
        m = batch == b
        sub = reps[m].ravel()
        pct = rankdata(sub, method="average") / len(sub)     # within-batch percentile (0,1]
        reps_n[m] = np.interp(pct, grid, pooled_sorted).reshape(-1, 3)
        print(f"  {b}: floor now maps to ~{np.interp((sub == FLOOR).mean()/2, grid, pooled_sorted):,.0f}")

    tmed_n = np.median(reps_n, axis=1)

    print("\ncomputing percentile-normalized within/across enrichment + Mann-Whitney ...")
    nt_w, nt_a, p_w, p_a = enrichment_and_p(reps_n, tmed_n, comp, batch, ndet)
    enr_w = tmed_n / nt_w
    enr_a = tmed_n / nt_a
    lab_w = ((enr_w >= ENR_THRESHOLD) & (p_w <= P_THRESHOLD)).astype(int)
    lab_a = ((enr_a >= ENR_THRESHOLD) & (p_a <= P_THRESHOLD)).astype(int)

    add = {
        "target_median_pctnorm": tmed_n,
        "nontarget_median_within_pctnorm": nt_w,
        "enrichment_within_pctnorm": enr_w,
        "pvalue_within_pctnorm": p_w,
        "label_within_pctnorm": lab_w,
        "nontarget_median_across_pctnorm": nt_a,
        "enrichment_across_pctnorm": enr_a,
        "pvalue_across_pctnorm": p_a,
        "label_across_pctnorm": lab_a,
    }
    for k, v in add.items():
        big[k] = v

    for name, df in dfs.items():
        mask = (big["_src"] == name).to_numpy()
        for c in add:
            df[c] = big.loc[mask, c].to_numpy()
        df.to_csv(os.path.join(OUT_DIR, name), index=False)
    print(f"\nSaved {len(dfs)} files to {OUT_DIR}/")

    # quick before/after across-hit summary vs the earlier regimes
    batches = sorted(set(batch))
    print("\nacross hits per batch  (raw -> median-norm -> percentile-norm):")
    for b in batches:
        m = big["ASMS_BATCH_NAME"] == b
        print(f"  {b}: {int(big.loc[m,'label_across'].sum())} -> "
              f"{int(big.loc[m,'label_across_norm'].sum())} -> "
              f"{int(big.loc[m,'label_across_pctnorm'].sum())}")
    print(f"  ALL: {int(big['label_across'].sum())} -> "
          f"{int(big['label_across_norm'].sum())} -> {int(big['label_across_pctnorm'].sum())}")


if __name__ == "__main__":
    main()
