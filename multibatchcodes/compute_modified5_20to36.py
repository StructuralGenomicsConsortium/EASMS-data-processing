# -*- coding: utf-8 -*-
"""Full modified5-style columns for MultiBatch_20to36 -> MultiBatch_20to36_modified1/.

Adds, in ONE output folder, the whole analysis stack we built for the 3-batch set
(modified2..modified5), now computed over all 136 proteins / 17 batches:

    per-molecule : target_median, signal_flag        (n_detected already present)
    raw median   : within & across (nontarget_median, enrichment, pvalue, label)
    median-norm  : same, on per-batch median-scaled intensities  (_norm)
    percentile-norm: same, on rank/quantile-normalized intensities (_pctnorm)
    bead         : bead_signal, bead_detected, bead_ratio, bead_binder_flag_r5/r10

within-batch = the 7 other proteins in the same ASMS_BATCH_NAME.
across-batch = all 135 other proteins.

Speed shortcuts (mathematically exact, so we skip redundant Mann-Whitney passes):
  * median-norm is a per-batch LINEAR scale -> every within-batch ratio/rank is
    unchanged, so the *_within_norm columns are derived from the raw ones.
  * percentile-norm applies ONE monotonic map per batch; median commutes with a
    monotonic map, so within-batch RANKS (hence pvalue_within) are unchanged.
    Only enrichment_within_pctnorm differs (ratio of transformed medians), which
    needs medians only, no Mann-Whitney. Across-batch always mixes batches -> must
    be recomputed for both norms.

Run from the repo root:  python multibatchcodes/compute_modified5_20to36.py [--validate]
"""

import os
import sys
import glob
import functools
import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.stats import mannwhitneyu, rankdata

print = functools.partial(print, flush=True)  # live progress on long runs

SRC_DIR = "MultiBatch_20to36"
OUT_DIR = "MultiBatch_20to36_modified1"
RESULTS_DIR = "MultiBatchResults"
BEADS = os.path.join(RESULTS_DIR, "beads_clean.xlsx")
CACHE = os.path.join(RESULTS_DIR, "_20to36_newcols.pkl")  # computed cols, so a
#            locked output file never forces a ~35-min recompute (see --write-only)
REP = ["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]
FLOOR = 3000
STRONG_MIN = 1e6
MODERATE_MIN = 1e4
ENR_THRESHOLD = 5.0
P_THRESHOLD = 0.05
C5 = ["5Beads_V1P0113_1", "5Beads_V1P0113_2", "5Beads_V1P0113_3"]
C10 = ["10Beads_V1P0113_1", "10Beads_V1P0113_2", "10Beads_V1P0113_3"]


def signal_flags(reps):
    n_strong = (reps >= STRONG_MIN).sum(axis=1)
    n_mod = ((reps >= MODERATE_MIN) & (reps < STRONG_MIN)).sum(axis=1)
    flag = np.full(len(reps), "none", dtype=object)
    m = n_mod > 0
    flag[m] = [f"{k} moderate" for k in n_mod[m]]
    s = n_strong > 0
    flag[s] = [f"{k} strong" for k in n_strong[s]]
    return flag


def enrichment_and_p(reps, tmed, idx_by_comp, batch, ndet,
                     do_within=True, within_pvalue=True,
                     do_across=True, across_pvalue=True):
    """Per-protein-median background enrichment + Mann-Whitney, within & across.

    Flags let callers skip work known to be redundant (see module docstring).
    Returns nt_within, nt_across, p_within, p_across (NaN where not computed).
    """
    n = len(tmed)
    nt_w = np.full(n, np.nan); nt_a = np.full(n, np.nan)
    p_w = np.full(n, np.nan); p_a = np.full(n, np.nan)
    for c, idxs in idx_by_comp.items():
        idxs = np.asarray(idxs)
        vals = tmed[idxs]
        bts = batch[idxs]
        for j, i in enumerate(idxs):
            if do_across:
                others = np.delete(vals, j)
                if others.size:
                    nt_a[i] = np.median(others)
                    if across_pvalue:
                        p_a[i] = 1.0 if ndet[i] == 0 else mannwhitneyu(
                            reps[i], others, alternative="two-sided").pvalue
            if do_within:
                same = (bts == bts[j]); same[j] = False
                within = vals[same]
                if within.size:
                    nt_w[i] = np.median(within)
                    if within_pvalue:
                        p_w[i] = 1.0 if ndet[i] == 0 else mannwhitneyu(
                            reps[i], within, alternative="two-sided").pvalue
    return nt_w, nt_a, p_w, p_a


def labels(enr, pval):
    return ((enr >= ENR_THRESHOLD) & (pval <= P_THRESHOLD)).astype(int)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    paths = sorted(glob.glob(os.path.join(SRC_DIR, "*.csv")))
    if "--write-only" in sys.argv:   # resume writing from cache (no recompute)
        print(f"write-only: loading cached columns from {CACHE}")
        report_write(paths, write_outputs(paths, pd.read_pickle(CACHE)))
        return
    print(f"Loading {len(paths)} files (lightweight) from {SRC_DIR}/ ...")
    parts = []
    for p in paths:
        d = pd.read_csv(p, usecols=["COMPOUND_ID", "ASMS_BATCH_NAME", "n_detected"] + REP,
                        low_memory=False)
        d["_src"] = os.path.basename(p)
        parts.append(d)
    big = pd.concat(parts, ignore_index=True)
    reps = big[REP].to_numpy(dtype=float)
    # A handful of rows have a NaN replicate. Treat a missing read as the floor
    # (not detected) for ALL derived metrics -- the original rep columns in the
    # output are left untouched. Critical: without this, a NaN flows through the
    # percentile transform (np.sort -> NaN at top -> np.interp) and then the
    # 135-protein np.median background, turning ~all of nontarget_median_across_*
    # into NaN.
    n_nan_rows = int(np.isnan(reps).any(axis=1).sum())
    reps = np.where(np.isnan(reps), FLOOR, reps)
    print(f"imputed {n_nan_rows} rows with a NaN replicate -> floor {FLOOR}")
    batch = big["ASMS_BATCH_NAME"].to_numpy()
    comp = big["COMPOUND_ID"].to_numpy()
    ndet = big["n_detected"].to_numpy()
    tmed = np.median(reps, axis=1)
    n = len(big)
    print(f"{n:,} rows | {big['COMPOUND_ID'].nunique():,} compounds | {len(set(batch))} batches")

    idx_by_comp = defaultdict(list)
    for i, c in enumerate(comp):
        idx_by_comp[c].append(i)

    # ---------- 1. raw median enrichment (within + across) ----------
    print("raw median enrichment (within + across Mann-Whitney) ...")
    nt_w, nt_a, p_w, p_a = enrichment_and_p(reps, tmed, idx_by_comp, batch, ndet)
    enr_w, enr_a = tmed / nt_w, tmed / nt_a

    # ---------- 2. median normalization ----------
    print("median-normalized (across only; within derived) ...")
    grand = np.median(reps[reps > FLOOR])
    frow = np.empty(n)
    for b in set(batch):
        m = batch == b
        frow[m] = grand / np.median(reps[m][reps[m] > FLOOR])
    reps_n = reps * frow[:, None]
    tmed_n = tmed * frow
    _, nt_a_n, _, p_a_n = enrichment_and_p(reps_n, tmed_n, idx_by_comp, batch, ndet,
                                           do_within=False)
    enr_a_n = tmed_n / nt_a_n
    nt_w_n = nt_w * frow                 # within background scales linearly
    enr_w_n = enr_w                      # ratio unchanged
    p_w_n = p_w                          # rank unchanged

    # ---------- 3. percentile normalization ----------
    print("percentile-normalized (rank -> pooled reference) ...")
    pooled_sorted = np.sort(reps.ravel())
    grid = np.linspace(0.0, 1.0, len(pooled_sorted))
    reps_p = np.empty_like(reps)
    for b in set(batch):
        m = batch == b
        sub = reps[m].ravel()
        pct = rankdata(sub, method="average") / len(sub)
        reps_p[m] = np.interp(pct, grid, pooled_sorted).reshape(-1, 3)
    tmed_p = np.median(reps_p, axis=1)
    # within needs medians (no MWU: pvalue == raw within); across needs full MWU
    nt_w_p, nt_a_p, _, p_a_p = enrichment_and_p(reps_p, tmed_p, idx_by_comp, batch, ndet,
                                                within_pvalue=False)
    enr_w_p, enr_a_p = tmed_p / nt_w_p, tmed_p / nt_a_p
    p_w_p = p_w                          # within ranks preserved by the per-batch monotonic map

    # ---------- 4. bead specificity ----------
    print("bead specificity ...")
    b = pd.read_excel(BEADS)
    b["bead_signal"] = np.maximum(b[C5].median(axis=1), b[C10].median(axis=1))
    bead = (b[["SGC ID for Component", "bead_signal"]]
            .rename(columns={"SGC ID for Component": "COMPOUND_ID"})
            .dropna(subset=["COMPOUND_ID"]).drop_duplicates("COMPOUND_ID"))
    bmap = dict(zip(bead["COMPOUND_ID"], bead["bead_signal"]))
    bead_signal = np.array([bmap.get(c, np.nan) for c in comp])
    bead_detected = (bead_signal > FLOOR).astype(int)
    bead_ratio = tmed / bead_signal
    bbr5 = ((bead_signal > FLOOR) & (bead_ratio < 5)).astype(int)
    bbr10 = ((bead_signal > FLOOR) & (bead_ratio < 10)).astype(int)

    # ---------- assemble new columns ----------
    cols = {
        "target_median": tmed, "signal_flag": signal_flags(reps),
        "nontarget_median_within": nt_w, "enrichment_within": enr_w,
        "pvalue_within": p_w, "label_within": labels(enr_w, p_w),
        "nontarget_median_across": nt_a, "enrichment_across": enr_a,
        "pvalue_across": p_a, "label_across": labels(enr_a, p_a),
        "norm_scale_factor": frow, "target_median_norm": tmed_n,
        "nontarget_median_within_norm": nt_w_n, "enrichment_within_norm": enr_w_n,
        "pvalue_within_norm": p_w_n, "label_within_norm": labels(enr_w_n, p_w_n),
        "nontarget_median_across_norm": nt_a_n, "enrichment_across_norm": enr_a_n,
        "pvalue_across_norm": p_a_n, "label_across_norm": labels(enr_a_n, p_a_n),
        "target_median_pctnorm": tmed_p,
        "nontarget_median_within_pctnorm": nt_w_p, "enrichment_within_pctnorm": enr_w_p,
        "pvalue_within_pctnorm": p_w_p, "label_within_pctnorm": labels(enr_w_p, p_w_p),
        "nontarget_median_across_pctnorm": nt_a_p, "enrichment_across_pctnorm": enr_a_p,
        "pvalue_across_pctnorm": p_a_p, "label_across_pctnorm": labels(enr_a_p, p_a_p),
        "bead_signal": bead_signal, "bead_detected": bead_detected,
        "bead_ratio": bead_ratio, "bead_binder_flag_r5": bbr5, "bead_binder_flag_r10": bbr10,
    }
    newdf = pd.DataFrame(cols)

    if "--validate" in sys.argv:
        validate(reps, tmed, comp, batch, ndet, nt_a, p_a)

    # NaN sanity on the across backgrounds (this was the pctnorm bug)
    for c in ["nontarget_median_across", "nontarget_median_across_norm",
              "nontarget_median_across_pctnorm"]:
        print(f"  NaN in {c}: {newdf[c].isna().sum():,} ({newdf[c].isna().mean():.1%})")
    for lab in ["label_within", "label_across", "label_across_norm", "label_across_pctnorm"]:
        print(f"  {lab}=1: {int(newdf[lab].sum()):,}")

    # cache the computed columns BEFORE writing, so a single locked output file
    # never forces a recompute -- resume with --write-only after closing it.
    newdf["_src"] = big["_src"].to_numpy()
    newdf.to_pickle(CACHE)
    print(f"cached computed columns to {CACHE}")
    report_write(paths, write_outputs(paths, newdf))


def write_outputs(paths, newdf):
    """Re-read each source file, append the computed columns, write to OUT_DIR.

    Resilient to locked output files (open in Excel): skips them and returns the
    list of names that could not be written.
    """
    src = newdf["_src"].to_numpy()
    newcols = [c for c in newdf.columns if c != "_src"]
    locked = []
    for p in paths:
        name = os.path.basename(p)
        add = newdf.loc[src == name, newcols].reset_index(drop=True)
        try:
            full = pd.read_csv(p, low_memory=False)
            for c in newcols:
                full[c] = add[c].to_numpy()
            full.to_csv(os.path.join(OUT_DIR, name), index=False)
        except PermissionError:
            locked.append(name)
    return locked


def report_write(paths, locked):
    print(f"\nwrote {len(paths) - len(locked)}/{len(paths)} files to {OUT_DIR}/")
    if locked:
        print(f"[LOCKED] {len(locked)} not written (open in Excel?): {locked}")
        print("  close them, then run:  "
              "python multibatchcodes/compute_modified5_20to36.py --write-only")
    print(f"  bead_binder_flag_r10=1 among label_within hits: "
          f"{int(((newdf['label_within']==1)&(newdf['bead_binder_flag_r10']==1)).sum()):,}")


def validate(reps, tmed, comp, batch, ndet, nt_a, p_a, n_check=200, seed=0):
    """Direct recompute of across nontarget-median + MWU on a random sample."""
    rng = np.random.default_rng(seed)
    idx_by_comp = defaultdict(list)
    for i, c in enumerate(comp):
        idx_by_comp[c].append(i)
    idx = rng.choice(len(reps), size=n_check, replace=False)
    max_nt = max_p = 0.0
    for i in idx:
        peers = np.array([k for k in idx_by_comp[comp[i]] if k != i])
        if not peers.size:
            continue
        med = np.median(tmed[peers])
        max_nt = max(max_nt, abs(med - nt_a[i]))
        if ndet[i] > 0:
            pv = mannwhitneyu(reps[i], tmed[peers], alternative="two-sided").pvalue
            max_p = max(max_p, abs(pv - p_a[i]))
    print(f"[validate] {n_check} rows | max abs-diff nt_across={max_nt:.2e} p_across={max_p:.2e}")


if __name__ == "__main__":
    main()
