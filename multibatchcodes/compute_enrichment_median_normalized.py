# -*- coding: utf-8 -*-
"""Median-normalized enrichment -> MultiBatch_20to22_modified3/.

modified3 = every column of MultiBatch_20to22_modified2/ PLUS a normalized set.

Batch effect: sgcto_22 runs ~2x hotter than sgcto_20/21, which inflates the
across-batch enrichment for that batch. We remove it with per-batch median
scaling, then recompute the across-batch metrics on the normalized intensities.

Normalization (per-batch multiplicative scaling to a common reference):

    batch_med_b  = median of DETECTED replicate values (>3000) in batch b
    grand_med    = median of DETECTED replicate values over ALL batches
    norm_scale_factor(b) = grand_med / batch_med_b        # >1 cold batch, <1 hot
    normalized replicate = raw replicate * norm_scale_factor(b)

After scaling, every batch's detected median = grand_med, so the across-batch
background is on a common scale.

Added columns (all suffixed _norm; computed exactly like the modified2 columns
but on normalized intensities):

    norm_scale_factor
    target_median_norm              = median of the 3 normalized reps
    nontarget_median_within_norm, enrichment_within_norm, pvalue_within_norm, label_within_norm
    nontarget_median_across_norm, enrichment_across_norm, pvalue_across_norm, label_across_norm

NOTE: within-batch scaling cancels in every ratio/rank, so the *_within_norm
columns are mathematically identical to the modified2 within columns (kept for a
clean parallel set; the run asserts the equality). Only the across-batch columns
change. Plots of intensities-after-normalization and before/after across hits are
saved to MultiBatchResults/.

Run from the repo root:  python multibatchcodes/compute_enrichment_median_normalized.py
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.stats import mannwhitneyu

SRC_DIR = "MultiBatch_20to22_modified2"
OUT_DIR = "MultiBatch_20to22_modified3"
RESULTS_DIR = "MultiBatchResults"
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
        vals = tmed[idxs]
        bts = batch[idxs]
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
    os.makedirs(RESULTS_DIR, exist_ok=True)
    paths = sorted(glob.glob(os.path.join(SRC_DIR, "*.csv")))
    print(f"Loading {len(paths)} files from {SRC_DIR}/ ...")
    dfs = {os.path.basename(p): pd.read_csv(p) for p in paths}

    parts = [df.assign(_src=name) for name, df in dfs.items()]
    big = pd.concat(parts, ignore_index=True)
    reps = big[REP].to_numpy(dtype=float)
    batch = big["ASMS_BATCH_NAME"].to_numpy()
    comp = big["COMPOUND_ID"].to_numpy()
    ndet = big["n_detected"].to_numpy()

    # ---- per-batch median scaling factors (from DETECTED replicate values) ----
    detected = reps[reps > FLOOR]
    grand_med = np.median(detected)
    factor = {}
    print("\nnorm_scale_factor per batch (grand_med / batch_med):")
    for b in sorted(set(batch)):
        bm = np.median(reps[batch == b][reps[batch == b] > FLOOR])
        factor[b] = grand_med / bm
        print(f"  {b}: batch_med={bm:,.0f}  factor={factor[b]:.3f}")
    print(f"  grand_med (reference) = {grand_med:,.0f}")

    frow = np.array([factor[b] for b in batch])            # per-row scale
    reps_n = reps * frow[:, None]                          # normalized replicates
    tmed_n = np.median(reps_n, axis=1)                     # target_median_norm

    print("\ncomputing normalized within/across enrichment + Mann-Whitney ...")
    nt_w, nt_a, p_w, p_a = enrichment_and_p(reps_n, tmed_n, comp, batch, ndet)

    enr_w = tmed_n / nt_w
    enr_a = tmed_n / nt_a
    lab_w = ((enr_w >= ENR_THRESHOLD) & (p_w <= P_THRESHOLD)).astype(int)
    lab_a = ((enr_a >= ENR_THRESHOLD) & (p_a <= P_THRESHOLD)).astype(int)

    add = {
        "norm_scale_factor": frow,
        "target_median_norm": tmed_n,
        "nontarget_median_within_norm": nt_w,
        "enrichment_within_norm": enr_w,
        "pvalue_within_norm": p_w,
        "label_within_norm": lab_w,
        "nontarget_median_across_norm": nt_a,
        "enrichment_across_norm": enr_a,
        "pvalue_across_norm": p_a,
        "label_across_norm": lab_a,
    }
    for k, v in add.items():
        big[k] = v

    # sanity: within-batch metrics must be unchanged by per-batch scaling
    ok = np.allclose(big["enrichment_within_norm"].fillna(-1),
                     big["enrichment_within"].fillna(-1), rtol=1e-9)
    print(f"[check] enrichment_within_norm == enrichment_within : {ok}")
    print(f"[check] within labels identical: "
          f"{int((big['label_within_norm'] == big['label_within']).all())}")

    # ---- write modified3 ----
    newcols = list(add.keys())
    for name, df in dfs.items():
        mask = (big["_src"] == name).to_numpy()
        for c in newcols:
            df[c] = big.loc[mask, c].to_numpy()
        df.to_csv(os.path.join(OUT_DIR, name), index=False)
    print(f"\nSaved {len(dfs)} files to {OUT_DIR}/")

    make_plots(big)
    return big


def make_plots(big):
    batches = sorted(big["ASMS_BATCH_NAME"].unique())

    # 1. intensities before vs after normalization (per-batch box, detected only)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for ax, col, title in [(axes[0], "target_median", "Before normalization"),
                           (axes[1], "target_median_norm", "After median normalization")]:
        data = [big.loc[(big["ASMS_BATCH_NAME"] == b) & (big[col] > FLOOR * 0.5), col]
                for b in batches]
        ax.boxplot(data, tick_labels=batches, showfliers=False, whis=(5, 95))
        meds = [d.median() for d in data]
        for i, m in enumerate(meds, 1):
            ax.text(i, m, f"{m:,.0f}", ha="center", va="bottom", fontsize=8)
        ax.set_yscale("log")
        ax.set_title(title)
        ax.set_ylabel("target intensity (detected, log)")
    fig.suptitle("Per-batch intensity — median normalization aligns the batches")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "modified3_intensity_before_after.png"), dpi=150)

    # 2. across-batch hits before vs after normalization
    before = big.groupby("ASMS_BATCH_NAME")["label_across"].sum().reindex(batches)
    after = big.groupby("ASMS_BATCH_NAME")["label_across_norm"].sum().reindex(batches)
    x = np.arange(len(batches)); w = 0.38
    fig, ax = plt.subplots(figsize=(8, 4.5))
    b1 = ax.bar(x - w/2, before.values, w, label="label_across (raw)")
    b2 = ax.bar(x + w/2, after.values, w, label="label_across_norm")
    ax.bar_label(b1, fmt="%d", fontsize=8); ax.bar_label(b2, fmt="%d", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(batches)
    ax.set_ylabel("across-batch hits")
    ax.set_title("Across-batch calls before vs after normalization\n"
                 "the sgcto_22 inflation collapses")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, "modified3_across_hits_before_after.png"), dpi=150)

    print("\ntotal across hits  raw:", int(big["label_across"].sum()),
          "| normalized:", int(big["label_across_norm"].sum()))
    print("across hits per batch (raw -> norm):")
    for b in batches:
        print(f"  {b}: {int(before[b])} -> {int(after[b])}")


if __name__ == "__main__":
    main()
