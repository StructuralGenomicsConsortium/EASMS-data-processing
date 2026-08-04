# -*- coding: utf-8 -*-
"""Bead-only replicate consistency analysis.

Bead-only (no protein) control data from MultiBatchResults/beads_clean.xlsx.
Each compound is measured under two conditions, each with three replicates:

    5-bead  : 5Beads_V1P0113_1 / _2 / _3
    10-bead : 10Beads_V1P0113_1 / _2 / _3

Question: how self-consistent are these replicates? We check, in order,
    (1) do the three 5-bead reps agree with each other,
    (2) do the three 10-bead reps agree with each other,
    (3) do the 5-bead and 10-bead conditions agree,
plus per-column and pooled intensity histograms.

3000 is a hard floor ("not detected"). Rows where a whole condition sits at the
floor carry no agreement signal, so they are excluded from the agreement metrics
(three identical 3000's would otherwise report perfect, meaningless agreement).

Each '# %%' section is a self-contained cell (Shift+Enter in VS Code) or run the
whole file top to bottom. Assumes it is run from the repo root.
"""

# %% Setup — load the bead data and define the replicate columns
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BEADS_XLSX = os.path.join("MultiBatchResults", "beads_clean.xlsx")
RESULTS_DIR = "MultiBatchResults"
FLOOR = 3000  # "not detected" placeholder / lower measurement floor

REP5 = ["5Beads_V1P0113_1", "5Beads_V1P0113_2", "5Beads_V1P0113_3"]
REP10 = ["10Beads_V1P0113_1", "10Beads_V1P0113_2", "10Beads_V1P0113_3"]
ALL6 = REP5 + REP10

beads = pd.read_excel(BEADS_XLSX)
print(f"{len(beads):,} rows loaded from {BEADS_XLSX}")
for c in ALL6:
    at_floor = (beads[c] == FLOOR).mean()
    print(f"  {c:20} {beads[c].isna().sum():>3} NaN | {at_floor:5.1%} at floor {FLOOR}")


# %% 1. Histograms — each of the 6 replicate columns + pooled total
# Log-x because intensities span ~5 orders of magnitude with a spike at the
# 3000 floor. The floor pile-up is the leftmost bar and is expected.
_pos = beads[ALL6].where(beads[ALL6] > 0)
lo = np.nanmin(_pos.values)
hi = np.nanmax(_pos.values)
bins = np.logspace(np.log10(lo), np.log10(hi), 60)

fig, axes = plt.subplots(2, 3, figsize=(15, 7), sharex=True, sharey=True)
for ax, c in zip(axes.ravel(), ALL6):
    ax.hist(beads[c].clip(lower=lo), bins=bins)
    ax.axvline(FLOOR, color="r", ls="--", lw=1, label=f"floor {FLOOR}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(f"{c}\n{(beads[c] == FLOOR).mean():.0%} at floor")
    ax.set_xlabel("intensity")
    ax.legend(loc="upper right", fontsize=8)
fig.suptitle("Bead-only intensity distribution per replicate column")
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "bead_hist_per_column.png"), dpi=150)
plt.show()

# Pooled "total" — every value from all 6 columns in one histogram
fig, ax = plt.subplots(figsize=(8, 4))
ax.hist(beads[ALL6].values.ravel(), bins=bins)
ax.axvline(FLOOR, color="r", ls="--", lw=1, label=f"floor {FLOOR}")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("intensity (all 6 columns pooled)")
ax.set_ylabel("count")
ax.set_title(f"Bead-only intensity — pooled across all 6 replicate columns "
             f"({6 * len(beads):,} values)")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "bead_hist_total.png"), dpi=150)
plt.show()


# %% Helper — within-condition replicate agreement (abstract summaries only)
# For a set of 3 replicate columns: drop rows where all 3 are at the floor
# (no signal), then summarise agreement with (a) a CV histogram and (b) a
# printed correlation matrix. No rep-vs-rep scatter grid by request.
def within_condition_agreement(df, cols, name):
    d = df[cols]
    keep = ~(d == FLOOR).all(axis=1)          # at least one rep above floor
    d = d[keep]
    n = len(d)
    print(f"\n[{name}] {n:,} rows with signal "
          f"({keep.mean():.1%} of {len(df):,}; rest all at floor)")

    # CV across the 3 reps: std/mean per row. Low CV -> reps agree.
    cv = d.std(axis=1, ddof=1) / d.mean(axis=1)
    print(f"[{name}] median CV {cv.median():.1%} | "
          f"CV<=20%: {(cv <= 0.20).mean():.1%} | CV>50%: {(cv > 0.50).mean():.1%}")

    # Correlation matrix on log10 intensities (abstract agreement summary)
    logd = np.log10(d.clip(lower=FLOOR))
    print(f"[{name}] pairwise Pearson(log10):")
    print(logd.corr().round(3).to_string())

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(cv.clip(upper=2), bins=60)
    ax.axvline(cv.median(), color="r", ls="--", label=f"median = {cv.median():.1%}")
    ax.set_xlabel("CV across the 3 replicates (std/mean)")
    ax.set_ylabel("number of compounds")
    ax.set_title(f"{name}: replicate consistency (lower CV = more consistent)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, f"bead_cv_{name}.png"), dpi=150)
    plt.show()
    return cv


# %% 2. Do the three 5-bead replicates agree?
cv5 = within_condition_agreement(beads, REP5, "5bead")


# %% 3. Do the three 10-bead replicates agree?
cv10 = within_condition_agreement(beads, REP10, "10bead")


# %% 4. Do the 5-bead and 10-bead conditions agree?
# Collapse each condition to its 3-rep mean, keep rows where BOTH means are
# above the floor, then compare. One summary scatter + a ratio histogram.
from scipy.stats import pearsonr, spearmanr

mean5 = beads[REP5].mean(axis=1)
mean10 = beads[REP10].mean(axis=1)
keep = (mean5 > FLOOR) & (mean10 > FLOOR)
m5, m10 = mean5[keep], mean10[keep]
print(f"\n[5 vs 10] {keep.sum():,} compounds above floor in BOTH conditions "
      f"({keep.mean():.1%} of {len(beads):,})")

pear = pearsonr(np.log10(m5), np.log10(m10))[0]
spear = spearmanr(m5, m10)[0]
ratio = np.log2(m10 / m5)
print(f"[5 vs 10] Pearson(log10) {pear:.3f} | Spearman {spear:.3f}")
print(f"[5 vs 10] median log2(10-bead / 5-bead) = {ratio.median():+.2f} "
      f"(= {2 ** ratio.median():.2f}x); within 2x: {(ratio.abs() <= 1).mean():.1%}")

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ax = axes[0]
ax.scatter(m5, m10, s=3, alpha=min(0.5, 5000 / max(len(m5), 1)), edgecolors="none")
lo = min(m5.min(), m10.min())
hi = max(m5.max(), m10.max())
ax.plot([lo, hi], [lo, hi], "r--", lw=1, label="y = x")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("mean 5-bead intensity")
ax.set_ylabel("mean 10-bead intensity")
ax.set_title(f"5-bead vs 10-bead (mean of 3 reps)\n"
             f"Pearson(log)={pear:.3f}  Spearman={spear:.3f}")
ax.legend(loc="upper left")

ax = axes[1]
ax.hist(ratio.clip(-6, 6), bins=60)
ax.axvline(0, color="k", lw=1)
ax.axvline(ratio.median(), color="r", ls="--",
           label=f"median = {ratio.median():+.2f}")
ax.set_xlabel("log2(10-bead / 5-bead)")
ax.set_ylabel("number of compounds")
ax.set_title("Condition ratio (0 = identical)")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "bead_5_vs_10.png"), dpi=150)
plt.show()


# %% 5. Bead binders by threshold — mean over all 6 columns
# For each compound take the mean of all 6 replicate values (both conditions),
# then count how many exceed a signal threshold ("bead binders"). Because the
# right threshold is a judgement call, we sweep it and draw an exceedance curve
# (binder count vs threshold) so any cutoff can be read off one plot. A single
# THRESHOLD (default 10,000) is highlighted and reported.
THRESHOLD = 10_000

mean6 = beads[ALL6].mean(axis=1)
n_total = len(mean6)
n_binders = int((mean6 > THRESHOLD).sum())
print(f"\nmean-of-6 over {n_total:,} compounds | "
      f"binders (mean > {THRESHOLD:,}): {n_binders:,} ({n_binders / n_total:.1%})")

# Table at round thresholds so the numbers are readable without the plot
grid = [3_000, 5_000, 10_000, 20_000, 50_000, 100_000, 500_000, 1_000_000]
print(f"{'threshold':>12} {'binders':>10} {'fraction':>10}")
for t in grid:
    c = int((mean6 > t).sum())
    print(f"{t:>12,} {c:>10,} {c / n_total:>9.1%}")

# Exceedance curve: for every threshold on a log grid, how many compounds have
# mean-of-6 above it. Monotonically decreasing — read the binder count for any
# cutoff straight off the line.
ts = np.logspace(np.log10(FLOOR), np.log10(mean6.max()), 200)
counts = [(mean6.values[:, None] > ts).sum(axis=0)][0]  # vectorised over ts

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(ts, counts, lw=2)
ax.axvline(THRESHOLD, color="r", ls="--",
           label=f"threshold {THRESHOLD:,} -> {n_binders:,} binders")
ax.scatter([THRESHOLD], [n_binders], color="r", zorder=5)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("threshold on mean-of-6 intensity")
ax.set_ylabel("number of bead binders (mean-of-6 > threshold)")
ax.set_title("Bead binders vs threshold (exceedance curve)")
ax.grid(True, which="both", alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "bead_binders_by_threshold.png"), dpi=150)
plt.show()


# %%
