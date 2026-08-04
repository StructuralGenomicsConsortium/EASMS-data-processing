# -*- coding: utf-8 -*-
"""Is retention time (RT) associated with binder calls / enrichment?

A bias check: if hits or high enrichment cluster at particular retention times,
that hints at a chromatographic / MS artifact rather than real binding.

RT (min) vs:
  - BINARY_LABEL   (original binder call)  -> point-biserial r, Mann-Whitney, hit-rate-by-RT
  - label_within   (our median call)       -> hit-rate-by-RT
  - enrichment_within / enrichment_across  -> Spearman, median-enrichment-by-RT

Run from the repo root:  python multibatchcodes/analyze_rt_correlation.py
"""

import glob
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pointbiserialr, spearmanr, mannwhitneyu

RESULTS_DIR = "MultiBatchResults"
cols = ["RT (min)", "BINARY_LABEL", "enrichment_within", "enrichment_across",
        "label_within", "label_across"]
df = pd.concat([pd.read_csv(p, usecols=cols)
                for p in glob.glob("MultiBatch_20to22_modified5/*.csv")], ignore_index=True)
df = df.dropna(subset=["RT (min)"])
rt = df["RT (min)"].to_numpy()
bl = pd.to_numeric(df["BINARY_LABEL"], errors="coerce")
lw = df["label_within"]
print(f"{len(df):,} rows | RT range {rt.min():.2f}-{rt.max():.2f} min | "
      f"BINARY_LABEL=1: {int((bl==1).sum()):,} | label_within=1: {int((lw==1).sum()):,}")

# --- correlations ---
r_pb, p_pb = pointbiserialr(bl.fillna(0), rt)
u, p_mwu = mannwhitneyu(rt[bl == 1], rt[bl == 0], alternative="two-sided")
print(f"\nRT vs BINARY_LABEL:")
print(f"  point-biserial r = {r_pb:+.3f} (p={p_pb:.1e})")
print(f"  median RT  binders {np.median(rt[bl==1]):.2f}  vs  non {np.median(rt[bl==0]):.2f} min "
      f"(Mann-Whitney p={p_mwu:.1e})")

for col in ["enrichment_within", "enrichment_across"]:
    e = df[col]
    m = e.notna() & (e > 0)
    rho, pr = spearmanr(rt[m], np.log10(e[m]))
    print(f"RT vs {col}: Spearman rho = {rho:+.3f} (p={pr:.1e}, n={int(m.sum()):,})")

# --- plots ---
nb = 24
edges = np.linspace(rt.min(), rt.max(), nb + 1)
mid = (edges[:-1] + edges[1:]) / 2
which = np.digitize(rt, edges) - 1
which = np.clip(which, 0, nb - 1)

fig, axes = plt.subplots(1, 3, figsize=(17, 5))

# (1) hit rate vs RT
ax = axes[0]
for lab, series, color in [("BINARY_LABEL", (bl == 1).to_numpy(), "#c1121f"),
                           ("label_within", (lw == 1).to_numpy(), "#2a7fb8")]:
    rate = np.array([series[which == i].mean() if (which == i).any() else np.nan for i in range(nb)])
    ax.plot(mid, rate, "-o", ms=4, color=color, label=lab)
ax.set_xlabel("RT (min)")
ax.set_ylabel("hit rate (fraction called)")
ax.set_title("Hit rate vs retention time")
ax.legend()

# (2) RT distribution by BINARY_LABEL
ax = axes[1]
bins = np.linspace(rt.min(), rt.max(), 40)
ax.hist(rt[bl == 0], bins=bins, density=True, histtype="step", lw=2, color="0.5", label="non-binder")
ax.hist(rt[bl == 1], bins=bins, density=True, histtype="step", lw=2, color="#c1121f", label="binder")
ax.set_xlabel("RT (min)")
ax.set_ylabel("density")
ax.set_title(f"RT distribution by BINARY_LABEL\n(point-biserial r={r_pb:+.3f})")
ax.legend()

# (3) median enrichment vs RT
ax = axes[2]
for col, color in [("enrichment_within", "#2a7fb8"), ("enrichment_across", "#e08a1e")]:
    e = df[col].to_numpy()
    med = np.array([np.nanmedian(e[which == i]) if (which == i).any() else np.nan for i in range(nb)])
    ax.plot(mid, med, "-o", ms=4, color=color, label=col)
ax.axhline(1, color="0.7", lw=1)
ax.set_xlabel("RT (min)")
ax.set_ylabel("median enrichment")
ax.set_title("Median enrichment vs retention time")
ax.legend()

fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "rt_correlation.png"), dpi=150)
print("\nsaved MultiBatchResults/rt_correlation.png")
