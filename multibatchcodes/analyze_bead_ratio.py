# -*- coding: utf-8 -*-
"""Bead-binder specificity ratio for the within-batch hits.

Instead of removing a compound just because it shows ANY bead-only signal, we
keep it unless its protein signal fails to clear the bead background by a margin:

    bead_signal = max( median(5-bead reps), median(10-bead reps) )   # conservative
    bead_ratio  = target_median (on protein) / bead_signal
    bead_binder_flag = (bead_signal > 3000)  AND  (bead_ratio < threshold)

Compounds not detected on beads (bead at 3000 floor) are never flagged (they have
no bead binding); their ratio = target/3000 is a harmless lower bound.

Matched by COMPOUND_ID (= "SGC ID for Component" in the bead file). Reports how
many within-batch hits are bead binders at 5x and 10x, and plots the target-vs-
bead scatter + bead_ratio distribution so the threshold can be picked from shape.

Run from the repo root:  python multibatchcodes/analyze_bead_ratio.py
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS_DIR = "MultiBatchResults"
FLOOR = 3000
THRESHOLDS = [5, 10]
C5 = ["5Beads_V1P0113_1", "5Beads_V1P0113_2", "5Beads_V1P0113_3"]
C10 = ["10Beads_V1P0113_1", "10Beads_V1P0113_2", "10Beads_V1P0113_3"]

# --- bead signal per compound (conservative: max of the two condition medians) ---
b = pd.read_excel(os.path.join(RESULTS_DIR, "beads_clean.xlsx"))
b["bead_signal"] = np.maximum(b[C5].median(axis=1), b[C10].median(axis=1))
b["bead_detected"] = b["bead_signal"] > FLOOR
bead = (b[["SGC ID for Component", "bead_signal", "bead_detected"]]
        .rename(columns={"SGC ID for Component": "COMPOUND_ID"})
        .dropna(subset=["COMPOUND_ID"]).drop_duplicates("COMPOUND_ID"))
print(f"bead compounds: {len(bead):,} | detected on beads: "
      f"{int(bead['bead_detected'].sum()):,} ({bead['bead_detected'].mean():.1%})")

# --- protein data (within-batch hits) ---
cols = ["COMPOUND_ID", "ASMS_BATCH_NAME", "target_median", "label_within"]
df = pd.concat([pd.read_csv(p, usecols=cols)
                for p in glob.glob("MultiBatch_20to22_modified2/*.csv")], ignore_index=True)
m = df.merge(bead, on="COMPOUND_ID", how="left")
m["bead_ratio"] = m["target_median"] / m["bead_signal"]

hits = m[m["label_within"] == 1].copy()
det = hits[hits["bead_detected"] == True]
print(f"\nwithin-batch hits: {len(hits):,} | matched to bead: {int(hits['bead_signal'].notna().sum()):,} "
      f"| detected on beads: {len(det):,}")
for t in THRESHOLDS:
    f = (hits["bead_detected"] == True) & (hits["bead_ratio"] < t)
    print(f"  bead_binder flag (ratio < {t}): {int(f.sum()):,}  ({f.mean():.1%} of hits)")

# --- plots ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# (1) target vs bead scatter (hits only)
ax = axes[0]
kept = hits[~((hits["bead_detected"] == True) & (hits["bead_ratio"] < 10))]
flag = hits[(hits["bead_detected"] == True) & (hits["bead_ratio"] < 10)]
ax.scatter(kept["bead_signal"].clip(lower=FLOOR), kept["target_median"], s=10,
           alpha=0.5, color="#2ca25f", label=f"keep (ratio>=10 or bead floor)  n={len(kept)}")
ax.scatter(flag["bead_signal"].clip(lower=FLOOR), flag["target_median"], s=12,
           alpha=0.6, color="#c1121f", label=f"bead binder (ratio<10)  n={len(flag)}")
lo, hi = FLOOR, hits["target_median"].max()
xs = np.array([lo, hi])
for t, ls in [(1, ":"), (5, "--"), (10, "-")]:
    ax.plot(xs, t * xs, "k", ls=ls, lw=1, label=f"ratio = {t}")
ax.axvline(FLOOR, color="0.7", lw=1)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("bead_signal (max of 5/10-bead medians; 3000 = not detected)")
ax.set_ylabel("target_median on protein")
ax.set_title("Within-batch hits: protein vs bead signal\n(points below a line fail that ratio)")
ax.legend(fontsize=8, loc="lower right")

# (2) bead_ratio distribution for hits detected on beads
ax = axes[1]
r = det["bead_ratio"]
bins = np.logspace(np.log10(max(r.min(), 0.1)), np.log10(r.max()), 40)
ax.hist(r, bins=bins, color="#2a7fb8")
for t, c in [(5, "orange"), (10, "red")]:
    ax.axvline(t, color=c, ls="--", label=f"ratio = {t} ({int((r<t).sum())} below)")
ax.set_xscale("log")
ax.set_xlabel("bead_ratio = target_median / bead_signal  (hits detected on beads)")
ax.set_ylabel("number of hits")
ax.set_title(f"Specificity of the {len(det)} bead-detected hits\n(left of a line = likely bead binder)")
ax.legend()

fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "bead_ratio_analysis.png"), dpi=150)
print("\nsaved MultiBatchResults/bead_ratio_analysis.png")
