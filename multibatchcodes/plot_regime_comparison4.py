# -*- coding: utf-8 -*-
"""Compare binder-calling regimes incl. percentile normalization (modified4).

Reads MultiBatch_20to22_modified4/ and draws one grouped bar chart of hit counts
per batch (+ ALL) for four regimes, log y:

    within-batch                (label_within)          -- batch-controlled baseline
    across-batch raw            (label_across)          -- inflated by the hot batch
    across-batch median-norm    (label_across_norm)     -- level corrected (modified3)
    across-batch percentile-norm(label_across_pctnorm)  -- level + detection-rate (modified4)

Run from the repo root:  python multibatchcodes/plot_regime_comparison4.py
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DIR = "MultiBatch_20to22_modified4"
RESULTS_DIR = "MultiBatchResults"
os.makedirs(RESULTS_DIR, exist_ok=True)

cols = ["ASMS_BATCH_NAME", "label_within", "label_across",
        "label_across_norm", "label_across_pctnorm"]
df = pd.concat([pd.read_csv(p, usecols=cols) for p in sorted(glob.glob(os.path.join(DIR, "*.csv")))],
               ignore_index=True)

batches = sorted(df["ASMS_BATCH_NAME"].unique())
groups = batches + ["ALL"]
regimes = [("label_within", "within-batch", "#2a7fb8"),
           ("label_across", "across raw", "#c1121f"),
           ("label_across_norm", "across median-norm", "#e08a1e"),
           ("label_across_pctnorm", "across percentile-norm", "#2ca25f")]


def count(col, g):
    return int((df[col] if g == "ALL" else df.loc[df["ASMS_BATCH_NAME"] == g, col]).sum())


x = np.arange(len(groups))
w = 0.2
fig, ax = plt.subplots(figsize=(11, 5.5))
for k, (col, label, color) in enumerate(regimes):
    vals = [count(col, g) for g in groups]
    bars = ax.bar(x + (k - 1.5) * w, vals, w, label=label, color=color)
    ax.bar_label(bars, fmt="%d", fontsize=7, padding=2)

ax.set_yscale("log")
ax.set_xticks(x)
ax.set_xticklabels(groups)
ax.set_ylabel("molecules called as binders (log scale)")
ax.set_xlabel("batch")
ax.set_title("Binder calls by regime — within vs across (raw / median-norm / percentile-norm)\n"
             "MultiBatch_20to22")
ax.legend(ncol=2, fontsize=9)
ax.axvline(len(batches) - 0.5, color="0.85", lw=1)
fig.tight_layout()
out = os.path.join(RESULTS_DIR, "regime_comparison_modified4.png")
fig.savefig(out, dpi=150)
plt.show()
print("saved", out)
print(pd.DataFrame({label: [count(col, g) for g in groups] for col, label, _ in regimes},
                   index=groups).to_string())
