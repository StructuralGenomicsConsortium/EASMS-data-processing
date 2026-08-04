# -*- coding: utf-8 -*-
"""Compare binder-calling regimes: within-batch vs across-batch vs across-normalized.

Reads MultiBatch_20to22_modified3/ (which carries label_within, label_across, and
label_across_norm) and draws one grouped bar chart of hit counts per batch (+ ALL)
for the three regimes, log y so all three are readable together.

Run from the repo root:  python multibatchcodes/plot_regime_comparison.py
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DIR = "MultiBatch_20to22_modified3"
RESULTS_DIR = "MultiBatchResults"
os.makedirs(RESULTS_DIR, exist_ok=True)

cols = ["ASMS_BATCH_NAME", "label_within", "label_across", "label_across_norm"]
df = pd.concat([pd.read_csv(p, usecols=cols) for p in sorted(glob.glob(os.path.join(DIR, "*.csv")))],
               ignore_index=True)

batches = sorted(df["ASMS_BATCH_NAME"].unique())
groups = batches + ["ALL"]
regimes = [("label_within", "within-batch", "#2a7fb8"),
           ("label_across", "across-batch (raw)", "#c1121f"),
           ("label_across_norm", "across-batch (normalized)", "#e08a1e")]

# counts[regime][group]
def count(col, g):
    return int((df[col] if g == "ALL" else df.loc[df["ASMS_BATCH_NAME"] == g, col]).sum())

x = np.arange(len(groups))
w = 0.26
fig, ax = plt.subplots(figsize=(10, 5.5))
for k, (col, label, color) in enumerate(regimes):
    vals = [count(col, g) for g in groups]
    bars = ax.bar(x + (k - 1) * w, vals, w, label=label, color=color)
    ax.bar_label(bars, fmt="%d", fontsize=8, padding=2)

ax.set_yscale("log")
ax.set_xticks(x)
ax.set_xticklabels(groups)
ax.set_ylabel("molecules called as binders (log scale)")
ax.set_xlabel("batch")
ax.set_title("Binder calls by regime: within-batch vs across-batch vs across-normalized\n"
             "MultiBatch_20to22  (within is balanced; raw-across inflates sgcto_22; "
             "normalization collapses it)")
ax.legend()
ax.axvline(len(batches) - 0.5, color="0.8", lw=1)  # separate ALL
fig.tight_layout()
out = os.path.join(RESULTS_DIR, "regime_comparison.png")
fig.savefig(out, dpi=150)
plt.show()
print("saved", out)
print(pd.DataFrame({label: [count(col, g) for g in groups] for col, label, _ in regimes},
                   index=groups).to_string())
