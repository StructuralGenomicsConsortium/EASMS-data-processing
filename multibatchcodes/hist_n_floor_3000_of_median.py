"""
Histogram of n_floor_3000_of_median (per-molecule count of proteins where the ASMS
median(REP1-3) sits at the 3000 floor) for the 1647 binder molecules.
Reads the exported matrix CSV; saves a PNG next to the other plots.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = r"d:/0000-UHN/EASMS-data-processing/EASMS-data-processing/MultiBatch_20to36_plots"
CSV = os.path.join(OUT, "MultiBatch_20to36_matrix_sorted_floor_full.csv")
PNG = os.path.join(OUT, "hist_n_floor_3000_of_median.png")

v = pd.read_csv(CSV, usecols=["n_floor_3000_of_median"])["n_floor_3000_of_median"].to_numpy()
print(f"n={len(v)} | min={v.min()} max={v.max()} mean={v.mean():.1f} median={np.median(v):.0f}")

fig, ax = plt.subplots(figsize=(10, 5))
bins = np.arange(v.min(), v.max() + 2) - 0.5          # one bin per integer count
ax.hist(v, bins=bins, color="#4C78A8", edgecolor="white", linewidth=0.3)
ax.axvline(v.mean(),   color="crimson", ls="--", lw=1, label=f"mean = {v.mean():.1f}")
ax.axvline(np.median(v), color="black", ls=":",  lw=1, label=f"median = {np.median(v):.0f}")
ax.set_xlabel("n_floor_3000_of_median  (# of 136 proteins at the 3000 floor)")
ax.set_ylabel("number of molecules")
ax.set_title(f"Distribution of n_floor_3000_of_median  (n={len(v)} binder molecules)")
ax.legend()
fig.tight_layout()
fig.savefig(PNG, dpi=150)
print("saved", PNG)
