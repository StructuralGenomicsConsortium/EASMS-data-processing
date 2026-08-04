# -*- coding: utf-8 -*-
"""Protein replicate-intensity range — MultiBatch_20to22.

What is the dynamic range of the mass-spec intensities in the protein data?
For every protein file in MultiBatch_20to22/ we pool the three replicate columns
(POS_INT_REP1/2/3) into one set of intensity values, DROP the 3000 floor
("not detected"), and then:

    1. print + save a per-protein min / max table (plus an "ALL" row),
    2. draw a per-protein box plot of intensity (log x, whiskers = full min..max),
    3. draw one pooled histogram over all proteins.

Mirrors the bead-only intensity look, but per protein and for the raw replicate
intensity (not the CV). Run from the repo root:

    python multibatchcodes/protein_intensity_range.py
"""

# %% Setup — pool replicate intensities per protein, drop the 3000 floor
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

MULTIBATCH_DIR = "MultiBatch_20to22"     # input: per-protein CSVs
RESULTS_DIR = "MultiBatchResults"        # output: plots + table
FLOOR = 3000                             # "not detected" placeholder, removed
REP = ["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]
os.makedirs(RESULTS_DIR, exist_ok=True)

# intensity[protein]     = all replicate values above the floor
# intensity_pos[protein] = the same, but only from BINARY_LABEL == 1 rows (binders)
intensity = {}
intensity_pos = {}
for path in sorted(glob.glob(os.path.join(MULTIBATCH_DIR, "*.csv"))):
    df = pd.read_csv(path, usecols=REP + ["BINARY_LABEL"])
    name = os.path.splitext(os.path.basename(path))[0]

    v = df[REP].to_numpy().ravel()
    intensity[name] = v[np.isfinite(v) & (v != FLOOR)]

    is_pos = pd.to_numeric(df["BINARY_LABEL"], errors="coerce") == 1
    vp = df.loc[is_pos, REP].to_numpy().ravel()
    intensity_pos[name] = vp[np.isfinite(vp) & (vp != FLOOR)]

    print(f"{name:12} {len(intensity[name]):>8,} above floor | "
          f"{len(intensity_pos[name]):>7,} from BINARY_LABEL=1")

all_vals = np.concatenate(list(intensity.values()))
all_pos = np.concatenate(list(intensity_pos.values()))
print(f"\n{len(intensity)} proteins | {len(all_vals):,} total above floor {FLOOR} | "
      f"{len(all_pos):,} from BINARY_LABEL=1")


# %% 1. Min / max table per protein (+ ALL), printed and saved to CSV
rows = []
for name, v in intensity.items():
    rows.append({
        "protein": name,
        "n": len(v),
        "min": v.min(),
        "median": np.median(v),
        "max": v.max(),
    })
rows.append({
    "protein": "ALL", "n": len(all_vals),
    "min": all_vals.min(), "median": np.median(all_vals), "max": all_vals.max(),
})
table = pd.DataFrame(rows)
with pd.option_context("display.float_format", lambda x: f"{x:,.1f}"):
    print(table.to_string(index=False))
_out = os.path.join(RESULTS_DIR, "protein_intensity_minmax.csv")
table.to_csv(_out, index=False)
print(f"\nSaved {_out}")


# %% 2. Per-protein box plot of replicate intensity (log x, whiskers = min..max)
# Proteins ordered by median so the spread is easy to scan. whis=(0,100) makes
# the whiskers reach the true min and max, i.e. the full dynamic range per
# protein (no points hidden as outliers).
order = sorted(intensity, key=lambda k: np.median(intensity[k]))
data = [intensity[k] for k in order]
overall_median = np.median(all_vals)

fig, ax = plt.subplots(figsize=(9, 0.4 * len(order) + 1))
ax.boxplot(data, vert=False, showfliers=False, whis=(0, 100))
ax.set_yticks(range(1, len(order) + 1))
ax.set_yticklabels(order)
ax.set_xscale("log")
ax.axvline(overall_median, color="r", ls="--",
           label=f"overall median = {overall_median:,.0f}")
ax.set_xlabel("replicate intensity (POS_INT_REP1/2/3, floor 3000 removed)")
ax.set_ylabel("protein (file)")
ax.set_title("Per-protein replicate intensity range — MultiBatch_20to22")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "protein_intensity_box.png"), dpi=150)
plt.show()


# %% 3. Pooled histogram of all replicate intensities (log x)
bins = np.logspace(np.log10(all_vals.min()), np.log10(all_vals.max()), 60)
fig, ax = plt.subplots(figsize=(8, 4))
ax.hist(all_vals, bins=bins)
ax.set_xscale("log")
ax.set_yscale("log")
ax.axvline(np.median(all_vals), color="r", ls="--",
           label=f"median = {np.median(all_vals):,.0f}")
ax.set_xlabel("replicate intensity (all proteins pooled, floor 3000 removed)")
ax.set_ylabel("count")
ax.set_title(f"Protein replicate intensity distribution — MultiBatch_20to22 "
             f"({len(all_vals):,} values)")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "protein_intensity_hist_all.png"), dpi=150)
plt.show()


# %% 4. Histogram of intensity for binders only (BINARY_LABEL == 1)
# Same log-x bins as cell 3 so it lines up with the full distribution, which is
# drawn faintly behind for context. Binders sit at the high-intensity end.
fig, ax = plt.subplots(figsize=(8, 4))
ax.hist(all_vals, bins=bins, color="0.8", label=f"all ({len(all_vals):,})")
ax.hist(all_pos, bins=bins, color="C1", alpha=0.8,
        label=f"BINARY_LABEL=1 ({len(all_pos):,})")
ax.set_xscale("log")
ax.set_yscale("log")
ax.axvline(np.median(all_pos), color="r", ls="--",
           label=f"binder median = {np.median(all_pos):,.0f}")
ax.set_xlabel("replicate intensity (floor 3000 removed)")
ax.set_ylabel("count")
ax.set_title("Binder (BINARY_LABEL=1) replicate intensity — MultiBatch_20to22")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "protein_intensity_hist_positives.png"), dpi=150)
plt.show()


# %%
