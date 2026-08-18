# -*- coding: utf-8 -*-
"""Plot: the per-protein intensity median is a BATCH effect — MultiBatch_20to22.

Redraws MultiBatchResults/protein_intensity_box.png (from protein_intensity_range.py:
pool POS_INT_REP1/2/3, drop the 3000 floor, order proteins by median) with every
protein COLOURED BY ITS ASMS_BATCH_NAME, which makes the answer visible: the top 8
proteins of that plot are exactly the 8 proteins of sgcto_22.

Three panels:
  A  per-protein box plot (log x, whiskers = full min..max), coloured by batch
  B  batch hotness vs censoring: median of detected values against the share of
     replicate values stuck at the 3000 floor, one point per protein
  C  per-batch floor fraction, with the detected median direct-labelled

B and C carry the second point: the hot batch is not just scaled up, it is also
far LESS censored (29% floor vs 82%), because a late read pushes borderline
values up over the 3000 "not detected" placeholder.

Colour: dataviz categorical slots 1-3 assigned by sorted batch name (by entity,
never by rank), validated all-pairs on the light surface.

Run from the repo root:  python multibatchcodes/plot_protein_median_by_batch.py
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

MULTIBATCH_DIR = "MultiBatch_20to22"
RESULTS_DIR = "MultiBatchResults"
FLOOR = 3000
REP = ["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]
os.makedirs(RESULTS_DIR, exist_ok=True)

# --- dataviz tokens (light surface) ------------------------------------------
SURFACE   = "#fcfcfb"
INK       = "#0b0b0b"
INK_2     = "#52514e"
MUTED     = "#898781"
GRID      = "#e1e0d9"
BASELINE  = "#c3c2b7"
SLOTS     = ["#2a78d6", "#eb6834", "#1baf7a"]      # blue, orange, aqua

plt.rcParams.update({
    "font.family": ["Segoe UI", "DejaVu Sans", "sans-serif"],
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "axes.edgecolor": BASELINE, "axes.linewidth": 0.8,
    "axes.labelcolor": INK_2, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK_2, "ytick.labelcolor": INK_2,
    "grid.color": GRID, "grid.linewidth": 0.8, "grid.linestyle": "-",
})

# %% Load — per protein: all replicate values, detected values, batch
recs = {}
for path in sorted(glob.glob(os.path.join(MULTIBATCH_DIR, "*.csv"))):
    name = os.path.splitext(os.path.basename(path))[0]
    df = pd.read_csv(path, usecols=REP + ["ASMS_BATCH_NAME"])
    v_all = df[REP].to_numpy(dtype=float).ravel()
    v_all = v_all[np.isfinite(v_all)]
    v_det = v_all[v_all != FLOOR]
    batch = str(df["ASMS_BATCH_NAME"].dropna().iloc[0])
    if df["ASMS_BATCH_NAME"].nunique() != 1:
        raise SystemExit(f"{name}: expected 1 batch, got {df['ASMS_BATCH_NAME'].unique()}")
    recs[name] = {"batch": batch, "det": v_det,
                  "median": float(np.median(v_det)),
                  "floor_frac": 1 - len(v_det) / len(v_all)}

batches = sorted({r["batch"] for r in recs.values()})
if len(batches) > len(SLOTS):
    raise SystemExit(f"{len(batches)} batches > {len(SLOTS)} validated colour slots")
COLOR = dict(zip(batches, SLOTS))                   # by entity (sorted name), not rank

tab = (pd.DataFrame([{"protein": k, "batch": v["batch"], "median_detected": v["median"],
                      "floor_frac": v["floor_frac"], "n_detected": len(v["det"])}
                     for k, v in recs.items()])
       .sort_values("median_detected", ascending=False).reset_index(drop=True))
tab.insert(0, "rank", tab.index + 1)

bstat = (tab.groupby("batch")
         .agg(n_proteins=("protein", "size"), mean_rank=("rank", "mean"))
         .reset_index())
for b in batches:                                   # pooled, not a mean of per-protein stats
    v_all = np.concatenate([np.append(recs[k]["det"], []) for k in recs if recs[k]["batch"] == b])
    n_all = sum(len(recs[k]["det"]) / (1 - recs[k]["floor_frac"]) for k in recs if recs[k]["batch"] == b)
    bstat.loc[bstat["batch"] == b, "median_detected"] = float(np.median(v_all))
    bstat.loc[bstat["batch"] == b, "floor_frac"] = 1 - len(v_all) / n_all
bstat = bstat.sort_values("median_detected", ascending=False).reset_index(drop=True)

_csv = os.path.join(RESULTS_DIR, "protein_median_by_batch_table.csv")
tab.to_csv(_csv, index=False)                        # table view (relief rule)

# %% Figure
fig = plt.figure(figsize=(13.5, 8.8))
gs = fig.add_gridspec(2, 2, width_ratios=[1.45, 1], height_ratios=[1.9, 1],
                      wspace=0.34, hspace=0.42,
                      left=0.085, right=0.975, top=0.855, bottom=0.075)
axA = fig.add_subplot(gs[:, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 1])

# --- A: per-protein box plot, ordered by median, coloured by batch -----------
order = list(tab["protein"])[::-1]                   # ascending -> highest at top
for i, name in enumerate(order, start=1):
    r = recs[name]
    c = COLOR[r["batch"]]
    bp = axA.boxplot([r["det"]], positions=[i], vert=False, widths=0.60,
                     showfliers=False, whis=(0, 100), patch_artist=True)
    bp["boxes"][0].set(facecolor=c, alpha=0.20, edgecolor=c, linewidth=1.3)
    for w in bp["whiskers"] + bp["caps"]:
        w.set(color=c, linewidth=1.0, alpha=0.85)
    bp["medians"][0].set(color=c, linewidth=2.2, solid_capstyle="butt")

overall = float(np.median(np.concatenate([r["det"] for r in recs.values()])))
axA.axvline(overall, color=MUTED, lw=1.0, zorder=0)
axA.text(overall * 1.06, 0.35, f"overall median {overall:,.0f}",
         color=INK_2, fontsize=8.5, va="bottom")

axA.set_yticks(range(1, len(order) + 1))
axA.set_yticklabels(order, fontsize=9)
for lbl in axA.get_yticklabels():                    # tick label inherits its batch colour? no —
    lbl.set_color(INK_2)                             # text wears text tokens, mark carries identity
axA.set_xscale("log")
_lo = min(r["det"].min() for r in recs.values())
_hi = max(r["det"].max() for r in recs.values())
axA.set_xlim(_lo * 0.88, _hi * 1.25)                 # never clip a min..max whisker
axA.set_ylim(0.3, len(order) + 0.7)
axA.xaxis.grid(True, which="major", zorder=0)
axA.set_axisbelow(True)
for s in ("top", "right", "left"):
    axA.spines[s].set_visible(False)
axA.tick_params(axis="y", length=0)
axA.set_xlabel("replicate intensity  (POS_INT_REP1/2/3, 3000 floor removed)", fontsize=9.5)
axA.set_title("Per-protein intensity range, ordered by median",
              fontsize=11, color=INK, loc="left", pad=8)

# direct label: the top 8 are exactly one batch (relief + the headline).
# Drawn just OUTSIDE the axes (x in axes fraction, y in data coords) so it never
# collides with the min..max whiskers, which run the full width of the panel.
hot = bstat.loc[0, "batch"]
n_hot = int(bstat.loc[0, "n_proteins"])
top = len(order)
tr = axA.get_yaxis_transform()
axA.plot([1.006, 1.020, 1.020, 1.006],
         [top - n_hot + 0.62, top - n_hot + 0.62, top + 0.38, top + 0.38],
         color=COLOR[hot], lw=1.3, transform=tr, clip_on=False)
axA.text(1.030, top - n_hot / 2 + 0.5, f"all {n_hot} proteins of {hot}",
         color=INK, fontsize=9.5, va="center", ha="center",
         rotation=270, transform=tr, clip_on=False)

# --- B: hotness vs censoring, one point per protein --------------------------
for b in batches:
    sub = tab[tab["batch"] == b]
    axB.scatter(sub["floor_frac"] * 100, sub["median_detected"],
                s=64, color=COLOR[b], alpha=0.95,
                edgecolor=SURFACE, linewidth=1.6, zorder=3)
# one direct label per cluster, nudged off the cluster it names
for b, dx, dy, ha in zip(bstat["batch"], (14, -14, 14), (11, 13, 13),
                         ("left", "right", "left")):
    sub = tab[tab["batch"] == b]
    peak = sub.loc[sub["median_detected"].idxmax()]        # anchor at the top-most dot
    axB.annotate(b, (peak["floor_frac"] * 100, peak["median_detected"]),
                 textcoords="offset points", xytext=(dx, dy),
                 ha=ha, fontsize=9, color=INK)
axB.set_yscale("log")
axB.set_xlim(10, 95)
axB.set_ylim(6.5e3, 4.4e4)
axB.grid(True, which="major", axis="both", zorder=0)
axB.set_axisbelow(True)
for s in ("top", "right"):
    axB.spines[s].set_visible(False)
axB.set_xlabel("share of replicate values at the 3000 floor  (%)", fontsize=9.5)
axB.set_ylabel("median detected intensity", fontsize=9.5)
axB.set_title("Hot batches are also the LEAST censored",
              fontsize=11, color=INK, loc="left", pad=8)

# --- C: per-batch floor fraction, median direct-labelled ---------------------
y = np.arange(len(bstat))[::-1]
axC.barh(y, bstat["floor_frac"] * 100, height=0.42,
         color=[COLOR[b] for b in bstat["batch"]], alpha=0.9, zorder=3)
for yy, (_, r) in zip(y, bstat.iterrows()):
    axC.text(r["floor_frac"] * 100 + 2.0, yy, f"{r['floor_frac']:.0%}",
             va="center", fontsize=9, color=INK_2)
axC.set_yticks(y)
axC.set_yticklabels([f"{r.batch}\nmedian {r.median_detected:,.0f}"
                     for r in bstat.itertuples()], fontsize=8.5)
axC.tick_params(axis="y", length=0)
axC.set_xlim(0, 100)
axC.xaxis.grid(True, which="major", zorder=0)
axC.set_axisbelow(True)
for s in ("top", "right", "left"):
    axC.spines[s].set_visible(False)
axC.set_xlabel("replicate values at the 3000 floor  (%)", fontsize=9.5)
axC.set_title("Censoring by batch", fontsize=11, color=INK, loc="left", pad=8)

fig.suptitle("The per-protein intensity median tracks the BATCH, not the protein",
             fontsize=14, color=INK, x=0.085, ha="left", y=0.975)
fig.text(0.085, 0.938,
         f"MultiBatch_20to22  ·  {len(recs)} proteins, {len(batches)} batches of 8  ·  "
         f"batch explains the top of the ranking entirely (ranks 1-{n_hot} = {hot})",
         fontsize=9.5, color=INK_2, ha="left")

# shared figure legend — one per figure, above everything it scopes
legend_h = [Line2D([0], [0], color=COLOR[b], lw=2.8, solid_capstyle="butt",
                   label=f"{b}    median {bstat.loc[bstat.batch == b, 'median_detected'].iloc[0]:,.0f}"
                         f"    ·    floor {bstat.loc[bstat.batch == b, 'floor_frac'].iloc[0]:.0%}")
            for b in bstat["batch"]]
fig.legend(handles=legend_h, loc="upper left", bbox_to_anchor=(0.083, 0.917),
           ncol=3, fontsize=9.5, frameon=False, labelcolor=INK_2,
           handlelength=1.7, columnspacing=3.4, handletextpad=0.7)

_png = os.path.join(RESULTS_DIR, "protein_median_by_batch.png")
fig.savefig(_png, dpi=200, facecolor=SURFACE)
print(bstat.to_string(index=False))
print(f"\nSaved {_png}\nSaved {_csv}")
