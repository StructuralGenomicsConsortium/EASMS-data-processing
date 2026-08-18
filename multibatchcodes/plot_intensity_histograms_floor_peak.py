# -*- coding: utf-8 -*-
"""Is there a low-band peak above the 3000 floor that could serve as the floor value?

Idea under test (user's): the intensity histogram should be a mixture -- a low peak
made of non-binders / background, plus a right tail of real binders. If that low peak
is identifiable, its location is a better substitute for a floored value than the
arbitrary 3000, AND -- crucially -- estimated per batch it should SCALE with the batch
factor, which is exactly the property that makes the within-batch ratio cancel.

The risk: values are CLAMPED onto 3000, not merely bounded by it, so the background
peak may sit underneath the spike and be invisible. These plots check that.

Two figures:
  ..._by_batch.png     one panel per batch, 8 proteins pooled
  ..._by_protein.png   small multiples, one panel per protein
Both log x (intensity) and log y (counts), matching protein_intensity_range.py, so
the point mass at the floor and the shape of the detected part are visible together.

The clamp is a POINT MASS, so it is drawn as its own bar just left of the detected
range rather than being folded into the first histogram bin.

Reported per group: floor share, modal bin of the detected values, whether that mode
is INTERIOR (a real peak) or sits at the left edge (background hidden below the
floor), and the low percentiles.

Run from the repo root:
    python multibatchcodes/plot_intensity_histograms_floor_peak.py
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MULTIBATCH_DIR = "MultiBatch_20to22"
RESULTS_DIR = "MultiBatchResults"
FLOOR = 3000.0
REP = ["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]
os.makedirs(RESULTS_DIR, exist_ok=True)

SURFACE, INK, INK_2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, BASELINE = "#e1e0d9", "#c3c2b7"
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a"]
plt.rcParams.update({
    "font.family": ["Segoe UI", "DejaVu Sans", "sans-serif"],
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "axes.edgecolor": BASELINE, "axes.linewidth": 0.8,
    "axes.labelcolor": INK_2, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK_2, "ytick.labelcolor": INK_2,
    "grid.color": GRID, "grid.linewidth": 0.8,
})

BINS = np.logspace(np.log10(FLOOR), np.log10(6e7), 70)

# %% Load
prot = {}
for path in sorted(glob.glob(os.path.join(MULTIBATCH_DIR, "*.csv"))):
    name = os.path.splitext(os.path.basename(path))[0]
    df = pd.read_csv(path, usecols=REP + ["ASMS_BATCH_NAME"])
    v = df[REP].to_numpy(dtype=float).ravel()
    v = v[np.isfinite(v)]
    prot[name] = {"batch": str(df["ASMS_BATCH_NAME"].dropna().iloc[0]),
                  "all": v, "det": v[v != FLOOR]}
batches = sorted({p["batch"] for p in prot.values()})
COLOR = dict(zip(batches, SLOTS))


def describe(det, all_v, per_decade=12):
    """Mode of log10(detected) from RAW bin counts, and whether it is interior.

    Deliberately NOT a KDE: a Gaussian kernel against the hard 3000 boundary leaks
    mass below the cutoff, depresses the density at the edge, and manufactures a
    spurious interior maximum. On this data that artifact reported a peak for every
    protein, including batches whose counts decay monotonically from the floor.
    Raw histogram counts have no boundary bias.

    interior=False means the observed distribution only DECAYS from the floor, i.e.
    the background peak lies BELOW 3000 and cannot be read off the observed data.
    """
    lg = np.log10(det)
    edges = np.arange(np.log10(FLOOR), lg.max() + 1 / per_decade, 1 / per_decade)
    cnt, _ = np.histogram(lg, bins=edges)
    j = int(np.argmax(cnt))
    return {"floor_frac": float((all_v == FLOOR).mean()),
            "mode": float(10 ** ((edges[j] + edges[j + 1]) / 2)),
            "interior": bool(j > 0),
            "rise_from_floor": float(cnt[j] / max(cnt[0], 1)),
            "peak_bin": j,
            "p1": float(np.percentile(det, 1)), "p5": float(np.percentile(det, 5)),
            "p25": float(np.percentile(det, 25)), "median": float(np.median(det))}


def draw(ax, det, all_v, color, title, stat, fs=9.5):
    ax.hist(det, bins=BINS, color=color, alpha=0.85, zorder=3)
    n_floor = int((all_v == FLOOR).sum())
    ax.bar([FLOOR * 0.93], [n_floor], width=FLOOR * 0.075, color=color, alpha=0.30,
           edgecolor=color, linewidth=1.1, zorder=3)          # the clamp point mass
    if stat["interior"]:
        ax.axvline(stat["mode"], color=INK, lw=1.1, ls=(0, (4, 2)), zorder=4)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(FLOOR * 0.845, 7e7)
    ax.xaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title(title, fontsize=fs + 1.5, color=INK, loc="left", pad=6)
    msg = (f"peak {stat['mode']:,.0f}" if stat["interior"]
           else "NO peak above floor\n(decays from 3000)")
    ax.text(0.975, 0.95, f"floor {stat['floor_frac']:.0%}\n{msg}",
            transform=ax.transAxes, ha="right", va="top", fontsize=fs,
            color=INK_2, linespacing=1.35)


# %% Numeric report
rows = []
for b in batches:
    members = [k for k, v in prot.items() if v["batch"] == b]
    det = np.concatenate([prot[k]["det"] for k in members])
    allv = np.concatenate([prot[k]["all"] for k in members])
    rows.append({"group": b, "kind": "batch", **describe(det, allv)})
for k, v in prot.items():
    rows.append({"group": k, "kind": "protein", "batch": v["batch"],
                 **describe(v["det"], v["all"])})
rep = pd.DataFrame(rows)

print("=" * 100)
print("IS THERE AN INTERIOR LOW-BAND PEAK ABOVE THE FLOOR?")
print("=" * 100)
show = ["group", "floor_frac", "mode", "interior", "rise_from_floor",
        "p1", "p5", "p25", "median"]
print(rep[rep.kind == "batch"][show].to_string(index=False,
                                               float_format=lambda x: f"{x:,.3f}"))
print()
print(rep[rep.kind == "protein"][["group", "batch"] + show[1:]]
      .to_string(index=False, float_format=lambda x: f"{x:,.3f}"))
n_int = int(rep[rep.kind == "protein"]["interior"].sum())
print(f"\nproteins with an interior peak: {n_int} of {(rep.kind == 'protein').sum()}")
print("rise_from_floor = density at the mode / density at the floor edge."
      "\n  <= 1.02  ->  monotone decay from the clamp: the background peak is BELOW"
      "\n              the floor and cannot be read off the observed histogram.")

_csv = os.path.join(RESULTS_DIR, "floor_peak_report.csv")
rep.to_csv(_csv, index=False)

# %% Figure 1 — per batch
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
for ax, b in zip(axes, batches):
    members = [k for k, v in prot.items() if v["batch"] == b]
    det = np.concatenate([prot[k]["det"] for k in members])
    allv = np.concatenate([prot[k]["all"] for k in members])
    st = rep[(rep.kind == "batch") & (rep.group == b)].iloc[0].to_dict()
    draw(ax, det, allv, COLOR[b], b, st)
    ax.set_xlabel("replicate intensity", fontsize=9.5)
axes[0].set_ylabel("count (log)", fontsize=9.5)
fig.suptitle("Intensity distribution per batch — is the low peak above the 3000 floor?",
             fontsize=14, color=INK, x=0.055, ha="left", y=0.985)
fig.text(0.055, 0.925,
         "faint bar at the far left = the point mass clamped onto 3000   ·   "
         "dashed line = modal bin of the detected values (raw counts, 12 bins/decade)",
         fontsize=9.5, color=INK_2)
fig.tight_layout(rect=(0, 0, 1, 0.88))
_p1 = os.path.join(RESULTS_DIR, "floor_peak_by_batch.png")
fig.savefig(_p1, dpi=200, facecolor=SURFACE)

# %% Figure 2 — small multiples per protein, grouped by batch
order = [k for b in batches for k in sorted(
    [k for k, v in prot.items() if v["batch"] == b])]
fig, axes = plt.subplots(4, 6, figsize=(19, 11), sharex=True, sharey=True)
for ax, name in zip(axes.ravel(), order):
    v = prot[name]
    st = rep[(rep.kind == "protein") & (rep.group == name)].iloc[0].to_dict()
    draw(ax, v["det"], v["all"], COLOR[v["batch"]],
         f"{name}   ·   {v['batch']}", st, fs=8)
for ax in axes[-1]:
    ax.set_xlabel("replicate intensity", fontsize=9)
for ax in axes[:, 0]:
    ax.set_ylabel("count (log)", fontsize=9)
fig.suptitle("Intensity distribution per protein — is the low peak above the 3000 floor?",
             fontsize=15, color=INK, x=0.045, ha="left", y=0.99)
fig.text(0.045, 0.955,
         "faint bar at the far left = the point mass clamped onto 3000   ·   "
         "dashed line = modal bin of the detected values (raw counts, 12 bins/decade)",
         fontsize=10, color=INK_2)
fig.tight_layout(rect=(0, 0, 1, 0.94))
_p2 = os.path.join(RESULTS_DIR, "floor_peak_by_protein.png")
fig.savefig(_p2, dpi=170, facecolor=SURFACE)

print(f"\nSaved {_p1}\nSaved {_p2}\nSaved {_csv}")


# %% Figure 3 — detected only, NO clamp bar, zoomed to the low decades
# The peak question lives just above the floor; on a 4-decade axis with the clamp
# spike drawn in, that region is unreadable. Here: linear y, no 3000s at all, and
# x stopping at 3e5, so the shape immediately above the floor is plain.
ZBINS = np.logspace(np.log10(FLOOR), np.log10(3e5), 46)
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
for ax, b in zip(axes, batches):
    members = [k for k, v in prot.items() if v["batch"] == b]
    det = np.concatenate([prot[k]["det"] for k in members])
    st = rep[(rep.kind == "batch") & (rep.group == b)].iloc[0].to_dict()
    ax.hist(det, bins=ZBINS, color=COLOR[b], alpha=0.85, zorder=3)
    if st["interior"]:
        ax.axvline(st["mode"], color=INK, lw=1.2, ls=(0, (4, 2)), zorder=4)
        ax.annotate(f"peak {st['mode']:,.0f}", (st["mode"], 0.965), xycoords=("data", "axes fraction"),
                    textcoords="offset points", xytext=(7, 0), ha="left", va="top",
                    fontsize=9.5, color=INK)
    else:
        ax.annotate("no peak — decays\nstraight from the floor", (FLOOR * 1.15, 0.955),
                    xycoords=("data", "axes fraction"), fontsize=9.5, color=INK,
                    ha="left", va="top", linespacing=1.35)
    ax.set_xscale("log")
    ax.set_xlim(FLOOR * 0.97, 3e5)
    ax.set_title(f"{b}   ·   floor {st['floor_frac']:.0%} (excluded here)",
                 fontsize=11, color=INK, loc="left", pad=6)
    ax.set_xlabel("replicate intensity", fontsize=9.5)
    ax.xaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
axes[0].set_ylabel("count", fontsize=9.5)
fig.suptitle("Detected values only, zoomed above the floor — where is the low peak?",
             fontsize=14, color=INK, x=0.055, ha="left", y=0.985)
fig.text(0.055, 0.925, "the 3000s are NOT drawn at all here   ·   linear y   ·   "
         "x stops at 3e5 so the region just above the floor is readable",
         fontsize=9.5, color=INK_2)
fig.tight_layout(rect=(0, 0, 1, 0.88))
_p3 = os.path.join(RESULTS_DIR, "floor_peak_zoom_no_clamp.png")
fig.savefig(_p3, dpi=200, facecolor=SURFACE)
print(f"Saved {_p3}")
