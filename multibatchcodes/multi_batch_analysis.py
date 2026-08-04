# -*- coding: utf-8 -*-
"""Multi-batch analysis scratchpad.

A collection of small, independent explorations of the MultiBatch/ data
(per-protein CSVs with the standard ASMS columns plus extra SPR gold-label
columns). Each section below is a self-contained '# %%' cell — in VS Code you
can run just one with Shift+Enter, or run the whole file top to bottom. Add new
explorations as new sections and keep the old ones around for reference.

Assumes it is run from the repo root (the default working directory in VS Code).
"""

# %% Setup — shared paths and imports
import os
import glob
import pandas as pd

# Absolute repo root so the paths resolve no matter the working directory (the
# VS Code interactive window's cwd is often NOT the repo root). Change this if
# you move the repo.
REPO_ROOT = r"d:\0000-UHN\EASMS-data-processing\EASMS-data-processing"
MULTIBATCH_DIR = os.path.join(REPO_ROOT, "MultiBatch")       # input: per-protein CSVs
RESULTS_DIR = os.path.join(REPO_ROOT, "MultiBatchResults")   # output: analysis results
os.makedirs(RESULTS_DIR, exist_ok=True)

# Quick diagnostic — confirm the folder and how many CSVs were found.
_files = sorted(glob.glob(os.path.join(MULTIBATCH_DIR, "*.csv")))
print(f"MultiBatch dir : {MULTIBATCH_DIR}")
print(f"CSV files found: {len(_files)}")


# %% 1. Row count per file -> MultiBatchResults/row_counts.csv
# Read every CSV in MultiBatch/ and record how many rows each has, then save
# the summary.
paths = sorted(glob.glob(os.path.join(MULTIBATCH_DIR, "*.csv")))
assert paths, f"No CSV files found in {MULTIBATCH_DIR} — check REPO_ROOT / that the data is there."

row_counts = []
for path in paths:
    n_rows = len(pd.read_csv(path))
    row_counts.append({"file": os.path.basename(path), "n_rows": n_rows})
    print(f"{os.path.basename(path):30} {n_rows:>8,} rows")

row_counts = pd.DataFrame(row_counts)
_out = os.path.join(RESULTS_DIR, "row_counts.csv")
row_counts.to_csv(_out, index=False)
print(f"\nSaved {_out}  ({len(row_counts)} files, {row_counts['n_rows'].sum():,} rows total)")


# %% 2. Replicate agreement — pool POS_INT_REP1/2/3 from every file
# Load just the three replicate columns from all files and keep rows where all
# three are present and > 0 (needed for the log axes below).
REP_COLS = ["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]
reps = pd.concat(
    [pd.read_csv(p, usecols=REP_COLS) for p in sorted(glob.glob(os.path.join(MULTIBATCH_DIR, "*.csv")))],
    ignore_index=True,
)
reps = reps.dropna(subset=REP_COLS)
reps = reps[(reps[REP_COLS] > 0).all(axis=1)]
print(f"{len(reps):,} measurements with all 3 positive replicates")


# %% 2a. Pairwise log-log scatter (rep vs rep) with y=x line + correlations
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import pearsonr, spearmanr

# Plot ALL measurements (no sampling). With millions of points this is slower
# and heavily overplotted, so markers are tiny and very transparent.
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for ax, (a, b) in zip(axes, combinations(REP_COLS, 2)):
    ax.scatter(reps[a], reps[b], s=1, alpha=0.05, edgecolors="none")
    ax.set_xscale("log")
    ax.set_yscale("log")
    lo = min(reps[a].min(), reps[b].min())
    hi = max(reps[a].max(), reps[b].max())
    ax.plot([lo, hi], [lo, hi], "r--", lw=1, label="y = x")
    pear = pearsonr(np.log10(reps[a]), np.log10(reps[b]))[0]
    spear = spearmanr(reps[a], reps[b])[0]
    ax.set_title(f"{a} vs {b}\nPearson(log)={pear:.3f}  Spearman={spear:.3f}")
    ax.set_xlabel(a)
    ax.set_ylabel(b)
    ax.legend(loc="upper left")

fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "replicate_scatter.png"), dpi=150)
plt.show()


# %% 2b. Coefficient of variation across the 3 reps (one summary of agreement)
# CV = std / mean per measurement. Low CV -> replicates agree.
cv = reps[REP_COLS].std(axis=1, ddof=1) / reps[REP_COLS].mean(axis=1)
print(f"median CV across replicates: {cv.median():.1%}")
print(f"fraction of measurements with CV > 30%: {(cv > 0.30).mean():.1%}")

fig, ax = plt.subplots(figsize=(7, 4))
ax.hist(cv.clip(upper=2), bins=60)
ax.axvline(cv.median(), color="r", ls="--", label=f"median = {cv.median():.1%}")
ax.set_xlabel("CV across REP1-3 (std/mean)")
ax.set_ylabel("number of measurements")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "replicate_cv_hist.png"), dpi=150)
plt.show()


# %% 3. Helper — replicate agreement (scatter + CV) for ANY subset of rows
# Reusable so we can run the same two plots on any slice (positives, negatives,
# one protein, ...). Pass a dataframe that has REP_COLS and a short `name` used
# in the titles and output filenames.
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import pearsonr, spearmanr


def plot_replicate_agreement(df, name, s=None, alpha=None):
    d = df.dropna(subset=REP_COLS)
    d = d[(d[REP_COLS] > 0).all(axis=1)]
    n = len(d)
    print(f"[{name}] {n:,} measurements with all 3 positive replicates")

    # Marker size/opacity scale with point count so any subset is visible:
    # few points -> bigger, more opaque; millions -> tiny, faint.
    if s is None:
        s = 1 if n > 200_000 else (3 if n > 20_000 else 10)
    if alpha is None:
        alpha = min(0.5, max(0.03, 5000 / max(n, 1)))

    # pairwise log-log scatter with y=x + correlations
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, (a, b) in zip(axes, combinations(REP_COLS, 2)):
        ax.scatter(d[a], d[b], s=s, alpha=alpha, edgecolors="none", zorder=1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        lo = min(d[a].min(), d[b].min())
        hi = max(d[a].max(), d[b].max())
        ax.plot([lo, hi], [lo, hi], "r--", lw=1, label="y = x")
        pear = pearsonr(np.log10(d[a]), np.log10(d[b]))[0]
        spear = spearmanr(d[a], d[b])[0]
        ax.set_title(f"{a} vs {b}\nPearson(log)={pear:.3f}  Spearman={spear:.3f}")
        ax.set_xlabel(a)
        ax.set_ylabel(b)
        ax.legend(loc="upper left")
    fig.suptitle(f"Replicate agreement — {name}")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, f"replicate_scatter_{name}.png"), dpi=150)
    plt.show()

    # CV histogram
    cv = d[REP_COLS].std(axis=1, ddof=1) / d[REP_COLS].mean(axis=1)
    print(f"[{name}] median CV: {cv.median():.1%} | fraction CV>30%: {(cv > 0.30).mean():.1%}")
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(cv.clip(upper=2), bins=60)
    ax.axvline(cv.median(), color="r", ls="--", label=f"median = {cv.median():.1%}")
    ax.set_xlabel("CV across REP1-3 (std/mean)")
    ax.set_ylabel("number of measurements")
    ax.set_title(f"Replicate CV — {name}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS_DIR, f"replicate_cv_hist_{name}.png"), dpi=150)
    plt.show()


# %% 4. Positives only (BINARY_LABEL == 1)
# Pool the REP columns + BINARY_LABEL from every file, keep the positive
# measurements, and run the same scatter + CV plots on them.
_cols = REP_COLS + ["BINARY_LABEL"]
labeled = pd.concat(
    [pd.read_csv(p, usecols=_cols) for p in sorted(glob.glob(os.path.join(MULTIBATCH_DIR, "*.csv")))],
    ignore_index=True,
)
# robust to BINARY_LABEL stored as 1 / "1" / True
is_pos = pd.to_numeric(labeled["BINARY_LABEL"], errors="coerce") == 1
positives = labeled[is_pos]
print(f"{len(positives):,} positive measurements (BINARY_LABEL == 1) "
      f"out of {len(labeled):,} total")
plot_replicate_agreement(positives, "positives")


# %% 5. All rows, excluding "flat 3000" measurements
# Some rows have all three replicates pinned at 3000 (a floor / placeholder, not
# a real measurement). Drop rows where all 3 reps == FLOOR and redo the plots.
# (Uses `reps` from cell 2 — run cell 2 first.)
FLOOR = 3000
all_at_floor = (reps[REP_COLS] == FLOOR).all(axis=1)
print(f"dropping {all_at_floor.sum():,} rows where all 3 reps == {FLOOR} "
      f"({all_at_floor.mean():.1%} of {len(reps):,})")
plot_replicate_agreement(reps[~all_at_floor], "all_no_floor")


# %% 6. Per-protein replicate noise — which protein (file) is problematic?
# For each file, compute the median CV across its measurements (excluding rows
# where all 3 reps are missing/<=0 or pinned at the 3000 floor) and plot them
# sorted, so noisy proteins stand out.
FLOOR = 3000
per_protein = []
for path in sorted(glob.glob(os.path.join(MULTIBATCH_DIR, "*.csv"))):
    d = pd.read_csv(path, usecols=REP_COLS).dropna(subset=REP_COLS)
    d = d[(d[REP_COLS] > 0).all(axis=1)]
    d = d[~(d[REP_COLS] == FLOOR).all(axis=1)]
    if d.empty:
        continue
    cv = d[REP_COLS].std(axis=1, ddof=1) / d[REP_COLS].mean(axis=1)
    per_protein.append({
        "protein": os.path.splitext(os.path.basename(path))[0],
        "n": len(d),
        "median_cv": cv.median(),
    })

per_protein = pd.DataFrame(per_protein).sort_values("median_cv", ascending=False)
per_protein.to_csv(os.path.join(RESULTS_DIR, "per_protein_cv.csv"), index=False)
overall_median = per_protein["median_cv"].median()
print(per_protein.to_string(index=False))

fig, ax = plt.subplots(figsize=(8, 0.35 * len(per_protein) + 1))
ax.barh(per_protein["protein"], per_protein["median_cv"])
ax.axvline(overall_median, color="r", ls="--", label=f"overall median = {overall_median:.1%}")
ax.invert_yaxis()  # highest CV (worst) on top
ax.set_xlabel("median CV across REP1-3")
ax.set_ylabel("protein (file)")
ax.set_title("Per-protein replicate noise (higher = worse)")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "per_protein_cv.png"), dpi=150)
plt.show()


# %% 7. Per-protein replicate noise — box plot of the CV distribution
# Like cell 6 but shows each protein's full CV distribution (median line, IQR
# box, 5-95% whiskers) instead of just the median, sorted by median CV.
FLOOR = 3000
cv_by_protein = {}
for path in sorted(glob.glob(os.path.join(MULTIBATCH_DIR, "*.csv"))):
    d = pd.read_csv(path, usecols=REP_COLS).dropna(subset=REP_COLS)
    d = d[(d[REP_COLS] > 0).all(axis=1)]
    d = d[~(d[REP_COLS] == FLOOR).all(axis=1)]
    if d.empty:
        continue
    cv = d[REP_COLS].std(axis=1, ddof=1) / d[REP_COLS].mean(axis=1)
    cv_by_protein[os.path.splitext(os.path.basename(path))[0]] = cv.values

# proteins ordered by median CV (worst ends up on top with horizontal boxes).
# `protein_order` is reused by cell 8 so the two box plots line up row-for-row.
protein_order = sorted(cv_by_protein, key=lambda k: np.median(cv_by_protein[k]))
data = [cv_by_protein[k] for k in protein_order]
overall_median = np.median(np.concatenate(list(cv_by_protein.values())))

fig, ax = plt.subplots(figsize=(8, 0.4 * len(protein_order) + 1))
ax.boxplot(data, vert=False, showfliers=False, whis=(5, 95))
ax.set_yticks(range(1, len(protein_order) + 1))
ax.set_yticklabels(protein_order)
ax.axvline(overall_median, color="r", ls="--", label=f"overall median = {overall_median:.1%}")
ax.set_xlabel("CV across REP1-3 (std/mean)")
ax.set_ylabel("protein (file)")
ax.set_title("Per protein coefficient of variation (excluding rows with all floor values)")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "per_protein_cv_box.png"), dpi=150)
plt.show()


# %% 8. Per-protein replicate CV — positives only (BINARY_LABEL == 1), box plot
# Same per-protein CV box plot as cell 7 but restricted to positive
# measurements (real binders).
FLOOR = 3000
cv_by_protein_pos = {}
for path in sorted(glob.glob(os.path.join(MULTIBATCH_DIR, "*.csv"))):
    d = pd.read_csv(path, usecols=REP_COLS + ["BINARY_LABEL"]).dropna(subset=REP_COLS)
    d = d[pd.to_numeric(d["BINARY_LABEL"], errors="coerce") == 1]
    d = d[(d[REP_COLS] > 0).all(axis=1)]
    d = d[~(d[REP_COLS] == FLOOR).all(axis=1)]
    if d.empty:
        continue
    cv = d[REP_COLS].std(axis=1, ddof=1) / d[REP_COLS].mean(axis=1)
    cv_by_protein_pos[os.path.splitext(os.path.basename(path))[0]] = cv.values

# Reuse the protein order from cell 7 so this plot lines up row-for-row with it
# (run cell 7 first). Falls back to sorting by its own median if not available.
try:
    order = [p for p in protein_order if p in cv_by_protein_pos]
except NameError:
    order = sorted(cv_by_protein_pos, key=lambda k: np.median(cv_by_protein_pos[k]))
data = [cv_by_protein_pos[k] for k in order]
overall_median = np.median(np.concatenate(list(cv_by_protein_pos.values())))

fig, ax = plt.subplots(figsize=(8, 0.4 * len(order) + 1))
ax.boxplot(data, vert=False, showfliers=False, whis=(5, 95))
ax.set_yticks(range(1, len(order) + 1))
ax.set_yticklabels(order)
ax.axvline(overall_median, color="r", ls="--", label=f"overall median = {overall_median:.1%}")
ax.set_xlabel("CV across REP1-3 (std/mean)")
ax.set_ylabel("protein (file)")
ax.set_title("Per-protein replicate CV — positives only (BINARY_LABEL=1)")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "per_protein_cv_box_positives.png"), dpi=150)
plt.show()


# %%


