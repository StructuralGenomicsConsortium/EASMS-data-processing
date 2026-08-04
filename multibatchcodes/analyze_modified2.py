# -*- coding: utf-8 -*-
"""Analyse the median-based columns in MultiBatch_20to22_modified2/.

Reads the modified2 files and summarises / plots:
    1. per-batch intensity level        -> tests "is one batch hotter?"
    2. within- vs across-batch hits per batch
    3. enrichment_within vs enrichment_across (why across inflates)
    4. signal_flag distribution

Run from the repo root:  python multibatchcodes/analyze_modified2.py
"""

# %% Setup — load the modified2 columns
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DIR = "MultiBatch_20to22_modified2"
RESULTS_DIR = "MultiBatchResults"
FLOOR = 3000
os.makedirs(RESULTS_DIR, exist_ok=True)

cols = ["ASMS_BATCH_NAME", "target_median", "n_detected", "signal_flag",
        "enrichment_within", "enrichment_across",
        "pvalue_within", "pvalue_across", "label_within", "label_across"]
frames = []
for p in sorted(glob.glob(os.path.join(DIR, "*.csv"))):
    d = pd.read_csv(p, usecols=cols)
    d["protein"] = os.path.splitext(os.path.basename(p))[0]
    frames.append(d)
df = pd.concat(frames, ignore_index=True)
batches = sorted(df["ASMS_BATCH_NAME"].unique())
print(f"{len(df):,} rows | batches: {batches}")


# %% 1. Is one batch hotter? per-batch intensity + detection
print("\n=== per-batch signal level ===")
for b in batches:
    sub = df[df["ASMS_BATCH_NAME"] == b]
    detected = sub["target_median"][sub["target_median"] > FLOOR]
    print(f"{b}: detected {len(detected)/len(sub):5.1%} of rows | "
          f"median(target_median | detected) = {detected.median():,.0f} | "
          f"90th pct = {detected.quantile(0.9):,.0f}")

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
# (a) box of detected target_median per batch (log)
data = [df.loc[(df["ASMS_BATCH_NAME"] == b) & (df["target_median"] > FLOOR), "target_median"]
        for b in batches]
axes[0].boxplot(data, labels=batches, showfliers=False, whis=(5, 95))
axes[0].set_yscale("log")
axes[0].set_ylabel("target_median (detected rows only, log)")
axes[0].set_title("Per-batch signal level\n(higher batch -> across-enrichment inflates)")
# (b) fraction detected per batch
frac = [ (df["ASMS_BATCH_NAME"] == b).pipe(lambda m: (df.loc[m, "target_median"] > FLOOR).mean())
         for b in batches]
axes[1].bar(batches, frac)
axes[1].set_ylabel("fraction of rows detected (target_median > 3000)")
axes[1].set_title("Per-batch detection rate")
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "modified2_batch_level.png"), dpi=150)
plt.show()


# %% 2. Within- vs across-batch hits per batch
hits = df.groupby("ASMS_BATCH_NAME")[["label_within", "label_across"]].sum().loc[batches]
print("\n=== hits per batch ===")
print(hits.to_string())

x = np.arange(len(batches)); w = 0.38
fig, ax = plt.subplots(figsize=(8, 4.5))
b1 = ax.bar(x - w/2, hits["label_within"], w, label="label_within=1")
b2 = ax.bar(x + w/2, hits["label_across"], w, label="label_across=1")
ax.bar_label(b1, fmt="%d", fontsize=8); ax.bar_label(b2, fmt="%d", fontsize=8)
ax.set_xticks(x); ax.set_xticklabels(batches)
ax.set_ylabel("number of molecules called (hits)")
ax.set_title("Binder calls per batch: within vs across\n"
             "across is concentrated in the hot batch")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "modified2_hits_per_batch.png"), dpi=150)
plt.show()


# %% 3. enrichment_within vs enrichment_across, coloured by batch
d = df[(df["enrichment_within"] > 0) & (df["enrichment_across"] > 0)].copy()
d = d.sample(min(40000, len(d)), random_state=0)
fig, ax = plt.subplots(figsize=(7, 7))
for b in batches:
    s = d[d["ASMS_BATCH_NAME"] == b]
    ax.scatter(s["enrichment_within"], s["enrichment_across"], s=4, alpha=0.25, label=b)
lim = [d[["enrichment_within", "enrichment_across"]].min().min(),
       d[["enrichment_within", "enrichment_across"]].max().max()]
ax.plot(lim, lim, "k--", lw=1, label="y = x")
ax.axhline(5, color="r", ls=":", lw=1); ax.axvline(5, color="r", ls=":", lw=1)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("enrichment_within"); ax.set_ylabel("enrichment_across")
ax.set_title("Within vs across enrichment (red = 5x threshold)\n"
             "hot-batch points sit above y=x -> across over-calls")
ax.legend(markerscale=3)
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "modified2_within_vs_across.png"), dpi=150)
plt.show()


# %% 4. signal_flag distribution
order = ["3 strong", "2 strong", "1 strong", "3 moderate", "2 moderate", "1 moderate", "none"]
counts = df["signal_flag"].value_counts().reindex(order).fillna(0).astype(int)
print("\n=== signal_flag counts ===")
print(counts.to_string())
fig, ax = plt.subplots(figsize=(8, 4.5))
colors = ["#8c1515"]*3 + ["#d98a00"]*3 + ["0.7"]
bars = ax.bar(order, counts.values, color=colors)
ax.bar_label(bars, fmt="%d", fontsize=8)
ax.set_yscale("log")
ax.set_ylabel("molecules (log)")
ax.set_title("signal_flag distribution (strong = red, moderate = orange)")
plt.setp(ax.get_xticklabels(), rotation=25, ha="right")
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "modified2_signal_flag.png"), dpi=150)
plt.show()


# %% 5. Focused batch-22 view: shifted distribution + median bar
# Left: overlaid intensity distributions per batch (detected rows, log-x) — the
# sgcto_22 curve sits visibly to the right. Right: median detected intensity per
# batch as a bar, sgcto_22 highlighted, with value labels.
COLOR = {"sgcto_20": "0.6", "sgcto_21": "0.4", "sgcto_22": "#c1121f"}
det = {b: df.loc[(df["ASMS_BATCH_NAME"] == b) & (df["target_median"] > FLOOR), "target_median"]
       for b in batches}
meds = {b: det[b].median() for b in batches}

lo = min(v.min() for v in det.values())
hi = max(v.max() for v in det.values())
bins = np.logspace(np.log10(lo), np.log10(hi), 50)

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ax = axes[0]
for b in batches:
    ax.hist(det[b], bins=bins, density=True, histtype="step", lw=2,
            color=COLOR[b], label=f"{b} (median {meds[b]:,.0f})")
    ax.axvline(meds[b], color=COLOR[b], ls="--", lw=1.2)
ax.set_xscale("log")
ax.set_xlabel("target_median (detected rows, log)")
ax.set_ylabel("density")
ax.set_title("Intensity distribution per batch\nsgcto_22 is shifted higher")
ax.legend()

ax = axes[1]
bars = ax.bar(batches, [meds[b] for b in batches],
              color=[COLOR[b] for b in batches])
ax.bar_label(bars, labels=[f"{meds[b]:,.0f}" for b in batches], padding=3)
ax.set_ylabel("median target_median (detected rows)")
ax.set_title("Median detected intensity per batch\nsgcto_22 ~2x the others")
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "modified2_batch22_median.png"), dpi=150)
plt.show()


# %%
