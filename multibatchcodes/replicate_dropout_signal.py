# -*- coding: utf-8 -*-
"""When replicates drop to the floor, how strong is the detected signal?

MultiBatch_20to22 protein data. Each molecule row has three replicate intensities
(POS_INT_REP1/2/3). A replicate is "low" when it sits at the 3000 floor (not
detected). We ask:

    Scenario A - exactly ONE replicate is low  -> classify the other TWO
    Scenario B - exactly TWO replicates are low -> classify the ONE that remains

Each detected (above-floor) replicate is binned by the standard MS intensity
table:

    weak            10^3 - 10^4      (here: 3000 < v < 10,000)
    moderate        10^4 - 10^6
    strong          10^6 - 10^7
    very abundant   10^7 +           (sometimes higher)

(The noise band 10^2-10^3 cannot occur: the floor is 3000.)

The point: if a molecule shows real signal in only one or two replicates, is that
detected signal typically weak (likely noise-ish / irreproducible) or genuinely
strong (a real hit the other replicates simply missed)?

Run from the repo root:  python multibatchcodes/replicate_dropout_signal.py
"""

# %% Setup — pool all replicate triples, tag floor drop-outs
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

MULTIBATCH_DIR = "MultiBatch_20to22"
RESULTS_DIR = "MultiBatchResults"
FLOOR = 3000
REP = ["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]
os.makedirs(RESULTS_DIR, exist_ok=True)

# standard-table bins for ABOVE-floor values -> 0 weak, 1 moderate, 2 strong, 3 abundant
CAT_EDGES = [1e4, 1e6, 1e7]
CATS = ["weak", "moderate", "strong", "very abundant"]
# x-tick labels that keep the category name AND its standard-table intensity range
CAT_RANGES = ["$10^3$–$10^4$", "$10^4$–$10^6$", "$10^6$–$10^7$", "$10^7$+"]
CAT_TICKLABELS = [f"{c}\n{r}" for c, r in zip(CATS, CAT_RANGES)]

reps = pd.concat(
    [pd.read_csv(p, usecols=REP) for p in sorted(glob.glob(os.path.join(MULTIBATCH_DIR, "*.csv")))],
    ignore_index=True,
).to_numpy(dtype=float)

is_low = reps == FLOOR
n_low = is_low.sum(axis=1)          # 0..3 floor replicates per molecule
print(f"{len(reps):,} molecule rows")
for k in range(4):
    print(f"  {int((n_low == k).sum()):>7,} rows with {k} replicate(s) at floor "
          f"({(n_low == k).mean():.1%})")


def classify(vals):
    """Bin above-floor intensities into CATS indices (0..3)."""
    return np.digitize(vals, CAT_EDGES)


def category_counts(vals):
    """Counts per category for a flat array of above-floor intensities."""
    idx = classify(vals)
    return np.array([(idx == c).sum() for c in range(len(CATS))])


# %% 1. Scenario A — exactly ONE replicate low: classify the other two
maskA = n_low == 1
detectedA = reps[maskA][~is_low[maskA]]     # the 2 above-floor values per row, flattened
countsA = category_counts(detectedA)
print(f"\n[A] {int(maskA.sum()):,} rows with exactly 1 low replicate "
      f"-> {len(detectedA):,} detected replicates")
for c, name in enumerate(CATS):
    print(f"    {name:14} {countsA[c]:>7,} ({countsA[c] / countsA.sum():.1%})")


# %% 2. Scenario B — exactly TWO replicates low: classify the one that remains
maskB = n_low == 2
detectedB = reps[maskB][~is_low[maskB]]     # the single above-floor value per row
countsB = category_counts(detectedB)
print(f"\n[B] {int(maskB.sum()):,} rows with exactly 2 low replicates "
      f"-> {len(detectedB):,} detected replicates")
for c, name in enumerate(CATS):
    print(f"    {name:14} {countsB[c]:>7,} ({countsB[c] / max(countsB.sum(),1):.1%})")


# %% 3. Summary table (saved) + grouped bar chart
summary = pd.DataFrame({
    "category": CATS,
    "A_one_low_count": countsA,
    "A_one_low_pct": countsA / countsA.sum(),
    "B_two_low_count": countsB,
    "B_two_low_pct": countsB / max(countsB.sum(), 1),
})
with pd.option_context("display.float_format", lambda x: f"{x:,.3f}"):
    print("\n" + summary.to_string(index=False))
_out = os.path.join(RESULTS_DIR, "replicate_dropout_signal.csv")
summary.to_csv(_out, index=False)
print(f"Saved {_out}")

x = np.arange(len(CATS))
w = 0.38
fig, ax = plt.subplots(figsize=(8, 4.5))
barsA = ax.bar(x - w / 2, summary["A_one_low_pct"], w,
               label=f"1 low, 2 detected (n={int(countsA.sum()):,})")
barsB = ax.bar(x + w / 2, summary["B_two_low_pct"], w,
               label=f"2 low, 1 detected (n={int(countsB.sum()):,})")
# annotate every bar with its raw count so the tiny strong/abundant bars
# (fraction ~0.001 or 0) are still readable, not invisible at this scale
ax.bar_label(barsA, labels=[f"{c:,}" for c in countsA], padding=2, fontsize=8)
ax.bar_label(barsB, labels=[f"{c:,}" for c in countsB], padding=2, fontsize=8)
ax.set_xticks(x)
ax.set_xticklabels(CAT_TICKLABELS)
ax.set_ylabel("fraction of detected replicates")
ax.set_xlabel("standard-table signal of the detected replicate(s)")
ax.set_title("Signal of the detected replicate(s) when the others drop to the 3000 floor\n"
             "MultiBatch_20to22")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "replicate_dropout_signal.png"), dpi=150)
plt.show()


# %% 4. Same data on a log scale — raw counts so the tiny strong bars are visible
# Log y makes the ~0.1% strong category readable next to the large weak/moderate
# bars. Plotted as absolute counts (log of a count is well-defined); "very
# abundant" is 0 so it has no bar on a log axis, but its 0 label still shows.
fig, ax = plt.subplots(figsize=(8, 4.5))
barsA = ax.bar(x - w / 2, countsA, w,
               label=f"1 low, 2 detected (n={int(countsA.sum()):,})")
barsB = ax.bar(x + w / 2, countsB, w,
               label=f"2 low, 1 detected (n={int(countsB.sum()):,})")
ax.bar_label(barsA, labels=[f"{c:,}" for c in countsA], padding=2, fontsize=8)
ax.bar_label(barsB, labels=[f"{c:,}" for c in countsB], padding=2, fontsize=8)
ax.set_yscale("log")
ax.set_xticks(x)
ax.set_xticklabels(CAT_TICKLABELS)
ax.set_ylabel("number of detected replicates (log scale)")
ax.set_xlabel("standard-table signal of the detected replicate(s)")
ax.set_title("Signal of the detected replicate(s) when the others drop to the 3000 floor\n"
             "MultiBatch_20to22 (log scale)")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "replicate_dropout_signal_log.png"), dpi=150)
plt.show()


# %%
