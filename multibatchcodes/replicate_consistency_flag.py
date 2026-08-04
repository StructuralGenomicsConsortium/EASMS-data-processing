# -*- coding: utf-8 -*-
"""Replicate-consistency QC flag — is a lone signal a real read or a bad replicate?

MultiBatch_20to22. 3000 = "detected but too low to quantify" (a floored low read),
not missing data. A molecule with ONE detected replicate and TWO floored (3000)
reads is anomalous: a genuinely strong analyte should not be near-zero in 2 of 3
injections.

Hypothesis (per-batch, replicates assumed positionally aligned across the 8
proteins of a batch): if the SAME compound shows this "1 detected + 2 floor" shape
across several proteins of a batch AND it is always the SAME replicate index that
carries the signal, then the two floored reads are a systematic replicate problem
-> the detected value is probably real and should be KEPT. If instead the signal
lands on DIFFERENT replicate indices across proteins, the lone high is the suspect
-> likely a spurious spike to DISCARD.

This script:
  1. tests replicate alignment empirically (do lone detections land on the same
     replicate index across a batch's proteins more than chance?),
  2. flags each (compound, batch): consistent_same_rep / inconsistent / isolated,
  3. plots example heatmaps (proteins x replicates) + a batch summary.

Run from the repo root:  python multibatchcodes/replicate_consistency_flag.py
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

DIR = "MultiBatch_20to22"
RESULTS_DIR = "MultiBatchResults"
FLOOR = 3000
STRONG_MIN = 1e6
REP = ["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]
CONSISTENT_FRAC = 2 / 3      # >= this share of lone proteins on one rep index -> consistent
os.makedirs(RESULTS_DIR, exist_ok=True)

# ---- load ----
frames = []
for p in sorted(glob.glob(os.path.join(DIR, "*.csv"))):
    d = pd.read_csv(p, usecols=["COMPOUND_ID", "ASMS_BATCH_NAME"] + REP)
    d["protein"] = os.path.splitext(os.path.basename(p))[0]
    frames.append(d)
df = pd.concat(frames, ignore_index=True)
reps = df[REP].to_numpy(float)
det = reps > FLOOR
n_det = det.sum(axis=1)
df["n_detected"] = n_det
df["is_lone"] = n_det == 1
df["lone_rep_index"] = np.where(n_det == 1, det.argmax(axis=1) + 1, 0)  # 1..3, else 0
df["max_intensity"] = reps.max(axis=1)
df["lone_is_strong"] = df["is_lone"] & (df["max_intensity"] >= STRONG_MIN)
batches = sorted(df["ASMS_BATCH_NAME"].unique())
print(f"{len(df):,} rows | lone detections (n_detected==1): {int(df['is_lone'].sum()):,} "
      f"| of those strong: {int(df['lone_is_strong'].sum()):,}")

# ---- per (compound, batch) lone-detection consistency ----
lone = df[df["is_lone"]]
grp = lone.groupby(["COMPOUND_ID", "ASMS_BATCH_NAME"])["lone_rep_index"]
stat = pd.DataFrame({
    "lone_count": grp.size(),
    "modal_index": grp.agg(lambda s: int(s.value_counts().index[0])),
    "modal_count": grp.agg(lambda s: int(s.value_counts().iloc[0])),
})
stat["consistency"] = stat["modal_count"] / stat["lone_count"]

def flag_of(r):
    if r["lone_count"] <= 1:
        return "isolated"
    return "consistent_same_rep" if (r["modal_count"] >= 2 and r["consistency"] >= CONSISTENT_FRAC) \
        else "inconsistent"
stat["flag"] = stat.apply(flag_of, axis=1)

df = df.merge(stat[["lone_count", "modal_index", "modal_count", "consistency", "flag"]],
              left_on=["COMPOUND_ID", "ASMS_BATCH_NAME"], right_index=True, how="left")
df["flag"] = df["flag"].fillna("none")        # compound not lone-detecting in this batch

# ---- 1. empirical replicate-alignment test (cleanest at lone_count == 2) ----
print("\n=== replicate-alignment test (compound-batches with exactly 2 lone proteins) ===")
print("under random replicate assignment, P(both on same index) = 1/3 = 0.333")
two = stat[stat["lone_count"] == 2]
rows = []
for b in batches:
    sub = two.xs(b, level="ASMS_BATCH_NAME") if b in two.index.get_level_values(1) else two.iloc[0:0]
    same = (sub["modal_count"] == 2).mean() if len(sub) else np.nan
    rows.append((b, len(sub), same))
    print(f"  {b}: n={len(sub):>5} pairs | same-index rate = {same:.1%}")
align = pd.DataFrame(rows, columns=["batch", "n_pairs", "same_rate"])

fig, ax = plt.subplots(figsize=(7, 4.5))
bars = ax.bar(align["batch"], align["same_rate"], color="#2a7fb8")
ax.bar_label(bars, labels=[f"{v:.0%}" for v in align["same_rate"]], padding=3)
ax.axhline(1/3, color="r", ls="--", label="chance (1/3)")
ax.set_ylabel("fraction of compound-pairs on the SAME replicate index")
ax.set_title("Do 2 lone detections of a compound share the same replicate?\n"
             "(well above 1/3 => replicates aligned & the effect is systematic)")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "replicate_alignment_test.png"), dpi=150)

# ---- 2. summary of flags + rescue counts ----
cb = stat.reset_index()   # one row per (compound, batch) that has >=1 lone detection
print("\n=== compound-batch flags (>=1 lone detection) ===")
print(cb["flag"].value_counts().to_string())

strong_lone = df[df["lone_is_strong"]]
print(f"\n=== the lone-STRONG rows ({len(strong_lone):,}) by flag ===")
print(strong_lone["flag"].value_counts().to_string())
print("  consistent_same_rep -> rescue candidates (keep the strong)")
print("  inconsistent        -> spike suspects (discard)")

sflag = strong_lone["flag"].value_counts().reindex(
    ["consistent_same_rep", "inconsistent", "isolated"]).fillna(0).astype(int)
fig, ax = plt.subplots(figsize=(7, 4.5))
colors = ["#2ca25f", "#c1121f", "0.6"]
bars = ax.bar(sflag.index, sflag.values, color=colors)
ax.bar_label(bars, fmt="%d")
ax.set_ylabel("lone-strong molecules")
ax.set_title("Lone-strong (1 strong + 2 floor) molecules by replicate-consistency flag\n"
             "green = keep candidate, red = spike suspect")
plt.setp(ax.get_xticklabels(), rotation=15, ha="right")
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "replicate_flag_summary.png"), dpi=150)

# ---- 3. example heatmaps (proteins x replicates) ----
def heatmap_matrix(cid, batch):
    sub = df[(df["COMPOUND_ID"] == cid) & (df["ASMS_BATCH_NAME"] == batch)].sort_values("protein")
    return sub["protein"].tolist(), sub[REP].to_numpy(float)

def pick(flag, n, strong=True, min_lone=3):
    cand = stat[(stat["flag"] == flag) & (stat["lone_count"] >= min_lone)].reset_index()
    if strong:
        s = strong_lone.groupby(["COMPOUND_ID", "ASMS_BATCH_NAME"]).size().reset_index()
        cand = cand.merge(s, on=["COMPOUND_ID", "ASMS_BATCH_NAME"])
    return cand.sort_values("lone_count", ascending=False).head(n)

examples = [("consistent_same_rep", "#2ca25f"), ("inconsistent", "#c1121f")]
picks = []
for flag, _ in examples:
    p = pick(flag, 3, strong=True, min_lone=3)
    if len(p) < 3:
        p = pick(flag, 3, strong=False, min_lone=(2 if flag == "inconsistent" else 3))
    picks.append(p)

fig, axes = plt.subplots(2, 3, figsize=(14, 8))
for r, (grp_df, (flag, col)) in enumerate(zip(picks, examples)):
    for c in range(3):
        ax = axes[r, c]
        if c < len(grp_df):
            row = grp_df.iloc[c]
            prots, mat = heatmap_matrix(row["COMPOUND_ID"], row["ASMS_BATCH_NAME"])
            im = ax.imshow(np.clip(mat, FLOOR, None), aspect="auto",
                           norm=LogNorm(vmin=FLOOR, vmax=max(mat.max(), FLOOR*2)),
                           cmap="viridis")
            ax.set_xticks(range(3)); ax.set_xticklabels(["REP1", "REP2", "REP3"])
            ax.set_yticks(range(len(prots))); ax.set_yticklabels(prots, fontsize=7)
            ax.set_title(f"{row['COMPOUND_ID']} [{row['ASMS_BATCH_NAME']}]\n"
                         f"{flag} (modal REP{int(row['modal_index'])}, "
                         f"{int(row['modal_count'])}/{int(row['lone_count'])})", fontsize=8, color=col)
            fig.colorbar(im, ax=ax, fraction=0.046)
        else:
            ax.axis("off")
fig.suptitle("Example compounds — intensity across a batch's proteins (rows) x replicates (cols)\n"
             "top: consistent vertical stripe (keep) | bottom: inconsistent (suspect)")
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "replicate_flag_heatmaps.png"), dpi=150)
print("\nsaved: replicate_alignment_test.png, replicate_flag_summary.png, replicate_flag_heatmaps.png")
