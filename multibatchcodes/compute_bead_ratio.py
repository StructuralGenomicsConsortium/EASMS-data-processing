# -*- coding: utf-8 -*-
"""Add bead-specificity columns -> MultiBatch_20to22_modified5/.

modified5 = every column of MultiBatch_20to22_modified4/ (raw median enrichment,
median-normalized, percentile-normalized) PLUS bead-only specificity columns.

A compound should NOT be dropped just because it shows some bead-only signal; we
compare its protein signal to its bead-only signal:

    bead_signal = max( median(5-bead reps), median(10-bead reps) )   # conservative
    bead_ratio  = target_median (on protein) / bead_signal

Matched by COMPOUND_ID (= "SGC ID for Component" in the bead file).

Thresholds are NOT yet decided, so the raw bead_signal and bead_ratio are stored
(apply any cutoff later); convenience flags are provided at 5x and 10x using the
provisional "detected on beads = bead_signal > 3000" basis. bead_signal is there
to recompute a stricter basis (e.g. >= 10,000) whenever the threshold is set.

Run from the repo root:  python multibatchcodes/compute_bead_ratio.py
"""

import os
import glob
import numpy as np
import pandas as pd

SRC_DIR = "MultiBatch_20to22_modified4"
OUT_DIR = "MultiBatch_20to22_modified5"
BEADS = os.path.join("MultiBatchResults", "beads_clean.xlsx")
FLOOR = 3000
C5 = ["5Beads_V1P0113_1", "5Beads_V1P0113_2", "5Beads_V1P0113_3"]
C10 = ["10Beads_V1P0113_1", "10Beads_V1P0113_2", "10Beads_V1P0113_3"]

# --- bead signal per compound (conservative: max of the two condition medians) ---
b = pd.read_excel(BEADS)
b["bead_signal"] = np.maximum(b[C5].median(axis=1), b[C10].median(axis=1))
bead = (b[["SGC ID for Component", "bead_signal"]]
        .rename(columns={"SGC ID for Component": "COMPOUND_ID"})
        .dropna(subset=["COMPOUND_ID"]).drop_duplicates("COMPOUND_ID"))
print(f"bead compounds: {len(bead):,} | detected on beads (>{FLOOR}): "
      f"{int((bead['bead_signal'] > FLOOR).sum()):,}")

os.makedirs(OUT_DIR, exist_ok=True)
paths = sorted(glob.glob(os.path.join(SRC_DIR, "*.csv")))
n_rows = matched = 0
hits_total = hits_r5 = hits_r10 = 0
for p in paths:
    df = pd.read_csv(p)
    df = df.merge(bead, on="COMPOUND_ID", how="left")
    df["bead_detected"] = (df["bead_signal"] > FLOOR).astype(int)   # NaN -> 0
    df["bead_ratio"] = df["target_median"] / df["bead_signal"]
    df["bead_binder_flag_r5"] = ((df["bead_signal"] > FLOOR) & (df["bead_ratio"] < 5)).astype(int)
    df["bead_binder_flag_r10"] = ((df["bead_signal"] > FLOOR) & (df["bead_ratio"] < 10)).astype(int)
    df.to_csv(os.path.join(OUT_DIR, os.path.basename(p)), index=False)

    n_rows += len(df)
    matched += int(df["bead_signal"].notna().sum())
    h = df["label_within"] == 1
    hits_total += int(h.sum())
    hits_r5 += int((h & (df["bead_binder_flag_r5"] == 1)).sum())
    hits_r10 += int((h & (df["bead_binder_flag_r10"] == 1)).sum())

print(f"\nSaved {len(paths)} files to {OUT_DIR}/")
print(f"rows: {n_rows:,} | matched to a bead entry: {matched/n_rows:.1%}")
print(f"within-batch hits: {hits_total:,} | bead_binder_r5: {hits_r5:,} "
      f"({hits_r5/hits_total:.1%}) | bead_binder_r10: {hits_r10:,} ({hits_r10/hits_total:.1%})")
print("\nnew columns: bead_signal, bead_detected, bead_ratio, "
      "bead_binder_flag_r5, bead_binder_flag_r10")
