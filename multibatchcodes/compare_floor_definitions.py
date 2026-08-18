# -*- coding: utf-8 -*-
"""Two different censoring rates: per-REPLICATE floor vs per-ROW median floor.

The plots (protein_intensity_range.py, plot_protein_median_by_batch.py) pool the
three replicate columns with .ravel(), so every individual REP value counts:
(3000, 4000, 6000) contributes one floored value out of three.

But the enrichment formula uses median(REP1..3) per row, and
median(3000, 4000, 6000) = 4000 -- no floor at all. Because 3000 is a hard clamp
(nothing in the data lies below it), median == 3000 requires >= 2 of the 3 reps
floored. So the row-level rate is strictly lower than the replicate-level rate.

This script reports both, per batch and per protein, plus how the protein
ranking changes when the box plot is built on per-row medians instead of pooled
replicates -- i.e. whether the "top 8 = sgcto_22" result survives the switch.

Run from the repo root:  python multibatchcodes/compare_floor_definitions.py
"""

import os
import glob
import numpy as np
import pandas as pd

MULTIBATCH_DIR = "MultiBatch_20to22"
RESULTS_DIR = "MultiBatchResults"
FLOOR = 3000
REP = ["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]
os.makedirs(RESULTS_DIR, exist_ok=True)

rows = []
for path in sorted(glob.glob(os.path.join(MULTIBATCH_DIR, "*.csv"))):
    name = os.path.splitext(os.path.basename(path))[0]
    df = pd.read_csv(path, usecols=REP + ["ASMS_BATCH_NAME"])
    batch = str(df["ASMS_BATCH_NAME"].dropna().iloc[0])

    R = df[REP].to_numpy(dtype=float)                    # n_rows x 3

    # --- definition 1: per individual replicate value (what the plots use)
    flat = R.ravel()
    flat = flat[np.isfinite(flat)]
    rep_floor = (flat == FLOOR).mean()
    med_rep_det = np.median(flat[flat != FLOOR])

    # --- definition 2: per row, on median(REP1..3)  (what enrichment uses)
    with np.errstate(all="ignore"):
        row_med = np.nanmedian(R, axis=1)
    row_med = row_med[np.isfinite(row_med)]
    row_floor = (row_med == FLOOR).mean()
    med_row_det = np.median(row_med[row_med != FLOOR])

    # how many reps are floored, per row -- explains the gap
    n_fl = (R == FLOOR).sum(axis=1)
    rows.append({
        "protein": name, "batch": batch,
        "floor_replicate": rep_floor,          # def 1
        "floor_rowmedian": row_floor,          # def 2
        "rows_0_floored": (n_fl == 0).mean(),
        "rows_1_floored": (n_fl == 1).mean(),
        "rows_2_floored": (n_fl == 2).mean(),
        "rows_3_floored": (n_fl == 3).mean(),
        "median_det_replicate": med_rep_det,
        "median_det_rowmedian": med_row_det,
    })

t = pd.DataFrame(rows)

print("=" * 92)
print("PER BATCH — the two censoring definitions")
print("=" * 92)
b = (t.groupby("batch")
     .agg(floor_replicate=("floor_replicate", "mean"),
          floor_rowmedian=("floor_rowmedian", "mean"),
          only_1_of_3=("rows_1_floored", "mean"),
          all_3_floored=("rows_3_floored", "mean"),
          med_det_rep=("median_det_replicate", "median"),
          med_det_row=("median_det_rowmedian", "median"))
     .sort_values("med_det_rep", ascending=False))
print((b.assign(**{c: (b[c] * 100).round(1) for c in
                   ["floor_replicate", "floor_rowmedian", "only_1_of_3", "all_3_floored"]})
       ).to_string(float_format=lambda x: f"{x:,.1f}"))
print("\nfloor_replicate / floor_rowmedian / only_1_of_3 / all_3_floored are %.")
print("The gap between the two definitions is exactly the 'only 1 of 3 floored'")
print("rows: they carry a 3000 into def 1 but their median is above the floor.")

# --- does the headline result survive the switch? ---------------------------
r1 = t.sort_values("median_det_replicate", ascending=False).reset_index(drop=True)
r2 = t.sort_values("median_det_rowmedian", ascending=False).reset_index(drop=True)
cmp = (r1[["protein", "batch"]].assign(rank_replicate=r1.index + 1)
       .merge(r2[["protein"]].assign(rank_rowmedian=r2.index + 1), on="protein"))
cmp["moved"] = cmp["rank_rowmedian"] - cmp["rank_replicate"]

print()
print("=" * 92)
print("PROTEIN RANKING — pooled replicates vs per-row medians")
print("=" * 92)
print(cmp.to_string(index=False))

top8_rep = set(r1.head(8)["protein"])
top8_row = set(r2.head(8)["protein"])
hot = b.index[0]
print(f"\ntop-8 identical under both definitions: {top8_rep == top8_row}")
print(f"top-8 all from {hot}: "
      f"{set(t.loc[t.protein.isin(top8_row), 'batch']) == {hot}}")
print(f"max rank movement: {cmp['moved'].abs().max()} position(s)")

_out = os.path.join(RESULTS_DIR, "floor_definitions_compare.csv")
t.to_csv(_out, index=False)
print(f"\nSaved {_out}")
