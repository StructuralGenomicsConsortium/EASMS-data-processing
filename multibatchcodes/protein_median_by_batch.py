# -*- coding: utf-8 -*-
"""Is the per-protein intensity median a BATCH effect? — MultiBatch_20to22.

Reproduces the ordering of MultiBatchResults/protein_intensity_box.png
(protein_intensity_range.py: pool POS_INT_REP1/2/3, drop the 3000 floor, sort by
median) and annotates every protein with its ASMS_BATCH_NAME, so we can see
whether the high-median proteins at the top of that box plot all come from the
same batch.

Also reports, per batch:
    - median of the detected (>floor) replicate values
    - the floor fraction (share of replicate values sitting at 3000)

because a "hot" batch is expected to show BOTH a higher detected median AND a
lower floor fraction (late reading pushes borderline values above the floor).

Run from the repo root:  python multibatchcodes/protein_median_by_batch.py
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
per_batch_vals = {}          # batch -> list of replicate arrays (all values, floor kept)
for path in sorted(glob.glob(os.path.join(MULTIBATCH_DIR, "*.csv"))):
    name = os.path.splitext(os.path.basename(path))[0]
    df = pd.read_csv(path, usecols=REP + ["ASMS_BATCH_NAME"])

    v_all = df[REP].to_numpy(dtype=float).ravel()
    v_all = v_all[np.isfinite(v_all)]
    v_det = v_all[v_all != FLOOR]                    # same filter as the box plot

    batches = sorted(df["ASMS_BATCH_NAME"].dropna().unique().tolist())
    rows.append({
        "protein": name,
        "batch": " | ".join(map(str, batches)),
        "n_batches": len(batches),
        "n_rows": len(df),
        "n_detected": len(v_det),
        "floor_frac": 1 - len(v_det) / len(v_all) if len(v_all) else np.nan,
        "median_detected": np.median(v_det) if len(v_det) else np.nan,
    })
    for b in batches:
        per_batch_vals.setdefault(str(b), []).append(
            df.loc[df["ASMS_BATCH_NAME"] == b, REP].to_numpy(dtype=float).ravel())

tab = pd.DataFrame(rows).sort_values("median_detected", ascending=False)

print("=" * 78)
print("PER PROTEIN — ordered by median of detected intensity (top = top of plot)")
print("=" * 78)
with pd.option_context("display.float_format", lambda x: f"{x:,.3f}"):
    print(tab.to_string(index=False))

print()
print("=" * 78)
print("PER BATCH")
print("=" * 78)
brows = []
for b, arrs in sorted(per_batch_vals.items()):
    v_all = np.concatenate(arrs)
    v_all = v_all[np.isfinite(v_all)]
    v_det = v_all[v_all != FLOOR]
    brows.append({
        "batch": b,
        "n_proteins": int((tab["batch"] == b).sum()),
        "n_values": len(v_all),
        "floor_frac": 1 - len(v_det) / len(v_all),
        "median_detected": np.median(v_det),
        "mean_detected": v_det.mean(),
    })
btab = pd.DataFrame(brows).sort_values("median_detected", ascending=False)
with pd.option_context("display.float_format", lambda x: f"{x:,.3f}"):
    print(btab.to_string(index=False))

# Rank check: are the top-median proteins concentrated in one batch?
print()
print("=" * 78)
print("RANK BY BATCH  (rank 1 = highest median protein)")
print("=" * 78)
tab = tab.reset_index(drop=True)
tab["rank"] = tab.index + 1
for b, g in tab.groupby("batch"):
    print(f"{b:16} ranks {sorted(g['rank'].tolist())}"
          f"   mean rank {g['rank'].mean():5.1f}")

_out = os.path.join(RESULTS_DIR, "protein_median_by_batch.csv")
tab.to_csv(_out, index=False)
print(f"\nSaved {_out}")
