# -*- coding: utf-8 -*-
"""Which statistic of the detected values actually SCALES with the batch factor?

The substitution value for a censored observation is only useful if it moves with the
batch: a value that stays put is exactly what makes the batch effect leak through the
within-batch ratio. So the test for any candidate estimator is not "is it a sensible
low value" but "does its ratio between a hot and a cold batch match the true batch
factor".

Every candidate here is computed on the DETECTED values (3000s removed). All of them
are biased by the floor in the SAME direction -- truncation from below compresses any
low statistic toward 3000, pulling its hot/cold ratio toward 1. So the estimator with
the LARGEST hot/cold ratio is the least floor-pinned, and that ratio is a lower bound
on the true batch factor.

Candidates: mean of the bottom 10% / 25%, low percentiles, the modal bin (peak),
the median, the mean.

Run from the repo root:  python multibatchcodes/compare_floor_estimators.py
"""

import os
import glob
import numpy as np
import pandas as pd

MULTIBATCH_DIR = "MultiBatch_20to22"
RESULTS_DIR = "MultiBatchResults"
FLOOR = 3000.0
REP = ["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]
os.makedirs(RESULTS_DIR, exist_ok=True)

det_by_batch, floor_by_batch = {}, {}
for path in sorted(glob.glob(os.path.join(MULTIBATCH_DIR, "*.csv"))):
    df = pd.read_csv(path, usecols=REP + ["ASMS_BATCH_NAME"])
    b = str(df["ASMS_BATCH_NAME"].dropna().iloc[0])
    v = df[REP].to_numpy(dtype=float).ravel()
    v = v[np.isfinite(v)]
    det_by_batch.setdefault(b, []).append(v[v != FLOOR])
    floor_by_batch.setdefault(b, []).append(v)


def modal_bin(det, per_decade=12):
    lg = np.log10(det)
    edges = np.arange(np.log10(FLOOR), lg.max() + 1 / per_decade, 1 / per_decade)
    cnt, _ = np.histogram(lg, bins=edges)
    j = int(np.argmax(cnt))
    return 10 ** ((edges[j] + edges[j + 1]) / 2), j > 0


def bottom_mean(det, frac):
    """Mean of the lowest `frac` of the detected values."""
    return det[det <= np.quantile(det, frac)].mean()


rows = []
for b in sorted(det_by_batch):
    det = np.concatenate(det_by_batch[b])
    allv = np.concatenate(floor_by_batch[b])
    peak, interior = modal_bin(det)
    rows.append({
        "batch": b,
        "floor%": 100 * (allv == FLOOR).mean(),
        "bottom10%mean": bottom_mean(det, 0.10),
        "bottom25%mean": bottom_mean(det, 0.25),
        "p5": np.percentile(det, 5),
        "p10": np.percentile(det, 10),
        "peak": peak if interior else np.nan,
        "median": np.median(det),
        "mean": det.mean(),
    })
t = pd.DataFrame(rows).set_index("batch")

print("=" * 96)
print("CANDIDATE SUBSTITUTION VALUES, per batch (all on DETECTED values)")
print("=" * 96)
print(t.to_string(float_format=lambda x: f"{x:,.0f}"))

hot, mid, cold = "sgcto_22", "sgcto_20", "sgcto_21"
cols = ["bottom10%mean", "bottom25%mean", "p5", "p10", "peak", "median", "mean"]
r = pd.DataFrame({
    "hot/cold  (22/21)": [t.loc[hot, c] / t.loc[cold, c] for c in cols],
    "hot/mid   (22/20)": [t.loc[hot, c] / t.loc[mid, c] for c in cols],
}, index=cols).sort_values("hot/cold  (22/21)", ascending=False)

print()
print("=" * 96)
print("DOES IT SCALE?  ratio between batches — bigger = less floor-pinned")
print("=" * 96)
print(r.to_string(float_format=lambda x: f"{x:,.2f}"))
print("""
A perfect estimator would reproduce the true batch factor. Every candidate is
compressed toward 1.00 by the floor, so the LARGEST ratio here is the least-biased
and is a lower bound on the truth. Ratios near 1.00 mean the statistic barely moves
between a 29%-censored and an 82%-censored batch -> it would NOT cancel the batch
effect if substituted.
""")
t.to_csv(os.path.join(RESULTS_DIR, "floor_estimator_compare.csv"))
print("Saved", os.path.join(RESULTS_DIR, "floor_estimator_compare.csv"))
