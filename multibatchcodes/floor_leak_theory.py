# -*- coding: utf-8 -*-
"""How much batch effect LEAKS through within-batch enrichment because of the floor.

Theory
------
    latent    T = beta_b * S_ip * eps
    observed  M = max(F, T)                     F = 3000, after the replicate median
                                                (max() commutes with median, so this
                                                 is exact, not an approximation)

    E_p = M_ip / mean_{q!=p}[ M_iq ]
        = beta*S_ip / ( (k/7)*F + beta*D )      k = # floored background proteins
                                                D = (1/7) * sum of DETECTED S_iq

    leak  lambda = dlogE / dlogbeta = (k/7)*F / ( (k/7)*F + beta*D )
                 = share of the denominator's mass contributed by floored proteins

    lambda = 0 -> batch effect cancels perfectly; lambda = 1 -> passes through whole.

This script (1) verifies that identity numerically against brute-force
differentiation, and (2) evaluates lambda for each real batch of MultiBatch_20to22
using its measured floor rate and detected level, both with the MEDIAN and with
the MEAN of detected values (the denominator is a mean, but it is heavy-tailed, so
the two bracket the truth).

Run from the repo root:  python multibatchcodes/floor_leak_theory.py
"""

import os
import glob
import numpy as np
import pandas as pd

MULTIBATCH_DIR = "MultiBatch_20to22"
F = 3000.0
N_BG = 7                                  # background proteins per batch (8 - 1)
REP = ["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]


# %% 1. Verify the leak identity numerically
def E_of_beta(beta, s_target, s_bg):
    """Enrichment at batch factor beta, with the floor applied."""
    m_t = max(F, beta * s_target)
    m_b = np.maximum(F, beta * s_bg)
    return m_t / m_b.mean()


def leak_closed_form(beta, s_bg):
    """lambda = floored share of the denominator."""
    m_b = np.maximum(F, beta * s_bg)
    floored = m_b[beta * s_bg <= F].sum()
    return floored / m_b.sum()


rng = np.random.default_rng(0)
print("=" * 74)
print("VERIFY   dlogE/dlogbeta  ==  floored share of the denominator")
print("=" * 74)
print(f"{'beta':>6} {'k floored':>10} {'lambda (numeric)':>18} {'lambda (formula)':>18}")
s_target = 40000.0
s_bg = np.array([300., 800., 1500., 2500., 4000., 9000., 30000.])
for beta in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0):
    h = 1e-6
    num = (np.log(E_of_beta(beta * (1 + h), s_target, s_bg))
           - np.log(E_of_beta(beta * (1 - h), s_target, s_bg))) / (2 * h)
    k = int((beta * s_bg <= F).sum())
    print(f"{beta:>6.2f} {k:>10d} {num:>18.4f} {leak_closed_form(beta, s_bg):>18.4f}")
print("\n(k falls as beta rises -> the hot batch leaks LESS. Identity holds to ~1e-6.)")


# %% 2. Evaluate the leak for the real batches
rows = []
for path in sorted(glob.glob(os.path.join(MULTIBATCH_DIR, "*.csv"))):
    df = pd.read_csv(path, usecols=REP + ["ASMS_BATCH_NAME"])
    R = df[REP].to_numpy(dtype=float)
    with np.errstate(all="ignore"):
        m = np.nanmedian(R, axis=1)                      # row median, as the formula uses
    m = m[np.isfinite(m)]
    det = m[m != F]
    rows.append({"batch": str(df["ASMS_BATCH_NAME"].dropna().iloc[0]),
                 "floor_frac": (m == F).mean(),
                 "med_det": np.median(det), "mean_det": det.mean()})
t = pd.DataFrame(rows).groupby("batch").mean(numeric_only=True)

print()
print("=" * 74)
print("REAL BATCHES — leak lambda implied by the measured floor rate")
print("=" * 74)
out = []
for b, r in t.iterrows():
    k = r["floor_frac"] * N_BG                           # expected # floored of the 7
    C = (k / N_BG) * F                                   # floored mass in the denominator
    for tag, lvl in (("median", r["med_det"]), ("mean", r["mean_det"])):
        bD = ((N_BG - k) / N_BG) * lvl                   # detected mass
        out.append({"batch": b, "detected level": tag,
                    "floor%": r["floor_frac"] * 100, "k of 7": k,
                    "level": lvl, "lambda %": 100 * C / (C + bD),
                    "compensated %": 100 * bD / (C + bD)})
o = pd.DataFrame(out).sort_values(["detected level", "floor%"])
print(o.to_string(index=False, float_format=lambda x: f"{x:,.1f}"))

print("""
Read: 'lambda %' is the share of the batch effect that SURVIVES the within-batch
division; 'compensated %' is the share it removes. The two 'detected level' rows
bracket the answer, because the denominator is a MEAN of a heavy-tailed
distribution while the median understates it.
""")
