# -*- coding: utf-8 -*-
"""Estimate the substitution value for a censored (3000) observation, per batch.

The number we want is  E[T | T < 3000]  -- the expected TRUE intensity of an
observation the instrument could not detect. It is below 3000 by definition and it
must differ per batch (a hot batch's undetected values sit nearer 3000 than a cold
batch's), which is the property that makes the within-batch ratio cancel.

Method: censored maximum likelihood on BINNED counts. Detected values are binned in
log space (400 bins); each bin contributes n_i * log[F(b_i) - F(a_i)], and the n
censored values contribute n_cens * log F(3000) -- all we know is that they fell
below it. Binning makes each likelihood evaluation ~400 CDF calls instead of ~200,000
density calls; with this many points the accuracy cost is negligible.

Two models, because the shape assumption is the weak point:

    1LN   one log-normal            -- simple, but the histograms are visibly a
                                       mixture (background peak + long binder tail),
                                       so this is expected to be biased
    2LN   two log-normals           -- background component + binder component, the
                                       shape the data actually shows

Validation (the part that decides whether any of this is usable): sgcto_22 is only
29% censored and its background peak is directly visible at ~12,651. So we CRIPPLE
it -- re-censor at the quantile that makes it 82% censored, matching sgcto_21 --
refit, and check whether it recovers the uncrippled answer. A method that cannot
recover a known answer at 82% censoring cannot be trusted on sgcto_21.

Run from the repo root:  python multibatchcodes/fit_censored_floor.py
"""

import os
import glob
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.stats import norm

MULTIBATCH_DIR = "MultiBatch_20to22"
RESULTS_DIR = "MultiBatchResults"
FLOOR = 3000.0
NBIN = 400
REP = ["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]
os.makedirs(RESULTS_DIR, exist_ok=True)

# %% Load per batch
raw = {}
for path in sorted(glob.glob(os.path.join(MULTIBATCH_DIR, "*.csv"))):
    df = pd.read_csv(path, usecols=REP + ["ASMS_BATCH_NAME"])
    b = str(df["ASMS_BATCH_NAME"].dropna().iloc[0])
    v = df[REP].to_numpy(dtype=float).ravel()
    raw.setdefault(b, []).append(v[np.isfinite(v)])
raw = {b: np.concatenate(a) for b, a in raw.items()}

TINY = 1e-300


def binned(v, floor):
    """Log-space bin counts of the detected part + the censored count."""
    det = v[v > floor]
    y = np.log(det)
    edges = np.linspace(np.log(floor), y.max() + 1e-9, NBIN + 1)
    cnt, _ = np.histogram(y, bins=edges)
    keep = cnt > 0
    return {"lo": edges[:-1][keep], "hi": edges[1:][keep], "n": cnt[keep],
            "n_cens": int((v <= floor).sum()), "lf": np.log(floor),
            "n_det": len(det)}


def cdf_1(x, p):
    return norm.cdf((x - p[0]) / np.exp(p[1]))


def cdf_2(x, p):
    w = 1 / (1 + np.exp(-p[0]))
    return (w * norm.cdf((x - p[1]) / np.exp(p[2]))
            + (1 - w) * norm.cdf((x - p[3]) / np.exp(p[4])))


def make_nll(d, cdf):
    def nll(p):
        try:
            mass = cdf(d["hi"], p) - cdf(d["lo"], p)
            tail = cdf(d["lf"], p)
        except FloatingPointError:
            return 1e18
        if not np.all(np.isfinite(mass)) or not np.isfinite(tail):
            return 1e18
        ll = np.sum(d["n"] * np.log(np.maximum(mass, TINY)))
        ll += d["n_cens"] * np.log(max(tail, TINY))
        return -ll if np.isfinite(ll) else 1e18
    return nll


def fit(d, which):
    lf = d["lf"]
    if which == 1:
        starts = [[lf + dm, np.log(s)] for dm in (-2, -1, 0, 1)
                  for s in (0.5, 1.0, 2.0)]
        cdf, bounds = cdf_1, [(lf - 12, lf + 12), (np.log(0.02), np.log(8))]
    else:
        starts = [[np.log(w / (1 - w)), lf + d1, np.log(s1), lf + d2, np.log(s2)]
                  for w in (0.5, 0.85) for d1 in (-1.5, -0.3, 0.5)
                  for d2 in (1.5, 3.0) for s1 in (0.4, 1.0) for s2 in (1.0,)]
        cdf = cdf_2
        bounds = [(-6, 6), (lf - 12, lf + 12), (np.log(0.02), np.log(8)),
                  (lf - 12, lf + 14), (np.log(0.05), np.log(8))]
    nll = make_nll(d, cdf)
    best = None
    for p0 in starts:
        r = minimize(nll, p0, method="L-BFGS-B", bounds=bounds)
        if best is None or r.fun < best.fun:
            best = r
    return best


def trunc_mean_ln(mu, s, lc):
    """E[X | X<c] for log X ~ N(mu,s); lc = log c. Natural-log units."""
    return np.exp(mu + s ** 2 / 2
                  + norm.logcdf((lc - mu - s ** 2) / s)
                  - norm.logcdf((lc - mu) / s))


def analyse(v, floor, label):
    d = binned(v, floor)
    r1, r2 = fit(d, 1), fit(d, 2)
    mu1, s1 = r1.x[0], np.exp(r1.x[1])

    w = 1 / (1 + np.exp(-r2.x[0]))
    a, sa, b, sb = r2.x[1], np.exp(r2.x[2]), r2.x[3], np.exp(r2.x[4])
    if a > b:                                   # component A = the lower = background
        w, a, sa, b, sb = 1 - w, b, sb, a, sa
    lf3 = np.log(FLOOR)
    la = np.log(max(w, TINY)) + norm.logcdf((lf3 - a) / sa)
    lb = np.log(max(1 - w, TINY)) + norm.logcdf((lf3 - b) / sb)
    z = logsumexp([la, lb])
    sub2 = (np.exp(la - z) * trunc_mean_ln(a, sa, lf3)
            + np.exp(lb - z) * trunc_mean_ln(b, sb, lf3))

    return {"group": label, "floor_used": floor,
            "cens%": 100 * d["n_cens"] / (d["n_cens"] + d["n_det"]),
            "1LN_mode": np.exp(mu1 - s1 ** 2),
            "1LN_sub": trunc_mean_ln(mu1, s1, lf3),
            "2LN_w_bg": w, "2LN_bg_mode": np.exp(a - sa ** 2),
            "2LN_bg_median": np.exp(a), "2LN_sub": sub2,
            "d_nll(2LN-1LN)": r2.fun - r1.fun}


# %% 1. Each real batch
t = pd.DataFrame([analyse(raw[b], FLOOR, b) for b in sorted(raw)]).set_index("group")
print("=" * 108)
print("SUBSTITUTION VALUE  E[T | T < 3000]  per batch    (the *_sub columns)")
print("=" * 108)
print(t.to_string(float_format=lambda x: f"{x:,.2f}"))
print()
for c in ("1LN_sub", "2LN_sub", "1LN_mode", "2LN_bg_mode"):
    lo, hi = t.loc["sgcto_21", c], t.loc["sgcto_22", c]
    print(f"   {c:<14} hot/cold ratio 22/21 = {hi / lo:>8,.2f}")

# %% 2. Validation
v22 = raw["sgcto_22"]
target = float((raw["sgcto_21"] <= FLOOR).mean())
c82 = float(np.quantile(v22, target))
print()
print("=" * 108)
print(f"VALIDATION — re-censor sgcto_22 at {c82:,.0f} so it becomes "
      f"{100 * target:.1f}% censored (as sgcto_21 is), then refit")
print("=" * 108)
truth = analyse(v22, FLOOR, "sgcto_22 TRUE (29% cens)")
crip = analyse(np.maximum(v22, c82), c82, f"sgcto_22 CRIPPLED ({100*target:.0f}% cens)")
print(pd.DataFrame([truth, crip]).set_index("group")
      .to_string(float_format=lambda x: f"{x:,.2f}"))
print(f"\ndirectly observed background peak of sgcto_22 = 12,651")
for k in ("1LN_mode", "2LN_bg_mode"):
    print(f"   {k:<12}  true {truth[k]:>12,.0f}   crippled {crip[k]:>12,.0f}"
          f"   ratio {crip[k]/truth[k]:>6,.2f}x")

pd.concat([t, pd.DataFrame([truth, crip]).set_index("group")]).to_csv(
    os.path.join(RESULTS_DIR, "censored_floor_fit.csv"))
print("\nSaved", os.path.join(RESULTS_DIR, "censored_floor_fit.csv"))
