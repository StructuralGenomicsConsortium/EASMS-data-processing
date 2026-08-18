# -*- coding: utf-8 -*-
"""Does a generally hotter batch look MORE enriched, just because the floor is 3000?

Short answer: yes for the rows that matter, and the amount is computable.

The within-batch fold divides one number by others measured in the same batch, so a
pure gain factor beta cancels EXACTLY -- but only for the part of the ratio that was
really measured. The 3000 placeholder does not scale with beta, so every value stuck
at the floor breaks the cancellation:

    observed   M = max(3000, beta * T)                T = latent true intensity
    fold       E = M_target / mean(M_background)

    sensitivity   s = dlogE / dlogbeta = lambda_bg - f_target

        lambda_bg  floored share of the background's mass  (mean of the other 7)
        f_target   floored share of the target row's own mean intensity

    s =  0  -> the batch effect cancels, the fold is honest
    s =  1  -> target measured, background all at the floor: the apparent fold is
               DIRECTLY PROPORTIONAL to how hot the batch ran (fold = target/3000)
    s = -1  -> target at the floor, background measured: the fold is squashed to 1

s is exact as a local sensitivity: while no value crosses the floor, a 1% hotter
batch reports a fold s% higher. (The closed form and its numerical check live in
floor_leak_theory.py, which treats the lambda_bg half.)

Four panels:
  A  two REAL compounds of one batch that look equally enriched, put through the
     same 2x gain: background measured -> fold unchanged; background at the floor
     -> fold x2. Pure arithmetic, nothing simulated.
  B  the controlled experiment: take the HOTTEST batch, treat its recorded values as
     the truth, and ask what the very same samples would have reported had the
     instrument run colder -- multiply everything by beta, re-apply max(3000, .),
     recompute. Turning the gain DOWN is the honest direction: every value it moves
     is one we actually measured, so no latent below 3000 has to be invented.
     (Treating the batch's own 3000s as if they were true 3000s makes this a LOWER
     bound: the real latents sit below 3000 and would floor even more readily.)
  C  mean s by fold band -- the leak is ~0 for the bulk and climbs to ~0.4 exactly
     in the hit range, because a hit is a compound its background never saw.
  D  the spread of s among hits: for a sixth of them s > 0.9, i.e. their fold is
     essentially target/3000 and moves one-for-one with batch gain.

What this does NOT show: a raw correlation between batch hotness and hit rate across
the 17 batches (Spearman 0.20, p = 0.43). Each batch holds different proteins and
compounds, so real biology swamps the artifact in a cross-batch comparison -- which
is why the controlled rescaling in panel B is the right way to see it.

Run from the repo root:  python multibatchcodes/floor_gain_demo.py
"""

import os
import glob

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.patches import Patch

DATA_DIR = "MultiBatch_20to36"
RESULTS_DIR = "MultiBatchResults"
FLOOR = 3000.0
N_PROT = 8                                   # proteins per batch -> N = 7 background
HIT_FOLD = 5.0                               # the pipeline's pass_fold cutoff
REP = ["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]
KEYS = ["ASMS_BATCH_NAME", "COMPOUND_ID"]
os.makedirs(RESULTS_DIR, exist_ok=True)

# --- dataviz tokens (light surface), as used by the other MultiBatchResults plots
SURFACE  = "#fcfcfb"
INK      = "#0b0b0b"
INK_2    = "#52514e"
MUTED    = "#898781"
GRID     = "#e1e0d9"
BASELINE = "#c3c2b7"
BLUE, ORANGE = "#2a78d6", "#eb6834"

plt.rcParams.update({
    "font.family": ["Segoe UI", "DejaVu Sans", "sans-serif"],
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "axes.edgecolor": BASELINE, "axes.linewidth": 0.8,
    "axes.labelcolor": INK_2, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK_2, "ytick.labelcolor": INK_2,
    "grid.color": GRID, "grid.linewidth": 0.8, "grid.linestyle": "-",
})


# %% 1. Load every protein file as one (batch, compound)-sorted table
def load():
    parts = []
    for path in sorted(glob.glob(os.path.join(DATA_DIR, "*.csv"))):
        d = pd.read_csv(path, usecols=KEYS + REP, low_memory=False)
        d["PROTEIN"] = os.path.splitext(os.path.basename(path))[0]
        parts.append(d)
    big = pd.concat(parts, ignore_index=True)
    # mergesort = stable, so the 8 proteins of a group stay in file order
    big = big.sort_values(KEYS + ["PROTEIN"], kind="mergesort").reset_index(drop=True)
    if len(big) % N_PROT:
        raise SystemExit(f"{len(big)} rows is not a multiple of {N_PROT}")
    first = big[KEYS].to_numpy()[::N_PROT]                   # first row of each block
    if not (np.repeat(first, N_PROT, axis=0) == big[KEYS].to_numpy()).all():
        raise SystemExit("a (batch, compound) group does not have exactly 8 rows")
    return big


big = load()
G = len(big) // N_PROT
R = big[REP].to_numpy(dtype=float)
print(f"{DATA_DIR}: {len(big):,} rows | {big['PROTEIN'].nunique()} proteins | "
      f"{big['ASMS_BATCH_NAME'].nunique()} batches | {G:,} (batch, compound) groups | "
      f"{(R == FLOOR).mean():.1%} of all replicate values sit at the {FLOOR:,.0f} floor")

# per row: A = mean of the valid replicates; floorA = the part of A contributed by
# replicates at the 3000 placeholder (NaN never equals FLOOR, so this is safe)
valid = np.isfinite(R)
n_valid = valid.sum(axis=1)
den = np.where(n_valid > 0, n_valid, 1)
A = np.where(n_valid > 0, np.nansum(R, axis=1) / den, np.nan)
floorA = np.where(n_valid > 0, FLOOR * (R == FLOOR).sum(axis=1) / den, np.nan)

A8, F8 = A.reshape(G, N_PROT), floorA.reshape(G, N_PROT)
SA = np.nansum(A8, axis=1, keepdims=True)
SF = np.nansum(F8, axis=1, keepdims=True)
K = np.isfinite(A8).sum(axis=1, keepdims=True)
bg_sum = SA - A8                                   # sum(A_i) over the other proteins
with np.errstate(divide="ignore", invalid="ignore"):
    fold = A8 * (K - 1.0) / bg_sum                 # article eq. 1, as the pipeline does it
    lam_bg = (SF - F8) / bg_sum                    # floored share of the background
    f_tgt = F8 / A8                                # floored share of the target itself
fold = np.where((K - 1 <= 0) | (bg_sum <= 0), np.nan, fold)
lam_bg = np.where(bg_sum <= 0, np.nan, lam_bg)

flat = pd.DataFrame({
    "batch": big["ASMS_BATCH_NAME"].to_numpy(),
    "compound": big["COMPOUND_ID"].to_numpy(),
    "protein": big["PROTEIN"].to_numpy(),
    "A": A, "fold": fold.ravel(),
    "lam_bg": lam_bg.ravel(), "f_tgt": f_tgt.ravel(),
}).dropna(subset=["fold", "lam_bg", "f_tgt"])
flat["s"] = flat["lam_bg"] - flat["f_tgt"]         # net gain sensitivity dlogE/dlogbeta
hits = flat[flat["fold"] >= HIT_FOLD]
print(f"usable rows: {len(flat):,} | hits (fold >= {HIT_FOLD:.0f}): {len(hits):,} | "
      f"mean s all rows {flat['s'].mean():+.3f} | mean s of hits {hits['s'].mean():+.3f} | "
      f"hits with s > 0.9: {(hits['s'] > 0.9).mean():.1%}")

# batch "hotness" = median of the DETECTED replicate values (using all values would
# just report how censored the batch is)
med_det = (pd.DataFrame({"batch": np.repeat(big["ASMS_BATCH_NAME"].to_numpy(), len(REP)),
                         "v": R.ravel()})
           .query("v > @FLOOR").groupby("batch")["v"].median())


# %% 2. Panel A — two real compounds that look equally enriched, same 2x gain
hot = str(med_det.idxmax())
band = flat[(flat["batch"] == hot) & flat["fold"].between(4.0, 8.0) & (flat["f_tgt"] == 0)]
floored = band[band["lam_bg"] > 0.999]                     # background all undetected
measured = band[band["lam_bg"] < 0.001]                    # background all detected
if floored.empty or measured.empty:
    raise SystemExit(f"no matched pair in batch {hot}")
# match the two as closely as possible on the fold they report today
floored = floored.iloc[[(floored["fold"] - measured["fold"].median()).abs().argmin()]]
measured = measured.iloc[[(measured["fold"] - floored["fold"].iloc[0]).abs().argmin()]]

ex = []
for tag, row in (("background at the 3000 floor\n(all 7 undetected)", floored.iloc[0]),
                 ("background really measured\n(all 7 detected)", measured.iloc[0])):
    t, b = row["A"], row["A"] / row["fold"]                 # target, mean background
    b2 = b if row["lam_bg"] > 0.5 else 2 * b                # the floor does not scale
    ex.append({"tag": tag, "protein": row["protein"], "compound": row["compound"],
               "target": t, "bg": b, "fold_1x": t / b, "fold_2x": 2 * t / b2})
ex = pd.DataFrame(ex)
print(f"\nPanel A — batch {hot}:")
print(ex.to_string(index=False, float_format=lambda x: f"{x:,.1f}"))


# %% 3. Panel B — run the hottest batch's own samples at a colder gain
cold_ratio = float(med_det.min() / med_det.max())          # coldest/hottest real batch
sel = (big["ASMS_BATCH_NAME"].to_numpy() == hot).reshape(G, N_PROT).all(axis=1)
Rg = R.reshape(G, N_PROT, len(REP))[sel]
print(f"\nPanel B: batch {hot} (hottest, median detected {med_det.max():,.0f}) | "
      f"{Rg.shape[0]:,} compounds | the coldest batch ran {1 / cold_ratio:.2f}x colder")

betas = np.unique(np.concatenate([np.logspace(np.log10(0.25), 0.0, 25),
                                  [0.25, 0.5, cold_ratio, 1.0]]))
curve = []
for beta in betas:
    X = np.maximum(FLOOR, beta * Rg)                        # the floor does not scale
    Ab = X.mean(axis=2)
    S = Ab.sum(axis=1, keepdims=True)
    fb = (Ab * (N_PROT - 1.0) / (S - Ab)).ravel()
    curve.append({"beta": beta, "hits": int((fb >= HIT_FOLD).sum()),
                  "top1pct_fold": float(fb[fb >= np.percentile(fb, 99)].mean()),
                  "max_fold": float(fb.max()), "floor_share": float((X == FLOOR).mean())})
curve = pd.DataFrame(curve)
base = curve.loc[np.isclose(curve["beta"], 1.0)].iloc[0]
curve["hits_idx"] = curve["hits"] / base["hits"]
curve["fold_idx"] = curve["top1pct_fold"] / base["top1pct_fold"]
at_cold = curve.loc[np.isclose(curve["beta"], cold_ratio)].iloc[0]
print(f"the same samples at the coldest batch's gain ({cold_ratio:.2f}x): "
      f"hits {int(base['hits'])} -> {int(at_cold['hits'])} ({at_cold['hits_idx'] - 1:+.0%}), "
      f"top-1% fold {base['top1pct_fold']:.2f} -> {at_cold['top1pct_fold']:.2f} "
      f"({at_cold['fold_idx'] - 1:+.0%})")


# %% 4. Panels C/D — where the leak sits, measured on the data as it stands
EDGES = [-np.inf, 0.999, 1.001, 1.5, 2, 3, HIT_FOLD, 10, 30, np.inf]
NAMES = ["< 1", "= 1\n(nothing\ndetected)", "1 - 1.5", "1.5 - 2", "2 - 3", "3 - 5",
         "5 - 10", "10 - 30", "> 30"]
flat["band"] = pd.cut(flat["fold"], EDGES, labels=NAMES)
bandtab = (flat.groupby("band", observed=True)
           .agg(n=("s", "size"), s=("s", "mean"), lam_bg=("lam_bg", "mean"),
                f_tgt=("f_tgt", "mean")).reset_index())
print("\nmean gain sensitivity s by fold band:")
print(bandtab.to_string(index=False, float_format=lambda x: f"{x:,.3f}"))

bat = (flat.groupby("batch").agg(n=("s", "size")).join(
       hits.groupby("batch").agg(n_hits=("s", "size"), s_hits=("s", "mean")))
       .reset_index())
bat["median_detected"] = bat["batch"].map(med_det)
bat["hit_rate_pct"] = 100 * bat["n_hits"] / bat["n"]
bat = bat.sort_values("median_detected").reset_index(drop=True)
print("\nper batch:")
print(bat.to_string(index=False, float_format=lambda x: f"{x:,.3f}"))
rho = bat["median_detected"].corr(bat["hit_rate_pct"], method="spearman")
print(f"\nSpearman(batch hotness, hit rate) = {rho:.2f} — the artifact is NOT visible as "
      f"a raw cross-batch\ncorrelation, because each batch holds different proteins; "
      f"that is what panel B controls for.")

_csv = os.path.join(RESULTS_DIR, "floor_gain_demo_batches.csv")
bat.to_csv(_csv, index=False)


# %% 5. Figure
fig = plt.figure(figsize=(13.8, 9.6))
gs = fig.add_gridspec(2, 2, wspace=0.26, hspace=0.60,
                      left=0.072, right=0.978, top=0.845, bottom=0.155)
axA, axB = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
axC, axD = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])

# --- A: same 2x gain on two real compounds -----------------------------------
x = np.arange(len(ex))
w = 0.34
for off, col, key, lab in ((-w / 2 - 0.012, BLUE, "fold_1x", "the batch as it ran"),
                           (+w / 2 + 0.012, ORANGE, "fold_2x", "same sample, 2x hotter batch")):
    axA.bar(x + off, ex[key], width=w, color=col, zorder=3, label=lab)
    for xi, v in zip(x + off, ex[key]):
        axA.text(xi, v + 0.25, f"{v:,.1f}", ha="center", va="bottom",
                 fontsize=9.5, color=INK_2)
axA.axhline(HIT_FOLD, color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=2)
axA.text(1.48, HIT_FOLD + 0.2, f"hit cutoff, fold {HIT_FOLD:.0f}", fontsize=8.5,
         color=INK_2, va="bottom", ha="right")
axA.set_xticks(x)
axA.set_xticklabels(ex["tag"], fontsize=9)
axA.set_xlim(-0.55, 1.55)
axA.set_ylim(0, ex[["fold_1x", "fold_2x"]].to_numpy().max() * 1.25)
axA.yaxis.grid(True, zorder=0)
axA.set_axisbelow(True)
for s in ("top", "right"):
    axA.spines[s].set_visible(False)
axA.tick_params(axis="x", length=0)
axA.set_ylabel("apparent enrichment fold", fontsize=9.5)
axA.set_title("A  Two real compounds that look equally enriched today",
              fontsize=11, color=INK, loc="left", pad=8)
axA.legend(loc="upper right", fontsize=9, frameon=False, labelcolor=INK_2,
           handlelength=1.1, borderaxespad=0.3)
axA.text(0.5, -0.30, f"batch {hot}  ·  {ex['protein'].iloc[0]} / {ex['compound'].iloc[0]}  "
                     f"and  {ex['protein'].iloc[1]} / {ex['compound'].iloc[1]}\n"
                     f"(a background that is undetected but not far below 3000 would land "
                     f"between the two bars)",
         transform=axA.transAxes, ha="center", fontsize=8.5, color=MUTED)

# --- B: the same samples run at a colder gain --------------------------------
for key, col, lab in (("hits_idx", BLUE, f"compounds called a hit (fold $\\geq$ {HIT_FOLD:.0f})"),
                      ("fold_idx", ORANGE, "size of the top 1% of folds")):
    axB.plot(curve["beta"], curve[key] * 100, color=col, lw=2.0, zorder=4, label=lab)
axB.axhline(100, color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1)
axB.axvline(cold_ratio, color=BASELINE, lw=1.0, zorder=1)
axB.text(cold_ratio * 1.04, 96, f"the coldest batch ran\n{1 / cold_ratio:.1f}x colder",
         fontsize=8.5, color=INK_2, va="top")
for key, col, dy, va in (("hits_idx", BLUE, -12, "top"), ("fold_idx", ORANGE, 12, "bottom")):
    axB.scatter([cold_ratio, 1.0], [at_cold[key] * 100, 100], s=42, color=col, zorder=5,
                edgecolor=SURFACE, linewidth=1.6)
    axB.annotate(f"{at_cold[key] * 100:.0f}%", (cold_ratio, at_cold[key] * 100),
                 textcoords="offset points", xytext=(-8, dy), ha="right", va=va,
                 fontsize=9.5, color=INK)
axB.set_xscale("log")
axB.set_xticks([0.25, 0.35, 0.5, 0.7, 1.0])
axB.set_xticklabels(["0.25x", "0.35x", "0.5x", "0.7x", "1x\n(as it ran)"], fontsize=9)
axB.xaxis.set_minor_formatter(NullFormatter())
axB.set_xlim(0.24, 1.07)
axB.set_ylim(40, 106)
axB.grid(True, which="major", zorder=0)
axB.set_axisbelow(True)
for s in ("top", "right"):
    axB.spines[s].set_visible(False)
axB.set_xlabel(f"gain applied to every value of batch {hot}", fontsize=9.5)
axB.set_ylabel("% of what the batch actually reported", fontsize=9.5)
axB.set_title("B  The same samples, run colder", fontsize=11, color=INK, loc="left", pad=8)
axB.legend(loc="lower right", fontsize=9, frameon=False, labelcolor=INK_2,
           handlelength=1.4, borderaxespad=0.5)
axB.text(0.5, -0.30, f"{Rg.shape[0]:,} compounds of the hottest batch, its own values taken "
                     f"as the truth. Turning the gain down\nonly moves numbers we really "
                     f"measured, so this is a lower bound on the effect.",
         transform=axB.transAxes, ha="center", fontsize=8.5, color=MUTED)

# --- C: mean sensitivity by fold band ----------------------------------------
xc = np.arange(len(bandtab))
cols = [BLUE if v >= 0 else ORANGE for v in bandtab["s"]]
axC.bar(xc, bandtab["s"], width=0.62, color=cols, zorder=3)
for xi, v, n in zip(xc, bandtab["s"], bandtab["n"]):
    axC.text(xi, v + (0.022 if v >= 0 else -0.022), f"{v:+.2f}", ha="center",
             va="bottom" if v >= 0 else "top", fontsize=8.5, color=INK_2)
axC.axhline(0, color=BASELINE, lw=1.0, zorder=2)
axC.set_xticks(xc)
axC.set_xticklabels(bandtab["band"], fontsize=8.5)
axC.set_ylim(-0.26, 0.52)
axC.yaxis.grid(True, zorder=0)
axC.set_axisbelow(True)
for s in ("top", "right"):
    axC.spines[s].set_visible(False)
axC.tick_params(axis="x", length=0)
axC.set_xlabel("apparent fold reported today", fontsize=9.5)
axC.set_ylabel("mean sensitivity  $s$ = dlog(fold) / dlog(gain)", fontsize=9.5)
axC.set_title("C  The leak climbs exactly into the hit range",
              fontsize=11, color=INK, loc="left", pad=8)
axC.text(0.5, -0.34, "$s$ = 0: the batch effect cancels  ·  $s$ = 1: the apparent fold "
                     "is proportional to batch gain\nnegative $s$: the target itself is at "
                     "the floor, so a cold batch squashes the fold toward 1",
         transform=axC.transAxes, ha="center", fontsize=8.5, color=MUTED)

# --- D: spread of s among the hits -------------------------------------------
edges = np.linspace(-1, 1, 41)
for pop, col, lab in ((flat["s"], MUTED, f"all {len(flat):,} rows"),
                      (hits["s"], BLUE, f"the {len(hits):,} hits (fold $\\geq$ {HIT_FOLD:.0f})")):
    axD.hist(pop, bins=edges, weights=np.full(len(pop), 100 / len(pop)),
             color=col, alpha=0.55 if col is MUTED else 0.9, zorder=3, label=lab)
share = (hits["s"] > 0.9).mean()
axD.axvline(0.9, color=ORANGE, lw=1.4, zorder=4)
axD.annotate(f"{share:.0%} of hits sit above 0.9 —\ntheir fold is essentially "
             f"target / 3000\nand tracks batch gain one-for-one",
             xy=(0.9, 7.5), xytext=(-14, 0), textcoords="offset points",
             ha="right", va="center", fontsize=9, color=INK,
             arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.2))
axD.set_xlim(-1.02, 1.02)
axD.set_yscale("log")
axD.yaxis.grid(True, which="major", zorder=0)
axD.set_axisbelow(True)
for s in ("top", "right"):
    axD.spines[s].set_visible(False)
axD.set_xlabel("gain sensitivity  $s$  of a single row", fontsize=9.5)
axD.set_ylabel("% of that population  (log)", fontsize=9.5)
axD.set_title("D  Hits are shifted into the leaking half", fontsize=11, color=INK,
              loc="left", pad=8)
axD.legend(loc="upper left", fontsize=9, frameon=False, labelcolor=INK_2,
           handlelength=1.1, borderaxespad=0.3,
           handles=[Patch(facecolor=MUTED, alpha=0.55, label=f"all {len(flat):,} rows"),
                    Patch(facecolor=BLUE, label=f"the {len(hits):,} hits "
                                                f"(fold $\\geq$ {HIT_FOLD:.0f})")])
axD.text(0.5, -0.30, f"per-batch mean $s$ among hits runs {bat['s_hits'].min():.2f} to "
                     f"{bat['s_hits'].max():.2f} — no batch is exempt",
         transform=axD.transAxes, ha="center", fontsize=8.5, color=MUTED)

fig.suptitle("Yes — with the background fixed at 3000, a hotter batch reports higher "
             "enrichment", fontsize=14.5, color=INK, x=0.072, ha="left", y=0.972)
fig.text(0.072, 0.930,
         "The within-batch fold cancels a batch gain only for the part of the ratio that "
         "was really measured. The 3000 placeholder does not scale, so the gain",
         fontsize=9.5, color=INK_2, ha="left")
fig.text(0.072, 0.905,
         f"passes through in proportion to $s$ = (floored share of the background) − "
         f"(floored share of the target): {hits['s'].mean():+.2f} on average for our hits, "
         f"~0 for everything else.",
         fontsize=9.5, color=INK_2, ha="left")
fig.text(0.072, 0.880,
         f"{DATA_DIR}  ·  {big['ASMS_BATCH_NAME'].nunique()} batches x {N_PROT} proteins  ·  "
         f"{(R == FLOOR).mean():.0%} of all replicate values are the 3000 placeholder",
         fontsize=9.5, color=MUTED, ha="left")

_png = os.path.join(RESULTS_DIR, "floor_gain_demo.png")
fig.savefig(_png, dpi=200, facecolor=SURFACE)
print(f"\nSaved {_png}\nSaved {_csv}")
