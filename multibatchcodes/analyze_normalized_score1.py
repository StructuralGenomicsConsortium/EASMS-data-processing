# -*- coding: utf-8 -*-
"""Characterise "Normalized score1": its range, its relation to fold_change, and to
BINARY_LABEL -- overall and PER PROTEIN.

Reads a *_modified folder produced by run_multibatch_pipeline.py.

    python multibatchcodes/analyze_normalized_score1.py --dir MultiBatch_20to22_modified

Reports
  1. range / distribution of the score (and how much of it is undefined, and why)
  2. vs the existing fold_change: Spearman + Pearson on log2 fold, overall and per protein
  3. vs BINARY_LABEL: AUC, medians per class, and the same per protein
  4. a comparison of the score's discrimination against fold_change's, per protein
Writes one summary CSV per section plus a 4-panel figure.
"""

import os
import argparse
import glob
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr, mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCORE = "Normalized score1"
SURFACE, INK, INK_2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, BASELINE = "#e1e0d9", "#c3c2b7"
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a"]
plt.rcParams.update({
    "font.family": ["Segoe UI", "DejaVu Sans", "sans-serif"],
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "axes.edgecolor": BASELINE, "axes.linewidth": 0.8,
    "axes.labelcolor": INK_2, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK_2, "ytick.labelcolor": INK_2,
    "grid.color": GRID, "grid.linewidth": 0.8,
})


def auc(pos, neg):
    """Mann-Whitney AUC = P(score of a positive > score of a negative)."""
    if len(pos) < 3 or len(neg) < 3:
        return np.nan
    u = mannwhitneyu(pos, neg, alternative="two-sided").statistic
    return u / (len(pos) * len(neg))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", default="MultiBatchResults")
    a = ap.parse_args()
    tag = os.path.basename(a.dir.rstrip("/\\")).replace("_modified", "")
    os.makedirs(a.out, exist_ok=True)

    cols = ["COMPOUND_ID", "ASMS_BATCH_NAME", "BINARY_LABEL", "fold_change",
            "n_detected", "target_median", SCORE, "ns1_enrichment_within",
            "ns1_log2_enrichment", "ns1_scale_factor",
            "SPR_LABEL", "pass_fold", "Pipeline_hits"]
    parts = []
    for p in sorted(glob.glob(os.path.join(a.dir, "*.csv"))):
        head = pd.read_csv(p, nrows=0).columns
        d = pd.read_csv(p, usecols=[c for c in cols if c in head], low_memory=False)
        d["protein"] = os.path.basename(p)[:-4]
        parts.append(d)
    df = pd.concat(parts, ignore_index=True)
    df["BINARY_LABEL"] = pd.to_numeric(df["BINARY_LABEL"], errors="coerce")
    s = df[SCORE]

    # ---------------------------------------------------------------- 1. range
    print("=" * 92)
    print(f"1. RANGE of '{SCORE}'   ({tag}, {len(df):,} rows)")
    print("=" * 92)
    fin = s.dropna()
    print(f"  defined      {len(fin):,} / {len(df):,}  ({len(fin)/len(df):.1%})")
    print(f"  undefined    {s.isna().sum():,}   (MAD=0, or <3 proteins per compound)")
    qs = [0, 0.1, 1, 5, 25, 50, 75, 95, 99, 99.9, 100]
    print("  percentiles: " + "  ".join(
        f"p{q:g}={np.percentile(fin, q):,.2f}" for q in qs))
    print(f"  mean {fin.mean():,.3f}   sd {fin.std():,.3f}   "
          f"min {fin.min():,.2f}   max {fin.max():,.2f}")
    for th in (3, 5, 10, 20):
        print(f"  |score| >= {th:<3}: {int((fin.abs() >= th).sum()):>8,} "
              f"({(fin.abs() >= th).mean():.2%})")

    # -------------------------------------------------- 2. vs existing fold_change
    print()
    print("=" * 92)
    print("2. RELATION TO fold_change")
    print("=" * 92)
    d2 = df.dropna(subset=[SCORE, "fold_change"]).copy()
    d2 = d2[d2["fold_change"] > 0]
    d2["log2_fold"] = np.log2(d2["fold_change"])
    rs = spearmanr(d2[SCORE], d2["fold_change"]).statistic
    rp = pearsonr(d2[SCORE], d2["log2_fold"]).statistic
    print(f"  overall  Spearman(score, fold_change) = {rs:.4f}   "
          f"Pearson(score, log2 fold) = {rp:.4f}   n = {len(d2):,}")
    rows = []
    for prot, g in d2.groupby("protein"):
        rows.append({"protein": prot, "batch": g["ASMS_BATCH_NAME"].iloc[0],
                     "n": len(g),
                     "spearman_vs_fold": spearmanr(g[SCORE], g["fold_change"]).statistic,
                     "pearson_vs_log2fold": pearsonr(g[SCORE], g["log2_fold"]).statistic})
    per_fold = pd.DataFrame(rows).sort_values("spearman_vs_fold")
    print(per_fold.to_string(index=False, float_format=lambda x: f"{x:,.4f}"))

    # ------------------------------------------------------- 3. vs BINARY_LABEL
    print()
    print("=" * 92)
    print("3. RELATION TO BINARY_LABEL")
    print("=" * 92)
    # BINARY_LABEL is not independent of fold_change -- quantify how circular it is
    # before using it to compare the two scores.
    for c in ("pass_fold", "Pipeline_hits"):
        if c in df.columns:
            x = pd.to_numeric(df[c], errors="coerce")
            print(f"  [circularity] BINARY_LABEL agrees with {c:<14}: "
                  f"{(x == df['BINARY_LABEL']).mean():.4%}")
    _p = df.loc[df.BINARY_LABEL == 1, "fold_change"].dropna()
    if len(_p):
        print(f"  [circularity] every BINARY_LABEL=1 row has fold_change >= {_p.min():.3f}"
              f"  (pipeline hit threshold is {5.0})")
        print("  => fold_change essentially DEFINES this label, so its AUC below is"
              "\n     near 1.0 by construction and is NOT a fair benchmark.")
    print()

    d3 = df.dropna(subset=[SCORE, "BINARY_LABEL"])
    pos, neg = d3.loc[d3.BINARY_LABEL == 1, SCORE], d3.loc[d3.BINARY_LABEL == 0, SCORE]
    print(f"  n(label=1) = {len(pos):,}   n(label=0) = {len(neg):,}")
    print(f"  median score  label=1: {pos.median():,.3f}   label=0: {neg.median():,.3f}")
    print(f"  AUC(score)        = {auc(pos, neg):.4f}")
    d3f = d3.dropna(subset=["fold_change"])
    print(f"  AUC(fold_change)  = "
          f"{auc(d3f.loc[d3f.BINARY_LABEL==1,'fold_change'], d3f.loc[d3f.BINARY_LABEL==0,'fold_change']):.4f}"
          "   <- the incumbent, same rows")
    rows = []
    for prot, g in d3.groupby("protein"):
        p1, p0 = g.loc[g.BINARY_LABEL == 1, SCORE], g.loc[g.BINARY_LABEL == 0, SCORE]
        gf = g.dropna(subset=["fold_change"])
        rows.append({
            "protein": prot, "batch": g["ASMS_BATCH_NAME"].iloc[0],
            "n_pos": len(p1), "n_neg": len(p0),
            "med_pos": p1.median(), "med_neg": p0.median(),
            "auc_score": auc(p1, p0),
            "auc_fold": auc(gf.loc[gf.BINARY_LABEL == 1, "fold_change"],
                            gf.loc[gf.BINARY_LABEL == 0, "fold_change"]),
        })
    per_lab = pd.DataFrame(rows)
    per_lab["auc_diff"] = per_lab["auc_score"] - per_lab["auc_fold"]
    per_lab = per_lab.sort_values("auc_score", ascending=False)
    print(per_lab.to_string(index=False, float_format=lambda x: f"{x:,.4f}"))
    ok = per_lab.dropna(subset=["auc_diff"])
    print(f"\n  score beats fold_change on {int((ok.auc_diff > 0).sum())} / {len(ok)} proteins"
          f"   mean AUC diff {ok.auc_diff.mean():+.4f}")

    # ------------------------------------- 4. vs SPR_LABEL (INDEPENDENT truth)
    print()
    print("=" * 92)
    print("4. RELATION TO SPR_LABEL  — the independent ground truth")
    print("=" * 92)
    if "SPR_LABEL" in df.columns and df["SPR_LABEL"].notna().sum() >= 20:
        d4 = df.dropna(subset=["SPR_LABEL"])
        sp, sn = d4.SPR_LABEL == 1, d4.SPR_LABEL == 0
        print(f"  rows with an SPR label: {len(d4):,}  "
              f"({int(sp.sum()):,} positive / {int(sn.sum()):,} negative)")
        cov = d4[SCORE].notna().mean()
        print(f"  of those, {SCORE} is defined for {cov:.1%}")
        both = d4.dropna(subset=[SCORE])
        b1, b0 = both.loc[both.SPR_LABEL == 1, SCORE], both.loc[both.SPR_LABEL == 0, SCORE]
        print(f"  median  SPR=1: {b1.median():,.3f}   SPR=0: {b0.median():,.3f}")
        print(f"  AUC({SCORE:<18}) = {auc(b1, b0):.4f}   n={len(both):,}")
        fo = d4.dropna(subset=["fold_change"])
        print(f"  AUC(fold_change            ) = "
              f"{auc(fo.loc[fo.SPR_LABEL==1,'fold_change'], fo.loc[fo.SPR_LABEL==0,'fold_change']):.4f}"
              f"   n={len(fo):,}")
        # head-to-head on exactly the rows where BOTH are defined
        hh = d4.dropna(subset=[SCORE, "fold_change"])
        if len(hh) >= 20:
            h1, h0 = (hh.SPR_LABEL == 1).to_numpy(), (hh.SPR_LABEL == 0).to_numpy()
            a_s = auc(hh.loc[h1, SCORE], hh.loc[h0, SCORE])
            a_f = auc(hh.loc[h1, "fold_change"], hh.loc[h0, "fold_change"])
            print(f"  head-to-head on the same {len(hh):,} rows:"
                  f"   score {a_s:.4f}   vs fold {a_f:.4f}")
            # paired bootstrap over rows -- is either AUC above chance, and is the
            # difference real? With a few hundred labelled rows it usually is not.
            rng = np.random.default_rng(0)
            sv = hh[SCORE].to_numpy()
            fv = hh["fold_change"].to_numpy()
            ds, df_, dd = [], [], []
            for _ in range(2000):
                k = rng.integers(0, len(hh), len(hh))
                y1, y0 = h1[k], h0[k]
                if y1.sum() < 3 or y0.sum() < 3:
                    continue
                x, y = auc(sv[k][y1], sv[k][y0]), auc(fv[k][y1], fv[k][y0])
                ds.append(x); df_.append(y); dd.append(x - y)
            q = lambda v: (np.percentile(v, 2.5), np.percentile(v, 97.5))
            print(f"    95% CI  score  [{q(ds)[0]:.3f}, {q(ds)[1]:.3f}]"
                  f"   fold [{q(df_)[0]:.3f}, {q(df_)[1]:.3f}]")
            lo, hi = q(dd)
            print(f"    95% CI  difference (score - fold) [{lo:+.3f}, {hi:+.3f}]"
                  f"  -> {'SIGNIFICANT' if lo > 0 or hi < 0 else 'NOT significant'}")
            print(f"    (AUC 0.5 = chance; CI covering 0.5 means no real"
                  " discrimination on this ground truth)")
    else:
        n_spr = int(df["SPR_LABEL"].notna().sum()) if "SPR_LABEL" in df.columns else 0
        print(f"  only {n_spr} rows carry an SPR label in this dataset — too few to"
              " compare. (20to36 has 605; 20to22 has 1.)")

    per_fold.to_csv(os.path.join(a.out, f"{tag}_ns1_vs_fold_per_protein.csv"), index=False)
    per_lab.to_csv(os.path.join(a.out, f"{tag}_ns1_vs_label_per_protein.csv"), index=False)

    # ------------------------------------------------------------------ figure
    batches = sorted(df["ASMS_BATCH_NAME"].dropna().unique())
    COLOR = dict(zip(batches, (SLOTS * 8)[:len(batches)]))
    fig, ax = plt.subplots(2, 2, figsize=(14, 9))

    a0 = ax[0, 0]
    a0.hist(fin, bins=140, color=SLOTS[0], alpha=0.85, zorder=3)
    a0.set_yscale("log")
    a0.axvline(0, color=MUTED, lw=1)
    a0.set_xlabel(SCORE); a0.set_ylabel("count (log)")
    a0.set_title(f"Distribution   ({len(fin):,} defined of {len(df):,})",
                 fontsize=11, loc="left", color=INK)

    a1 = ax[0, 1]
    sm = d2.sample(min(60000, len(d2)), random_state=0)
    a1.scatter(sm["log2_fold"], sm[SCORE], s=3, alpha=0.18,
               c=[COLOR.get(b, MUTED) for b in sm["ASMS_BATCH_NAME"]],
               linewidths=0, zorder=3)
    a1.set_xlabel("log2( fold_change )  [existing]"); a1.set_ylabel(SCORE)
    a1.set_title(f"vs fold_change   Spearman {rs:.3f}", fontsize=11, loc="left", color=INK)

    a2 = ax[1, 0]
    parts_ = [neg.dropna().to_numpy(), pos.dropna().to_numpy()]
    bp = a2.boxplot(parts_, vert=True, widths=0.5, showfliers=False, patch_artist=True)
    for i, b_ in enumerate(bp["boxes"]):
        b_.set(facecolor=SLOTS[i], alpha=0.25, edgecolor=SLOTS[i], linewidth=1.3)
    for m_ in bp["medians"]:
        m_.set(color=INK, linewidth=2)
    a2.set_xticks([1, 2]); a2.set_xticklabels([f"label 0\n(n={len(neg):,})",
                                               f"label 1\n(n={len(pos):,})"])
    a2.set_ylabel(SCORE)
    a2.set_title(f"vs BINARY_LABEL   AUC {auc(pos, neg):.3f}",
                 fontsize=11, loc="left", color=INK)

    a3 = ax[1, 1]
    o = per_lab.dropna(subset=["auc_score", "auc_fold"])
    a3.scatter(o["auc_fold"], o["auc_score"], s=70, zorder=3,
               c=[COLOR.get(b, MUTED) for b in o["batch"]],
               edgecolor=SURFACE, linewidth=1.5)
    lim = [min(o[["auc_fold", "auc_score"]].min()) - .02,
           max(o[["auc_fold", "auc_score"]].max()) + .02]
    a3.plot(lim, lim, color=MUTED, lw=1, zorder=2)
    a3.set_xlim(lim); a3.set_ylim(lim)
    a3.set_xlabel("AUC of fold_change"); a3.set_ylabel(f"AUC of {SCORE}")
    a3.set_title("Per protein: which separates binders better?",
                 fontsize=11, loc="left", color=INK)

    for x in ax.ravel():
        x.grid(True, zorder=0); x.set_axisbelow(True)
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
    fig.suptitle(f"'{SCORE}' — {tag}", fontsize=14, color=INK, x=0.04, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    png = os.path.join(a.out, f"{tag}_normalized_score1.png")
    fig.savefig(png, dpi=180, facecolor=SURFACE)
    print(f"\nSaved {png}")


if __name__ == "__main__":
    main()
