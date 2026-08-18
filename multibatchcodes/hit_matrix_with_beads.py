"""
Same ASMS median-intensity matrix as hit_matrix_intensity.py, PLUS a 'Beads' column
appended on the right (mean of the 6 bead replicate columns from beads_clean.xlsx),
drawn wider and separated by a gap so it stands out. New figure - does not touch the
existing hit_matrix_intensity.png.
"""
import csv, glob, os, statistics
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

BASE  = r"d:/0000-UHN/EASMS-data-processing/EASMS-data-processing/MultiBatch_20to36"
OUT   = r"d:/0000-UHN/EASMS-data-processing/EASMS-data-processing/MultiBatch_20to36_plots"
BEADS = r"d:/0000-UHN/EASMS-data-processing/EASMS-data-processing/MultiBatchResults/beads_clean.xlsx"
os.makedirs(OUT, exist_ok=True)
INT_COLS  = ['POS_INT_REP1','POS_INT_REP2','POS_INT_REP3']
BEAD_COLS = ['5Beads_V1P0113_1','5Beads_V1P0113_2','5Beads_V1P0113_3',
             '10Beads_V1P0113_1','10Beads_V1P0113_2','10Beads_V1P0113_3']
FLOOR = 3000.0

def med_int(row, I):
    ii = []
    for c in INT_COLS:
        try: ii.append(float(row[I[c]]))
        except: pass
    return statistics.median(ii) if ii else np.nan

files = sorted(glob.glob(os.path.join(BASE, "*.csv")))

# ---- pass 1: hit compounds + per-protein hit counts ----
hit_compounds = set(); prot_hits = {}
for f in files:
    prot = os.path.basename(f)[:-4]
    r = csv.reader(open(f, encoding="utf-8", errors="replace")); hdr = next(r)
    I = {c: hdr.index(c) for c in ['COMPOUND_ID','BINARY_LABEL']}
    prot_hits[prot] = 0
    for row in r:
        if row[I['BINARY_LABEL']].strip() == '1':
            hit_compounds.add(row[I['COMPOUND_ID']].strip()); prot_hits[prot] += 1

# ---- pass 2: intensities + flags ----
val = {}; is_hit = {}; has_spr = {}
for f in files:
    prot = os.path.basename(f)[:-4]
    r = csv.reader(open(f, encoding="utf-8", errors="replace")); hdr = next(r)
    I = {c: hdr.index(c) for c in ['COMPOUND_ID','BINARY_LABEL','SPR_CATEGORY']+INT_COLS if c in hdr}
    for row in r:
        comp = row[I['COMPOUND_ID']].strip()
        if comp not in hit_compounds: continue
        val[(comp, prot)] = med_int(row, I)
        is_hit[(comp, prot)] = row[I['BINARY_LABEL']].strip() == '1'
        if 'SPR_CATEGORY' in I:
            has_spr[(comp, prot)] = row[I['SPR_CATEGORY']].strip() not in ('not_found','')

# ---- bead mean (mean of the 6 bead replicates) per compound ----
b = pd.read_excel(BEADS)
b["bead_mean"] = b[BEAD_COLS].mean(axis=1)
bead_map = (b.rename(columns={"SGC ID for Component": "COMPOUND_ID"})
              .dropna(subset=["COMPOUND_ID"])
              .groupby("COMPOUND_ID")["bead_mean"].mean().to_dict())
n_bead_matched = sum(1 for c in hit_compounds if c in bead_map)
print(f"hit compounds: {len(hit_compounds)} | matched to bead file: {n_bead_matched}")

# ---- ordering (same as intensity plot) ----
proteins = sorted(prot_hits, key=lambda p: -prot_hits[p])
pidx = {p: j for j, p in enumerate(proteins)}
comp_primary = {}
for (comp, prot), h in is_hit.items():
    if h:
        j = pidx[prot]
        if comp not in comp_primary or j < comp_primary[comp][0]:
            comp_primary[comp] = (j, -(val.get((comp,prot)) or 0))
compounds = sorted(hit_compounds, key=lambda c: comp_primary.get(c, (999, 0)))
cidx = {c: i for i, c in enumerate(compounds)}

nR, nC = len(compounds), len(proteins)
M = np.full((nR, nC), np.nan)
for (comp, prot), v in val.items():
    M[cidx[comp], pidx[prot]] = v
M[~np.isfinite(M)] = FLOOR

# bead column vector (missing -> floor, same as matrix convention)
bead_col = np.array([[bead_map.get(c, FLOOR)] for c in compounds])
bead_col[~np.isfinite(bead_col)] = FLOOR

vmin = FLOOR
vmax = np.nanpercentile(np.concatenate([M.ravel(), bead_col.ravel()]), 99.9)
norm = LogNorm(vmin=vmin, vmax=vmax)
cmap = plt.get_cmap('viridis').copy()

# ---- plot: main heatmap + gap + wide bead column ----
fig_h = max(8, nR * 0.026)
fig = plt.figure(figsize=(17, fig_h))
# width ratios: main matrix (136) | gap | beads (wider per-cell) | colorbar
gs = fig.add_gridspec(1, 4, width_ratios=[nC, 6, 8, 2], wspace=0.06)
ax   = fig.add_subplot(gs[0, 0])
axb  = fig.add_subplot(gs[0, 2], sharey=ax)
axcb = fig.add_subplot(gs[0, 3])

im = ax.imshow(M, aspect='auto', cmap=cmap, norm=norm, interpolation='nearest')
for (comp, prot), h in is_hit.items():
    if h:
        ax.add_patch(Rectangle((pidx[prot]-0.5, cidx[comp]-0.5), 1, 1, fill=False,
                               edgecolor='red', linewidth=0.7))
sx = [pidx[p] for (c,p),s in has_spr.items() if s]
sy = [cidx[c] for (c,p),s in has_spr.items() if s]
ax.scatter(sx, sy, s=3, c='white', marker='.', linewidths=0)
ax.set_xticks(range(nC)); ax.set_xticklabels(proteins, rotation=90, fontsize=6)
ax.set_yticks([]); ax.set_ylabel(f"{nR} molecules selected for >=1 protein (grouped by hit protein)")
ax.set_xlabel("Protein target (all 136)")

# bead column (wide, gapped)
axb.imshow(bead_col, aspect='auto', cmap=cmap, norm=norm, interpolation='nearest')
axb.set_xticks([0]); axb.set_xticklabels(['Beads\n(mean of 6)'], fontsize=9, fontweight='bold')
axb.tick_params(axis='y', labelleft=False)
axb.set_ylabel("")

fig.colorbar(im, cax=axcb).set_label("ASMS median intensity (log)")
fig.suptitle("ASMS median-intensity matrix + bead-only signal column\n"
             "viridis: dark=3000 floor -> yellow=high | red outline = ASMS hit | "
             "white dot = SPR-tested | rightmost = mean of 6 bead reps (bead-only binders)",
             fontsize=11)
p = os.path.join(OUT, "hit_matrix_with_beads.png")
fig.savefig(p, dpi=150, bbox_inches='tight')
print("saved", p)
