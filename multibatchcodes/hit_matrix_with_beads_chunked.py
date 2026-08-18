"""
Same figure as hit_matrix_with_beads.py, but split into chunks of 200 molecules,
one PNG per chunk (1647 molecules -> 9 plots). Each chunk shows 200 rows x 136
proteins + gapped wide 'Beads' column, with COMPOUND_ID y-labels (readable at 200/plot).
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
OUT   = r"d:/0000-UHN/EASMS-data-processing/EASMS-data-processing/MultiBatch_20to36_plots/beads_chunks_50"
BEADS = r"d:/0000-UHN/EASMS-data-processing/EASMS-data-processing/MultiBatchResults/beads_clean.xlsx"
os.makedirs(OUT, exist_ok=True)
INT_COLS  = ['POS_INT_REP1','POS_INT_REP2','POS_INT_REP3']
BEAD_COLS = ['5Beads_V1P0113_1','5Beads_V1P0113_2','5Beads_V1P0113_3',
             '10Beads_V1P0113_1','10Beads_V1P0113_2','10Beads_V1P0113_3']
FLOOR = 3000.0
TARGET = 50   # ~molecules per plot; chunks are evenly sized so none end up tiny

def med_int(row, I):
    ii = []
    for c in INT_COLS:
        try: ii.append(float(row[I[c]]))
        except: pass
    return statistics.median(ii) if ii else np.nan

files = sorted(glob.glob(os.path.join(BASE, "*.csv")))

# pass 1: hit compounds + per-protein hit counts
hit_compounds = set(); prot_hits = {}
for f in files:
    prot = os.path.basename(f)[:-4]
    r = csv.reader(open(f, encoding="utf-8", errors="replace")); hdr = next(r)
    I = {c: hdr.index(c) for c in ['COMPOUND_ID','BINARY_LABEL']}
    prot_hits[prot] = 0
    for row in r:
        if row[I['BINARY_LABEL']].strip() == '1':
            hit_compounds.add(row[I['COMPOUND_ID']].strip()); prot_hits[prot] += 1

# pass 2: intensities + flags
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

# bead mean per compound
b = pd.read_excel(BEADS)
b["bead_mean"] = b[BEAD_COLS].mean(axis=1)
bead_map = (b.rename(columns={"SGC ID for Component": "COMPOUND_ID"})
              .dropna(subset=["COMPOUND_ID"])
              .groupby("COMPOUND_ID")["bead_mean"].mean().to_dict())

# ordering (same as full plot)
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
bead_col = np.array([[bead_map.get(c, FLOOR)] for c in compounds])
bead_col[~np.isfinite(bead_col)] = FLOOR

vmin = FLOOR
vmax = np.nanpercentile(np.concatenate([M.ravel(), bead_col.ravel()]), 99.9)
norm = LogNorm(vmin=vmin, vmax=vmax)
cmap = plt.get_cmap('viridis').copy()

# even-sized chunks near TARGET so none is tiny (e.g. 1647 -> 33 chunks of 49-50)
n_chunks = max(1, round(nR / TARGET))
base, rem = divmod(nR, n_chunks)          # first `rem` chunks get base+1
bounds = []
start = 0
for k in range(n_chunks):
    size = base + (1 if k < rem else 0)
    bounds.append((start, start + size)); start += size
print(f"{nR} molecules -> {n_chunks} plots of {base}-{base+1} each")

# precompute hit / spr coordinates for fast per-chunk filtering
hit_xy = [(cidx[c], pidx[p]) for (c,p),h in is_hit.items() if h]
spr_xy = [(cidx[c], pidx[p]) for (c,p),s in has_spr.items() if s]

for k, (r0, r1) in enumerate(bounds):
    rows = r1 - r0
    sub = M[r0:r1]; subb = bead_col[r0:r1]
    labels = compounds[r0:r1]

    fig_h = max(6, rows * 0.06)
    fig = plt.figure(figsize=(17, fig_h))
    gs = fig.add_gridspec(1, 4, width_ratios=[nC, 6, 8, 2], wspace=0.06)
    ax   = fig.add_subplot(gs[0, 0])
    axb  = fig.add_subplot(gs[0, 2], sharey=ax)
    axcb = fig.add_subplot(gs[0, 3])

    im = ax.imshow(sub, aspect='auto', cmap=cmap, norm=norm, interpolation='nearest')
    for (yy, xx) in hit_xy:
        if r0 <= yy < r1:
            ax.add_patch(Rectangle((xx-0.5, yy-r0-0.5), 1, 1, fill=False,
                                   edgecolor='red', linewidth=0.8))
    sx = [xx for (yy,xx) in spr_xy if r0 <= yy < r1]
    sy = [yy-r0 for (yy,xx) in spr_xy if r0 <= yy < r1]
    ax.scatter(sx, sy, s=6, c='white', marker='.', linewidths=0)

    ax.set_xticks(range(nC)); ax.set_xticklabels(proteins, rotation=90, fontsize=6)
    ax.set_yticks(range(rows)); ax.set_yticklabels(labels, fontsize=4)
    ax.set_ylabel(f"molecules {r0+1}-{r1} (grouped by hit protein)")
    ax.set_xlabel("Protein target (all 136)")

    axb.imshow(subb, aspect='auto', cmap=cmap, norm=norm, interpolation='nearest')
    axb.set_xticks([0]); axb.set_xticklabels(['Beads\n(mean of 6)'], fontsize=9, fontweight='bold')
    axb.tick_params(axis='y', labelleft=False)

    fig.colorbar(im, cax=axcb).set_label("ASMS median intensity (log)")
    fig.suptitle(f"ASMS intensity + beads  -  molecules {r0+1}-{r1} of {nR}  (part {k+1}/{n_chunks})\n"
                 "red outline = ASMS hit | white dot = SPR-tested | rightmost = mean of 6 bead reps",
                 fontsize=11)
    p = os.path.join(OUT, f"hit_matrix_beads_part{k+1:02d}.png")
    fig.savefig(p, dpi=150, bbox_inches='tight'); plt.close(fig)
    print("saved", p)

print("all chunks in:", OUT)
