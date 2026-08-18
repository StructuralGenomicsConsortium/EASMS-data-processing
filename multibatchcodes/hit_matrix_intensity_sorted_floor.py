"""
Same full-intensity hit matrix as hit_matrix_intensity.py (rows = molecules with
BINARY_LABEL=1 for >=1 protein, cols = all 136 proteins, every cell colored by ASMS
median intensity), BUT rows are sorted by how many cells sit at the 3000 floor:
molecules with the MOST 3000-values (least signal across proteins) go to the TOP,
molecules with the fewest 3000s (most signal) go to the bottom.

A cell counts as "3000" if its ASMS median intensity is at the 3000 floor OR the
molecule was not measured in that protein file (missing -> treated as floor).
"""
import csv, glob, os, statistics
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

BASE = r"d:/0000-UHN/EASMS-data-processing/EASMS-data-processing/MultiBatch_20to36"
OUT  = r"d:/0000-UHN/EASMS-data-processing/EASMS-data-processing/MultiBatch_20to36_plots"
os.makedirs(OUT, exist_ok=True)
INT_COLS = ['POS_INT_REP1','POS_INT_REP2','POS_INT_REP3']
FLOOR = 3000.0

def med_int(row, I):
    ints = []
    for c in INT_COLS:
        try: ints.append(float(row[I[c]]))
        except: pass
    return statistics.median(ints) if ints else np.nan

files = sorted(glob.glob(os.path.join(BASE, "*.csv")))

# ---- pass 1: find hit compounds and per-protein hit counts ----
hit_compounds = set()
prot_hits = {}
for f in files:
    prot = os.path.basename(f)[:-4]
    r = csv.reader(open(f, encoding="utf-8", errors="replace"))
    hdr = next(r)
    I = {c: hdr.index(c) for c in ['COMPOUND_ID','BINARY_LABEL'] if c in hdr}
    prot_hits[prot] = 0
    for row in r:
        if row[I['BINARY_LABEL']].strip() == '1':
            hit_compounds.add(row[I['COMPOUND_ID']].strip())
            prot_hits[prot] += 1
print(f"hit compounds: {len(hit_compounds)} | proteins: {len(files)}")

# ---- pass 2: median intensity of every hit-compound in every protein + flags ----
proteins = [os.path.basename(f)[:-4] for f in files]
val = {}    # (comp, prot) -> med intensity
is_hit = {} # (comp, prot) -> bool
has_spr = {}
for f in files:
    prot = os.path.basename(f)[:-4]
    r = csv.reader(open(f, encoding="utf-8", errors="replace"))
    hdr = next(r)
    I = {c: hdr.index(c) for c in
         ['COMPOUND_ID','BINARY_LABEL','SPR_CATEGORY'] + INT_COLS if c in hdr}
    for row in r:
        comp = row[I['COMPOUND_ID']].strip()
        if comp not in hit_compounds:
            continue
        val[(comp, prot)] = med_int(row, I)
        is_hit[(comp, prot)] = row[I['BINARY_LABEL']].strip() == '1'
        if 'SPR_CATEGORY' in I:
            has_spr[(comp, prot)] = row[I['SPR_CATEGORY']].strip() not in ('not_found','')

# ---- proteins ordered by hit count desc (same as before) ----
proteins = sorted(proteins, key=lambda p: -prot_hits[p])
pidx = {p: j for j, p in enumerate(proteins)}

# ---- build matrix keyed by a stable compound list first ----
compounds = sorted(hit_compounds)
cidx = {c: i for i, c in enumerate(compounds)}
nR, nC = len(compounds), len(proteins)
M = np.full((nR, nC), np.nan)
for (comp, prot), v in val.items():
    M[cidx[comp], pidx[prot]] = v

# molecules not present in a protein file -> treat as floor (not detected)
n_missing = int((~np.isfinite(M)).sum())
M[~np.isfinite(M)] = FLOOR

# ---- SORT rows by number of 3000(floor) cells, descending (most floors on top) ----
floor_count = (M <= FLOOR).sum(axis=1)          # per-row count of 3000/floor cells
order = np.argsort(-floor_count, kind='stable')  # most floors first
M = M[order]
compounds = [compounds[i] for i in order]
cidx = {c: i for i, c in enumerate(compounds)}
print(f"row floor-counts: top={floor_count[order][0]} .. bottom={floor_count[order][-1]} "
      f"(out of {nC} proteins)")

vmin = np.nanmin(M)
vmax_full = np.nanmax(M)
vmax = np.nanpercentile(M, 99.9)
print(f"matrix {nR} x {nC}; missing->floor(3000)={n_missing}; "
      f"intensity {vmin:.0f}..{vmax_full:.0f}; color-cap vmax={vmax:.0f}")

# ---- plot ----
fig_h = max(8, nR * 0.026)
fig, ax = plt.subplots(figsize=(16, fig_h))
cmap = plt.get_cmap('viridis').copy()
im = ax.imshow(M, aspect='auto', cmap=cmap,
               norm=LogNorm(vmin=max(vmin, 1), vmax=vmax),
               interpolation='nearest')

# outline ASMS hits
for (comp, prot), h in is_hit.items():
    if h:
        x, y = pidx[prot], cidx[comp]
        ax.add_patch(Rectangle((x-0.5, y-0.5), 1, 1, fill=False,
                               edgecolor='red', linewidth=0.7))
# small dot on SPR-tested cells
sx, sy = [], []
for (comp, prot), s in has_spr.items():
    if s:
        sx.append(pidx[prot]); sy.append(cidx[comp])
ax.scatter(sx, sy, s=3, c='white', marker='.', linewidths=0)

ax.set_xticks(range(nC)); ax.set_xticklabels(proteins, rotation=90, fontsize=6)
ax.set_yticks([]); ax.set_ylabel(f"{nR} molecules - sorted by #3000(floor) cells (most floors on top)")
ax.set_xlabel("Protein target (all 136)")
ax.set_title("ASMS median intensity matrix - rows sorted by number of 3000(floor) values\n"
             "top = most 3000s (least signal) -> bottom = fewest 3000s (most signal) | "
             "viridis dark=3000 floor -> yellow=high | red outline = ASMS hit | white dot = SPR-tested")
cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
cb.set_label("ASMS median intensity (log)")
fig.tight_layout()
p = os.path.join(OUT, "hit_matrix_intensity_sorted_floor.png")
fig.savefig(p, dpi=150, bbox_inches='tight')
print("saved", p)
