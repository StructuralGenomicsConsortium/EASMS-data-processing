import csv, glob, os, statistics, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

BASE = r"d:/0000-UHN/EASMS-data-processing/EASMS-data-processing/MultiBatch_20to36"
OUT  = r"d:/0000-UHN/EASMS-data-processing/EASMS-data-processing/MultiBatch_20to36_plots"
os.makedirs(OUT, exist_ok=True)

# ---- collect ASMS=1 cells: (compound, protein) -> {intensity, has_spr} ----
cell = {}                      # (comp, prot) -> dict
prot_hits = {}                 # protein -> hit count
for f in glob.glob(os.path.join(BASE, "*.csv")):
    prot = os.path.basename(f)[:-4]
    r = csv.reader(open(f, encoding="utf-8", errors="replace"))
    hdr = next(r)
    I = {c: hdr.index(c) for c in
         ['COMPOUND_ID','BINARY_LABEL','SPR_CATEGORY',
          'POS_INT_REP1','POS_INT_REP2','POS_INT_REP3'] if c in hdr}
    for row in r:
        g = lambda c: row[I[c]].strip() if c in I and I[c] < len(row) else ''
        if g('BINARY_LABEL') != '1':
            continue
        comp = g('COMPOUND_ID')
        ints = []
        for c in ['POS_INT_REP1','POS_INT_REP2','POS_INT_REP3']:
            try: ints.append(float(g(c)))
            except: pass
        med = statistics.median(ints) if ints else None
        has_spr = g('SPR_CATEGORY') not in ('not_found', '')
        cell[(comp, prot)] = {'int': med, 'spr': has_spr}
        prot_hits[prot] = prot_hits.get(prot, 0) + 1

# floor value = min observed intensity (for hits lacking intensity)
all_int = [v['int'] for v in cell.values() if v['int'] and v['int'] > 0]
floor = min(all_int)
vmax  = max(all_int)
print(f"cells (ASMS=1): {len(cell)}  |  intensity floor={floor:.3g}  vmax={vmax:.3g}")

# ---- axis ordering: proteins by hit count desc; compounds grouped by primary protein ----
proteins = sorted(prot_hits, key=lambda p: -prot_hits[p])
pidx = {p: i for i, p in enumerate(proteins)}

comp_cells = {}
for (comp, prot), v in cell.items():
    comp_cells.setdefault(comp, []).append((prot, v['int'] or floor))
def primary(comp):
    # protein with highest intensity for this compound
    return min(pidx[p] for p, _ in comp_cells[comp]), \
           -max(i for _, i in comp_cells[comp])
compounds = sorted(comp_cells, key=primary)
cidx = {c: i for i, c in enumerate(compounds)}

nR, nC = len(compounds), len(proteins)
M = np.full((nR, nC), np.nan)
spr_mask = np.zeros((nR, nC), dtype=bool)
for (comp, prot), v in cell.items():
    i, j = cidx[comp], pidx[prot]
    M[i, j] = v['int'] if (v['int'] and v['int'] > 0) else floor
    spr_mask[i, j] = v['spr']

print(f"matrix: {nR} molecules x {nC} proteins ; SPR-tested cells outlined: {spr_mask.sum()}")

# ---- plot ----
fig_h = max(8, nR * 0.026)          # ~0.026 in per row
fig, ax = plt.subplots(figsize=(15, fig_h))
cmap = plt.get_cmap('viridis').copy()
im = ax.imshow(M, aspect='auto', cmap=cmap,
               norm=LogNorm(vmin=floor, vmax=vmax),
               interpolation='nearest')

# outline SPR-tested cells
ys, xs = np.where(spr_mask)
for y, x in zip(ys, xs):
    ax.add_patch(Rectangle((x-0.5, y-0.5), 1, 1, fill=False,
                           edgecolor='red', linewidth=0.6))

ax.set_xticks(range(nC))
ax.set_xticklabels(proteins, rotation=90, fontsize=6)
ax.set_yticks([]); ax.set_ylabel(f"{nR} selected molecules (ASMS=1), grouped by protein")
ax.set_xlabel("Protein target")
ax.set_title("ASMS hit matrix - molecules x proteins\n"
             "color = ASMS median intensity (log; floor for hits w/o intensity) | "
             "red outline = SPR-tested cell")
cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
cb.set_label("ASMS median intensity")
fig.tight_layout()
p = os.path.join(OUT, "hit_matrix_overview.png")
fig.savefig(p, dpi=150, bbox_inches='tight')
print("saved", p)
