import csv, glob, os
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = r"d:/0000-UHN/EASMS-data-processing/EASMS-data-processing/MultiBatch_20to36"
OUT  = r"d:/0000-UHN/EASMS-data-processing/EASMS-data-processing/MultiBatch_20to36_plots"
os.makedirs(OUT, exist_ok=True)

asms1 = defaultdict(int)     # protein -> # ASMS positive rows
spr   = defaultdict(int)     # protein -> # of those with SPR value
proteins_all = []
for f in glob.glob(os.path.join(BASE, "*.csv")):
    prot = os.path.basename(f)[:-4]
    proteins_all.append(prot)
    r = csv.reader(open(f, encoding="utf-8", errors="replace"))
    hdr = next(r)
    I = {c: hdr.index(c) for c in ['BINARY_LABEL','SPR_CATEGORY'] if c in hdr}
    asms1.setdefault(prot, 0); spr.setdefault(prot, 0)
    for row in r:
        g = lambda c: row[I[c]].strip() if c in I and I[c] < len(row) else ''
        if g('BINARY_LABEL') == '1':
            asms1[prot] += 1
            if g('SPR_CATEGORY') not in ('not_found', ''):
                spr[prot] += 1

proteins_all = sorted(set(proteins_all))
counts = {p: asms1[p] for p in proteins_all}
n_total = len(proteins_all)
n_empty = sum(1 for p in proteins_all if counts[p] == 0)
n_hit   = n_total - n_empty
print(f"{n_total} proteins | {n_hit} with >=1 ASMS positive | {n_empty} empty")

# ================= PLOT 1: ranked horizontal bars, all proteins with hits =================
hit_prots = sorted([p for p in proteins_all if counts[p] > 0], key=lambda p: counts[p])
y = np.arange(len(hit_prots))
spr_part = np.array([spr[p] for p in hit_prots])
rest     = np.array([counts[p] - spr[p] for p in hit_prots])

fig_h = max(8, len(hit_prots) * 0.16)
fig, ax = plt.subplots(figsize=(11, fig_h))
ax.barh(y, spr_part, color='#c0392b', label='ASMS+ with SPR data')
ax.barh(y, rest, left=spr_part, color='#2c7fb8', label='ASMS+ no SPR')
for i, p in enumerate(hit_prots):
    ax.text(counts[p] + max(counts.values())*0.005, i, str(counts[p]),
            va='center', fontsize=6)
ax.set_yticks(y); ax.set_yticklabels(hit_prots, fontsize=6)
ax.set_xlabel("# ASMS-positive molecules (BINARY_LABEL = 1)")
ax.set_title(f"ASMS positives per protein  —  {n_hit} of {n_total} proteins have hits\n"
             f"({n_empty} proteins have 0 ASMS positives, listed separately)")
ax.legend(loc='lower right', frameon=False, fontsize=9)
ax.margins(y=0.005)
fig.tight_layout()
p1 = os.path.join(OUT, "asms_positives_per_protein.png")
fig.savefig(p1, dpi=150, bbox_inches='tight'); print("saved", p1)

# write the empty-protein list to a small text panel image
empties = sorted([p for p in proteins_all if counts[p] == 0])
fig2, ax2 = plt.subplots(figsize=(11, 3.2)); ax2.axis('off')
ax2.set_title(f"{n_empty} proteins with ZERO ASMS positives", fontsize=12, loc='left')
cols = 5
rows = -(-len(empties)//cols)
for k, name in enumerate(empties):
    c, rr = divmod(k, rows)
    ax2.text(c/cols, 1 - (rr+1)/(rows+1), "• "+name, fontsize=7,
             transform=ax2.transAxes, va='top')
fig2.tight_layout()
p2 = os.path.join(OUT, "empty_proteins.png")
fig2.savefig(p2, dpi=150, bbox_inches='tight'); print("saved", p2)

# ================= PLOT 2: distribution — how many proteins have how many positives =======
vals = [counts[p] for p in proteins_all]
bins = [0,1,2,3,6,11,21,51,101, max(vals)+1]
labels = ['0','1','2','3-5','6-10','11-20','21-50','51-100','100+']
hist = [0]*len(labels)
for v in vals:
    for bi in range(len(bins)-1):
        if bins[bi] <= v < bins[bi+1]:
            hist[bi]+=1; break
fig3, ax3 = plt.subplots(figsize=(9, 5))
colors = ['#bbbbbb'] + ['#2c7fb8']*(len(labels)-1)
bars = ax3.bar(labels, hist, color=colors, edgecolor='white')
for b, h in zip(bars, hist):
    if h: ax3.text(b.get_x()+b.get_width()/2, h+0.5, str(h), ha='center', fontsize=9)
ax3.set_xlabel("# ASMS-positive molecules per protein")
ax3.set_ylabel("# proteins")
ax3.set_title(f"Distribution of ASMS positives across {n_total} proteins\n"
              f"(grey = {n_empty} empty proteins)")
ax3.spines[['top','right']].set_visible(False)
fig3.tight_layout()
p3 = os.path.join(OUT, "asms_positive_distribution.png")
fig3.savefig(p3, dpi=150, bbox_inches='tight'); print("saved", p3)
