import csv, glob, statistics, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = r"d:/0000-UHN/EASMS-data-processing/EASMS-data-processing/MultiBatch_20to36"
OUT = r"d:/0000-UHN/EASMS-data-processing/EASMS-data-processing/MultiBatch_20to36_plots"
import os; os.makedirs(OUT, exist_ok=True)

cols = ['SPR_CATEGORY','SPR_KD_M','SPR_Binding_Affinity',
        'POS_INT_REP1','POS_INT_REP2','POS_INT_REP3']
data = []  # (med_int, kd, aff, cat)
for f in glob.glob(os.path.join(BASE, "*.csv")):
    r = csv.reader(open(f, encoding="utf-8", errors="replace"))
    hdr = next(r)
    I = {c: hdr.index(c) for c in cols if c in hdr}
    if 'SPR_CATEGORY' not in I:
        continue
    for row in r:
        g = lambda c: row[I[c]].strip() if c in I and I[c] < len(row) else ''
        cat = g('SPR_CATEGORY')
        if cat in ('not_found', ''):
            continue
        ints = []
        for c in ['POS_INT_REP1','POS_INT_REP2','POS_INT_REP3']:
            try: ints.append(float(g(c)))
            except: pass
        if not ints:
            continue
        def num(c):
            try: return float(g(c))
            except: return None
        data.append((statistics.median(ints), num('SPR_KD_M'), num('SPR_Binding_Affinity'), cat))

print("rows with SPR + intensity:", len(data))

COLORS = {
    'Confirmed hit':            '#1b7837',
    'Binder':                   '#4393c3',
    'Weak hit / unreliable KD': '#d6604d',
    'Unconfirmed':              '#999999',
}
CATS = ['Confirmed hit','Binder','Weak hit / unreliable KD','Unconfirmed']

def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i]); rk=[0]*len(v)
        for pos,i in enumerate(order): rk[i]=pos
        return rk
    rx, ry = rank(xs), rank(ys)
    mx, my = statistics.mean(rx), statistics.mean(ry)
    sx, sy = statistics.pstdev(rx), statistics.pstdev(ry)
    return sum((a-mx)*(b-my) for a,b in zip(rx,ry))/(len(rx)*sx*sy)

def scatter(idx, xlabel, fname, invert=False):
    fig, ax = plt.subplots(figsize=(7.5, 6))
    xs_all, ys_all = [], []
    for cat in CATS:
        pts = [(d[idx], d[0]) for d in data if d[3]==cat and d[idx] and d[idx]>0 and d[0]>0]
        if not pts: continue
        xs=[p[0] for p in pts]; ys=[p[1] for p in pts]
        xs_all += xs; ys_all += ys
        ax.scatter(xs, ys, s=28, alpha=0.6, label=f"{cat} (n={len(pts)})",
                   color=COLORS[cat], edgecolors='none')
    rho = spearman([math.log10(x) for x in xs_all], [math.log10(y) for y in ys_all])
    ax.set_xscale('log'); ax.set_yscale('log')
    if invert: ax.invert_xaxis()
    ax.set_xlabel(xlabel); ax.set_ylabel("ASMS median intensity (POS_INT_REP1-3)")
    ax.set_title(f"ASMS intensity vs {xlabel}\nSpearman rho = {rho:.3f}  (n={len(xs_all)})")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(True, which='both', ls=':', alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, fname), dpi=150)
    print("saved", fname, "rho=", round(rho,3))

# KD: lower KD = stronger binding, so invert x so "stronger" is to the right
scatter(1, "SPR_KD_M (M)", "asms_vs_KD.png", invert=True)
scatter(2, "SPR_Binding_Affinity", "asms_vs_affinity.png")

# Boxplot of ASMS intensity by category
fig, ax = plt.subplots(figsize=(8, 6))
box_data = [[d[0] for d in data if d[3]==cat and d[0]>0] for cat in CATS]
bp = ax.boxplot(box_data, tick_labels=[f"{c}\n(n={len(b)})" for c,b in zip(CATS,box_data)],
                showfliers=False, patch_artist=True)
for patch, cat in zip(bp['boxes'], CATS):
    patch.set_facecolor(COLORS[cat]); patch.set_alpha(0.6)
ax.set_yscale('log'); ax.set_ylabel("ASMS median intensity")
ax.set_title("ASMS median intensity by SPR category")
plt.xticks(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "asms_intensity_by_category.png"), dpi=150)
print("saved asms_intensity_by_category.png")
print("OUTPUT DIR:", OUT)
