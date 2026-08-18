"""
Export the molecule x protein matrices to an Excel workbook with 3 sheets:
  1) ASMS_median_intensity  - median(POS_INT_REP1-3) per molecule x protein
  2) SPR_KD_M               - median SPR_KD_M where SPR-tested (blank otherwise)
  3) SPR_Binding_Affinity   - median SPR_Binding_Affinity where SPR-tested

Rows = the 1647 molecules selected (BINARY_LABEL=1) for >=1 protein.
First column = SMILES (COMPOUND_ID kept as a helper column), then one column per protein.
Not-measured ASMS cells are left blank.
"""
import csv, glob, os, statistics
from collections import defaultdict
import pandas as pd

BASE = r"d:/0000-UHN/EASMS-data-processing/EASMS-data-processing/MultiBatch_20to36"
OUT  = r"d:/0000-UHN/EASMS-data-processing/EASMS-data-processing/MultiBatch_20to36_plots"
os.makedirs(OUT, exist_ok=True)
XLSX = os.path.join(OUT, "MultiBatch_20to36_matrix.xlsx")
INT_COLS = ['POS_INT_REP1','POS_INT_REP2','POS_INT_REP3']

files = sorted(glob.glob(os.path.join(BASE, "*.csv")))
proteins = [os.path.basename(f)[:-4] for f in files]

# pass 1: hit compounds + SMILES + per-protein hit count (for row ordering)
hit_compounds = set()
smiles = {}
prot_hits = defaultdict(int)
for f in files:
    prot = os.path.basename(f)[:-4]
    r = csv.reader(open(f, encoding="utf-8", errors="replace")); hdr = next(r)
    I = {c: hdr.index(c) for c in ['COMPOUND_ID','SMILES','BINARY_LABEL'] if c in hdr}
    for row in r:
        comp = row[I['COMPOUND_ID']].strip()
        if 'SMILES' in I and comp not in smiles:
            smiles[comp] = row[I['SMILES']].strip()
        if row[I['BINARY_LABEL']].strip() == '1':
            hit_compounds.add(comp); prot_hits[prot] += 1
print(f"hit compounds: {len(hit_compounds)} | proteins: {len(proteins)}")

# pass 2: collect the three metrics per (compound, protein)
intensity = defaultdict(dict)  # comp -> {prot: value}
kd        = defaultdict(dict)
aff       = defaultdict(dict)
primary   = {}                 # comp -> (protein rank, -intensity) for ordering
prank = {p: i for i, p in enumerate(sorted(proteins, key=lambda p: -prot_hits[p]))}
for f in files:
    prot = os.path.basename(f)[:-4]
    r = csv.reader(open(f, encoding="utf-8", errors="replace")); hdr = next(r)
    I = {c: hdr.index(c) for c in
         ['COMPOUND_ID','BINARY_LABEL','SPR_CATEGORY','SPR_KD_M','SPR_Binding_Affinity']
         + INT_COLS if c in hdr}
    for row in r:
        comp = row[I['COMPOUND_ID']].strip()
        if comp not in hit_compounds:
            continue
        ii = []
        for c in INT_COLS:
            try: ii.append(float(row[I[c]]))
            except: pass
        if ii:
            intensity[comp][prot] = statistics.median(ii)
        if 'SPR_CATEGORY' in I and row[I['SPR_CATEGORY']].strip() not in ('not_found', ''):
            for col, store in [('SPR_KD_M', kd), ('SPR_Binding_Affinity', aff)]:
                if col in I:
                    try: store[comp][prot] = float(row[I[col]])
                    except: pass
        if row[I['BINARY_LABEL']].strip() == '1':
            v = intensity[comp].get(prot, 0)
            key = (prank[prot], -v)
            if comp not in primary or key < primary[comp]:
                primary[comp] = key

# row order: group molecules by their primary (hit) protein; col order: proteins by hit count
compounds = sorted(hit_compounds, key=lambda c: primary.get(c, (9999, 0)))
prot_order = sorted(proteins, key=lambda p: -prot_hits[p])

def build_df(store):
    data = {p: [store.get(c, {}).get(p) for c in compounds] for p in prot_order}
    df = pd.DataFrame(data, index=compounds)
    df.insert(0, "SMILES", [smiles.get(c, "") for c in compounds])
    df.index.name = "COMPOUND_ID"
    return df.reset_index()   # COMPOUND_ID as first real column, SMILES second

df_int = build_df(intensity)
df_kd  = build_df(kd)
df_aff = build_df(aff)

# make SMILES the very first column, then COMPOUND_ID, then proteins
def reorder(df):
    cols = ['SMILES', 'COMPOUND_ID'] + [c for c in df.columns if c not in ('SMILES','COMPOUND_ID')]
    return df[cols]
df_int, df_kd, df_aff = reorder(df_int), reorder(df_kd), reorder(df_aff)

with pd.ExcelWriter(XLSX, engine="openpyxl") as xw:
    df_int.to_excel(xw, sheet_name="ASMS_median_intensity", index=False)
    df_kd.to_excel(xw,  sheet_name="SPR_KD_M",              index=False)
    df_aff.to_excel(xw, sheet_name="SPR_Binding_Affinity",  index=False)

print(f"saved {XLSX}")
print(f"  ASMS sheet: {df_int.shape[0]} molecules x {len(prot_order)} proteins")
print(f"  SPR_KD_M non-blank cells: {sum(len(v) for v in kd.values())}")
print(f"  SPR_Affinity non-blank cells: {sum(len(v) for v in aff.values())}")
