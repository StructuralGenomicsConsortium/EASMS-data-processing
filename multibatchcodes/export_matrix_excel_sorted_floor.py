"""
Excel version of the floor-sorted heatmap (hit_matrix_intensity_sorted_floor.py).

Rows = the 1647 molecules with BINARY_LABEL=1 for >=1 protein.
Cols = all 136 proteins. Cell = median(POS_INT_REP1-3) ASMS intensity.
To match the heatmap exactly: not-measured cells are filled with the 3000 floor,
and rows are sorted by number of 3000(floor) cells DESCENDING (most floors on top).

A helper column 'n_floor_3000' (count of 3000 cells in that row) is added right
after SMILES so the sort is visible/auditable in the sheet.

Sheets:
  1) ASMS_median_intensity  - the sorted intensity matrix (missing -> 3000)
  2) SPR_KD_M               - median SPR_KD_M where SPR-tested (blank otherwise)
  3) SPR_Binding_Affinity   - median SPR_Binding_Affinity where SPR-tested
"""
import csv, glob, os, statistics
from collections import defaultdict
import numpy as np
import pandas as pd

BASE = r"d:/0000-UHN/EASMS-data-processing/EASMS-data-processing/MultiBatch_20to36"
OUT  = r"d:/0000-UHN/EASMS-data-processing/EASMS-data-processing/MultiBatch_20to36_plots"
os.makedirs(OUT, exist_ok=True)
XLSX = os.path.join(OUT, "MultiBatch_20to36_matrix_sorted_floor.xlsx")
INT_COLS = ['POS_INT_REP1','POS_INT_REP2','POS_INT_REP3']
FLOOR = 3000.0

files = sorted(glob.glob(os.path.join(BASE, "*.csv")))
proteins = [os.path.basename(f)[:-4] for f in files]

# pass 1: hit compounds + SMILES + per-protein hit count (protein column order)
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

prot_order = sorted(proteins, key=lambda p: -prot_hits[p])

# build the intensity matrix with missing -> FLOOR, count floors per row, sort desc
compounds = sorted(hit_compounds)
M = np.array([[intensity.get(c, {}).get(p, np.nan) for p in prot_order] for c in compounds],
             dtype=float)
n_missing = int((~np.isfinite(M)).sum())
M[~np.isfinite(M)] = FLOOR
floor_count = (M <= FLOOR).sum(axis=1)
order = np.argsort(-floor_count, kind='stable')     # most floors first
compounds = [compounds[i] for i in order]
M = M[order]; floor_count = floor_count[order]
print(f"missing->floor(3000)={n_missing}; row floor-counts top={floor_count[0]} "
      f".. bottom={floor_count[-1]} (out of {len(prot_order)})")

# intensity sheet (sorted, missing filled with floor to match heatmap)
df_int = pd.DataFrame(M, columns=prot_order)
df_int.insert(0, "n_floor_3000", floor_count)
df_int.insert(0, "SMILES", [smiles.get(c, "") for c in compounds])
df_int.insert(0, "COMPOUND_ID", compounds)

# SPR sheets: same row order, blanks where not SPR-tested
def build_spr(store):
    data = {p: [store.get(c, {}).get(p) for c in compounds] for p in prot_order}
    df = pd.DataFrame(data)
    df.insert(0, "SMILES", [smiles.get(c, "") for c in compounds])
    df.insert(0, "COMPOUND_ID", compounds)
    return df
df_kd  = build_spr(kd)
df_aff = build_spr(aff)

with pd.ExcelWriter(XLSX, engine="openpyxl") as xw:
    df_int.to_excel(xw, sheet_name="ASMS_median_intensity", index=False)
    df_kd.to_excel(xw,  sheet_name="SPR_KD_M",              index=False)
    df_aff.to_excel(xw, sheet_name="SPR_Binding_Affinity",  index=False)

print(f"saved {XLSX}")
print(f"  ASMS sheet: {df_int.shape[0]} molecules x {len(prot_order)} proteins "
      f"(sorted by n_floor_3000 desc)")
print(f"  SPR_KD_M non-blank cells: {sum(len(v) for v in kd.values())}")
print(f"  SPR_Affinity non-blank cells: {sum(len(v) for v in aff.values())}")
