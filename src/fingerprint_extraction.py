# -*- coding: utf-8 -*-
"""Molecular fingerprint + property extraction (Step 7).

Computes `MW`, `ALOGP`, and the nine HitGen fingerprints per row from the
`SMILES` column. Optimized to avoid per-row / per-fingerprint redundancy:

  * Parse each SMILES to an RDKit Mol **once**, then hand the Mol to every
    fingerprint generator (a naive version re-parses the SMILES ~10 times per
    row — once inside each of the 9 generators, plus once for MW/ALOGP).
  * Call each generator **once over the whole column** instead of once per row,
    removing per-call numpy wrapping overhead.
  * Compute on the set of **unique** SMILES and map results back to rows, so
    duplicate SMILES (negatives, repeated compounds, isomer expansion) are
    computed only once.
"""

import numpy as np
from rdkit.Chem import Descriptors
import pandas as pd
from fingerprints import (
    HitGenMACCS, HitGenECFP4, HitGenECFP6, HitGenFCFP4, HitGenFCFP6,
    HitGenRDK, HitGenAvalon, HitGenTopTor, HitGenAtomPair,
)
from utils import to_mol

# Order defines the output column order for the fingerprint columns.
FINGERPRINT_CLASSES = {
    "ECFP4": HitGenECFP4(),
    "ECFP6": HitGenECFP6(),
    "FCFP4": HitGenFCFP4(),
    "FCFP6": HitGenFCFP6(),
    "MACCS": HitGenMACCS(),
    "RDK": HitGenRDK(),
    "AVALON": HitGenAvalon(),
    "TOPTOR": HitGenTopTor(),
    "ATOMPAIR": HitGenAtomPair(),
}


def _row_to_string(row):
    """Format one fingerprint row as a comma-separated string.

    Legacy "string" format: integer counts/bits as plain ints (`'1'`, not
    `'1.0'`), and `'nan'` for failed molecules.
    """
    return ",".join("nan" if v != v else str(int(v)) for v in row)


def extract_fingerprints(df, fp_format="array"):
    """Extracts molecular fingerprints and molecular properties for a DataFrame.

    Args:
        df (pd.DataFrame): Input DataFrame containing a "SMILES" column.
        fp_format (str): "array" (default) stores each fingerprint as a numpy
            float32 array; "string" stores it as a comma-separated string
            (legacy format from earlier pipeline versions).

    Returns:
        pd.DataFrame: Input columns + MW, ALOGP, and one column per fingerprint.
    """
    if "SMILES" not in df.columns:
        raise ValueError("Input DataFrame must contain a 'SMILES' column")

    # Unique SMILES → per-row code mapping. use_na_sentinel=False keeps NaN as
    # its own unique entry (rather than code -1).
    codes, uniques = pd.factorize(df["SMILES"], use_na_sentinel=False)

    # Parse each unique SMILES once. to_mol returns a Mol for valid SMILES,
    # passes existing Mols through, and yields None for invalid/NaN input.
    mols = [to_mol(s) for s in uniques]

    # Molecular properties per unique molecule.
    mw_u = np.array([Descriptors.MolWt(m) if m is not None else np.nan for m in mols])
    alogp_u = np.array([Descriptors.MolLogP(m) if m is not None else np.nan for m in mols])

    out = df.copy()
    # Map unique-level results back to every row via the factorize codes
    # (positional assignment — independent of the DataFrame index).
    out["MW"] = mw_u[codes]
    out["ALOGP"] = alogp_u[codes]

    # Each generator runs once over all unique mols → (U, dimension) array.
    # None mols come back as a row of NaN (handled inside generate_fps).
    for name, fp in FINGERPRINT_CLASSES.items():
        fps_u = fp.generate_fps(mols)
        if fp_format == "string":
            per_unique = [_row_to_string(r) for r in fps_u]
        else:
            per_unique = [r.astype(np.float32) for r in fps_u]
        out[name] = [per_unique[c] for c in codes]

    return out
