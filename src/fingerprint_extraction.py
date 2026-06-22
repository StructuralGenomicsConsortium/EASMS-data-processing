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
  * Fan the whole per-molecule cost (parse + MW/ALOGP + all 9 fingerprints) out
    across a single process pool, so parsing and featurization run in parallel
    on every core (RDKit holds the GIL, so processes, not threads). One pool for
    the whole step — not one per fingerprint.
"""

from concurrent.futures import ProcessPoolExecutor
from functools import partial

import numpy as np
from rdkit.Chem import Descriptors
import pandas as pd
from fingerprints import (
    HitGenMACCS, HitGenECFP4, HitGenECFP6, HitGenFCFP4, HitGenFCFP6,
    HitGenRDK, HitGenAvalon, HitGenTopTor, HitGenAtomPair,
    _chunk, _default_n_jobs, _PARALLEL_MIN_ITEMS,
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


def _featurize_chunk(smiles_chunk, fp_format):
    """Worker: featurize one chunk of UNIQUE SMILES.

    Parses each SMILES in the chunk once, then computes MW, ALOGP and every
    fingerprint from that single Mol. Generators run serially here (n_jobs=1) so
    there is no nested parallelism inside a worker. Defined at module level so
    ProcessPoolExecutor can pickle it.

    Returns (mw_array, alogp_array, {fp_name: [per-molecule value, ...]}) for the
    chunk, in order. Per-molecule values are float32 arrays ("array" format) or
    comma-separated strings ("string" format).
    """
    mols = [to_mol(s) for s in smiles_chunk]
    mw = np.array([Descriptors.MolWt(m) if m is not None else np.nan for m in mols])
    alogp = np.array([Descriptors.MolLogP(m) if m is not None else np.nan for m in mols])

    fps = {}
    for name, fp in FINGERPRINT_CLASSES.items():
        arr = fp.generate_fps(mols, n_jobs=1)   # serial within the worker
        if fp_format == "string":
            fps[name] = [_row_to_string(r) for r in arr]
        else:
            fps[name] = [r.astype(np.float32) for r in arr]
    return mw, alogp, fps


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
    uniques = list(uniques)

    # Featurize each UNIQUE molecule once, fanned out across one process pool:
    # each worker parses its chunk's SMILES and computes MW/ALOGP + all 9
    # fingerprints. Below the threshold (or on any failure) we run serially.
    n_jobs = _default_n_jobs()
    n_jobs = min(n_jobs, len(uniques)) if uniques else 1

    results = None
    if n_jobs > 1 and len(uniques) >= _PARALLEL_MIN_ITEMS:
        try:
            worker = partial(_featurize_chunk, fp_format=fp_format)
            with ProcessPoolExecutor(max_workers=n_jobs) as ex:
                results = list(ex.map(worker, _chunk(uniques, n_jobs)))
        except Exception as e:
            # Never let parallelism break the pipeline -- fall back to serial.
            print(f"  [fingerprints] parallel featurization failed "
                  f"({type(e).__name__}: {e}); falling back to serial.")
            results = None

    if results is not None:
        # Reassemble per-unique results in chunk order.
        mw_u = np.concatenate([r[0] for r in results])
        alogp_u = np.concatenate([r[1] for r in results])
        per_unique = {name: [] for name in FINGERPRINT_CLASSES}
        for r in results:
            for name in FINGERPRINT_CLASSES:
                per_unique[name].extend(r[2][name])
    else:
        mw_u, alogp_u, per_unique = _featurize_chunk(uniques, fp_format)

    out = df.copy()
    # Map unique-level results back to every row via the factorize codes
    # (positional assignment — independent of the DataFrame index).
    out["MW"] = mw_u[codes]
    out["ALOGP"] = alogp_u[codes]
    for name in FINGERPRINT_CLASSES:
        pu = per_unique[name]
        out[name] = [pu[c] for c in codes]

    return out
