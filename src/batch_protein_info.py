# -*- coding: utf-8 -*-
"""Per-batch protein summary.

Writes one CSV per input file listing the distinct proteins (one row per
TARGET_ID) in that batch, with the protein-identity columns. For a typical
batch of 8 proteins this produces an 8-row (plus header) file.
"""

import pandas as pd
from io_utils import pjoin, makedirs

# Columns included in the summary, in output order.
BATCH_PROTEIN_COLUMNS = ["PROTEIN_NUMBER", "UNIPROT_ID", "PROTEIN_SEQ", "TARGET_ID"]
OUTPUT_FILENAME = "proteins_in_batch.csv"


def write_batch_protein_info(file_path, output_dir):
    """Read the raw CSV at `file_path`, keep one row per distinct TARGET_ID, and
    write the protein summary to `<output_dir>/proteins_in_batch.csv`.

    Only the columns in BATCH_PROTEIN_COLUMNS that exist in the file are
    included. Rows are deduplicated by TARGET_ID and ordered by PROTEIN_NUMBER
    when available. Returns the output path, or None if it could not be produced
    (file unreadable or no TARGET_ID column).
    """
    try:
        header = pd.read_csv(file_path, nrows=0).columns.tolist()
    except Exception as e:
        print(f"  Warning: could not read {file_path} for batch protein info: {e}")
        return None

    cols = [c for c in BATCH_PROTEIN_COLUMNS if c in header]
    if "TARGET_ID" not in cols:
        print(f"  Warning: no TARGET_ID column in {file_path}; skipping batch protein info.")
        return None

    df = pd.read_csv(file_path, usecols=cols)
    df = df.drop_duplicates(subset=["TARGET_ID"])
    if "PROTEIN_NUMBER" in cols:
        df = df.sort_values("PROTEIN_NUMBER", kind="stable")
    df = df[cols]  # enforce the requested column order

    makedirs(output_dir, exist_ok=True)
    out_path = pjoin(output_dir, OUTPUT_FILENAME)
    df.to_csv(out_path, index=False)
    print(f"  Wrote batch protein info ({len(df)} proteins) -> {OUTPUT_FILENAME}")
    return out_path
