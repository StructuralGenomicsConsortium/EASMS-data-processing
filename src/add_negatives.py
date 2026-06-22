import os
import pandas as pd
import numpy as np
from io_utils import pjoin, exists


def add_negative_samples_from_masterlist(df, file_name, masterlist_path):
    """
    Adds negative samples from the master list that are not present in the input DataFrame.
    Copies specific additional columns from the master list.

    The master list (library) is resolved directly from the LIBRARY_NAME column
    of the input data, then loaded from a `<library>.{xlsx,csv}` or
    `<library>_lib.{xlsx,csv}` file under `masterlist_path`. No separate mapping
    file is needed.

    Args:
        df (pd.DataFrame): The input DataFrame.
        file_name (str): The name of the processed file (used only for log messages).
        masterlist_path (str): A folder of master-list files OR a single
            master-list file, local or remote (`gs://`).

    Returns:
        pd.DataFrame: Updated DataFrame with added negative samples.
    """

    # Resolve the library name straight from the LIBRARY_NAME column
    # (one file = one library).
    if "LIBRARY_NAME" not in df.columns:
        print(
            f"Warning: No LIBRARY_NAME column in {file_name}. Skipping negative sample addition.")
        return df

    library_values = df["LIBRARY_NAME"].dropna()
    if library_values.empty:
        print(
            f"Warning: LIBRARY_NAME column is empty in {file_name}. Skipping negative sample addition.")
        return df

    masterlist_name = str(library_values.iloc[0]).strip()

    # Resolve the master list file. `masterlist_path` may be a folder of library
    # files OR a single library file, local or remote (gs://). Library files are
    # named "<library>_lib.<ext>" (e.g. EASMS12kV1_lib.csv) or, in older sets, a
    # bare "<library>.<ext>"; .xlsx wins over .csv.
    _exts = (".xlsx", ".xls", ".csv")
    path = masterlist_path.rstrip("/")
    masterlist_file = None
    if os.path.splitext(path)[1].lower() in _exts:
        # A single library file was given directly; strip a trailing "_lib"
        # marker and require its library name to match LIBRARY_NAME.
        stem = os.path.splitext(os.path.basename(path))[0]
        if stem.lower().endswith("_lib"):
            stem = stem[:-len("_lib")]
        if exists(path) and stem == masterlist_name:
            masterlist_file = path
    else:
        for stem in (masterlist_name, f"{masterlist_name}_lib"):
            for ext in _exts:
                candidate = pjoin(path, f"{stem}{ext}")
                if exists(candidate):
                    masterlist_file = candidate
                    break
            if masterlist_file:
                break

    if masterlist_file is None:
        print(
            f"Warning: Master list file for '{masterlist_name}' (.xlsx/.csv, optional '_lib' suffix) not found at {masterlist_path}. Skipping negative sample addition.")
        return df

    # Load the master list file with the reader matching its extension
    if masterlist_file.lower().endswith(".csv"):
        master_df = pd.read_csv(masterlist_file)
    else:
        master_df = pd.read_excel(masterlist_file)

    # Ensure 'SMILES' column exists in the master list
    if "SMILES" not in master_df.columns:
        raise ValueError(
            f"Master list file {masterlist_name} must contain a 'SMILES' column")

    # Identify SMILES that are NOT present in df
    existing_smiles = set(df["SMILES"].dropna())
    new_entries = master_df[~master_df["SMILES"].isin(existing_smiles)].copy()

    if new_entries.empty:
        print(f"No new negative samples found for {file_name}.")
        return df

    # Assign BINARY_LABEL = 0 for these new negative samples (integer 0/1 convention)
    new_entries["BINARY_LABEL"] = 0
    new_entries["ENRICHMENT"] = np.nan
    new_entries["PVALUE"] = np.nan
    new_entries["MassSpec_Detected"] = "N"
    new_entries["TARGET_ID"] = df["TARGET_ID"][0]
    df["MassSpec_Detected"] = "Y"

    # Rename master-list source columns to the df target names, but only when
    # the source column is present (master lists differ in which columns they
    # carry). 'formula' -> 'COMPOUND_FORMULA' is the common case; 'COMPOUND_ID'
    # is usually already named that way in newer master lists, so it carries
    # over directly. Legacy 'SGC ID for Component' is mapped if present.
    column_mapping = {
        "SGC ID for Component": "COMPOUND_ID",
        "formula": "COMPOUND_FORMULA",
        # POOL_NAME is intentionally NOT mapped: negative samples were not
        # tested, so they have no pool and POOL_NAME is left empty.
    }
    for master_col, df_col in column_mapping.items():
        if master_col in new_entries.columns:
            new_entries[df_col] = new_entries[master_col]

    # Carry over only the columns that actually exist on the negative rows, so a
    # master list missing an optional column never breaks the concat. Columns
    # not listed here (e.g. POOL_NAME) stay empty for negatives after concat.
    desired_cols = ["SMILES", "BINARY_LABEL", "ENRICHMENT", "PVALUE",
                    "MassSpec_Detected", "TARGET_ID", "COMPOUND_ID", "COMPOUND_FORMULA"]
    columns_to_add = [c for c in desired_cols if c in new_entries.columns]
    df = pd.concat([df, new_entries[columns_to_add]], ignore_index=True)

    print(f"Added {len(new_entries)} negative samples to {file_name} from {masterlist_name}, including columns: {columns_to_add}")

    return df
