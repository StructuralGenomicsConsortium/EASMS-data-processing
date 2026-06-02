import os
import pandas as pd
import numpy as np

def add_negative_samples_from_masterlist(df, file_name, masterlist_path):
    """
    Adds negative samples from the master list that are not present in the input DataFrame.
    Copies specific additional columns from the master list.

    The master list (library) is resolved directly from the LIBRARY_NAME column
    of the input data, then loaded from `<library>.xlsx` or `<library>.csv` in
    `masterlist_path`. No separate mapping file is needed.

    Args:
        df (pd.DataFrame): The input DataFrame.
        file_name (str): The name of the processed file (used only for log messages).
        masterlist_path (str): The directory containing master list files.

    Returns:
        pd.DataFrame: Updated DataFrame with added negative samples.
    """

    # Resolve the library name straight from the LIBRARY_NAME column
    # (one file = one library).
    if "LIBRARY_NAME" not in df.columns:
        print(f"Warning: No LIBRARY_NAME column in {file_name}. Skipping negative sample addition.")
        return df

    library_values = df["LIBRARY_NAME"].dropna()
    if library_values.empty:
        print(f"Warning: LIBRARY_NAME column is empty in {file_name}. Skipping negative sample addition.")
        return df

    masterlist_name = str(library_values.iloc[0]).strip()

    # Resolve the master list file, accepting either .xlsx or .csv (xlsx wins
    # if both exist).
    masterlist_file = None
    for ext in (".xlsx", ".csv"):
        candidate = os.path.join(masterlist_path, f"{masterlist_name}{ext}")
        if os.path.exists(candidate):
            masterlist_file = candidate
            break

    if masterlist_file is None:
        print(f"Warning: Master list file '{masterlist_name}' (.xlsx or .csv) not found in {masterlist_path}. Skipping negative sample addition.")
        return df

    # Load the master list file with the reader matching its extension
    if masterlist_file.lower().endswith(".csv"):
        master_df = pd.read_csv(masterlist_file)
    else:
        master_df = pd.read_excel(masterlist_file)

    # Ensure 'SMILES' column exists in the master list
    if "SMILES" not in master_df.columns:
        raise ValueError(f"Master list file {masterlist_name} must contain a 'SMILES' column")

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
    
    
    # Define column mapping (Master List → df)
    column_mapping = {
        "SGC ID for Component": "COMPOUND_ID",  # Mapping 'SGC ID for Component' from master list to 'COMPOUND_ID' in df
        "SGC ID for Pool": "POOL_NAME",
        "formula": "COMPOUND_FORMULA"
        # Add more mappings if needed
    }

    # Apply column mappings
    for master_col, df_col in column_mapping.items():
        print("here")
        if master_col in new_entries.columns:
            new_entries[df_col] = new_entries[master_col]

    # Select only relevant columns to append
    columns_to_add = ["SMILES", "BINARY_LABEL","ENRICHMENT","PVALUE","MassSpec_Detected", "TARGET_ID"] + list(column_mapping.values())
    df = pd.concat([df, new_entries[columns_to_add]], ignore_index=True)

    print(f"Added {len(new_entries)} negative samples to {file_name} from {masterlist_name}, including columns: {list(column_mapping.values())}")

    return df
