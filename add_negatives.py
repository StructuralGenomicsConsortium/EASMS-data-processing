import os
import pandas as pd
import numpy as np

def add_negative_samples_from_masterlist(df, file_name, masterlist_path, MasterList_Information):
    """
    Adds negative samples from the master list that are not present in the input DataFrame.
    Copies specific additional columns from the master list.

    Args:
        df (pd.DataFrame): The input DataFrame.
        file_name (str): The name of the processed file to match with the master list.
        masterlist_path (str): The directory containing master list files.
        MasterList_Information (str): The path to the Excel file mapping file names to master lists.

    Returns:
        pd.DataFrame: Updated DataFrame with added negative samples.
    """

    # Load the master list mapping file
    masterlist_info = pd.read_excel(MasterList_Information)

    # Ensure the necessary columns exist
    if not {"FileName", "MaterListName"}.issubset(masterlist_info.columns):
        raise ValueError("MasterList_Information.xlsx must contain 'FileName' and 'MaterListName' columns")

    # Get the corresponding master list name for the given file
    masterlist_name = masterlist_info.loc[masterlist_info["FileName"] == file_name, "MaterListName"]

    if masterlist_name.empty:
        print(f"Warning: No master list found for {file_name}. Skipping negative sample addition.")
        return df

    masterlist_name = masterlist_name.values[0]  # Extract string value

    # Construct the full path to the master list file
    masterlist_file = os.path.join(masterlist_path, f"{masterlist_name}.xlsx")

    # Check if the master list file exists
    if not os.path.exists(masterlist_file):
        print(f"Warning: Master list file {masterlist_file} not found. Skipping negative sample addition.")
        return df

    # Load the master list file
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

    # Assign BINARY_LABEL = "N" for these new negative samples
    new_entries["BINARY_LABEL"] = "N"
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
