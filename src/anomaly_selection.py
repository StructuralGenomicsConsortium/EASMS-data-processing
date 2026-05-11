import os
import pandas as pd
import warnings

# Suppress FutureWarnings for pandas operations
# warnings.simplefilter(action='ignore', category=FutureWarning)

def filter_anomalous_data(df, sep_file_name):
    """
    Filters duplicate rows and processes SMILES with different ENRICHMENT values:
    - If all rows for a SMILES have ENRICHMENT < 1, keeps only the row with the smallest ENRICHMENT.
    - If all rows for a SMILES have ENRICHMENT > 10, keeps only the row with the highest ENRICHMENT.
    - If the subset contains mixed enrichment values, removes all its rows as it is confusing.
    - Logs all conflicting SMILES (before filtering) with all columns.
    - Logs removed SMILES due to mixed enrichment.

    Args:
        df (pd.DataFrame): The input DataFrame.
        sep_file_name (str): The name of the separated CSV file being processed.

    Returns:
        pd.DataFrame: Cleaned DataFrame with anomalies handled.
    """

    # Ensure the necessary columns exist
    required_columns = {"SMILES", "ENRICHMENT"}
    if not required_columns.issubset(df.columns):
        raise ValueError(f"Missing required columns: {required_columns - set(df.columns)}")

    # Step 1: Remove fully duplicate rows
    df_cleaned = df.drop_duplicates()

    # Step 2: Identify SMILES that have multiple ENRICHMENT values
    enrichment_groups = df_cleaned.groupby("SMILES")["ENRICHMENT"].nunique()
    conflicting_smiles = enrichment_groups[enrichment_groups > 1].index.tolist()

    # Prepare log dataframe with all conflicting SMILES before filtering
    conflict_log_df = df_cleaned[df_cleaned["SMILES"].isin(conflicting_smiles)].copy()

    # Prepare DataFrames for keeping and removing rows
    rows_to_keep_df = pd.DataFrame(columns=df_cleaned.columns)
    removed_df = pd.DataFrame(columns=df_cleaned.columns)  # Store removed conflicting SMILES

    '''for smiles in conflicting_smiles:
        subset = df_cleaned[df_cleaned["SMILES"] == smiles]  # Get all rows for this SMILES

        if subset["ENRICHMENT"].max() < 1:  
            # All values are < 1 → Keep the row with the lowest ENRICHMENT
            best_row = subset.loc[[subset["ENRICHMENT"].idxmin()]]
            if not best_row.empty:
                rows_to_keep_df = pd.concat([rows_to_keep_df, best_row], ignore_index=True)
        elif subset["ENRICHMENT"].min() > 5:  
            # All values are > 5 → Keep the row with the highest ENRICHMENT
            best_row = subset.loc[[subset["ENRICHMENT"].idxmax()]]
            if not best_row.empty:
                rows_to_keep_df = pd.concat([rows_to_keep_df, best_row], ignore_index=True)
        else:
            # Mixed values → Remove entire subset & log it
            if not subset.empty:
                removed_df = pd.concat([removed_df, subset], ignore_index=True)

    # Step 3: Merge back with non-conflicting SMILES
    final_df = pd.concat([df_cleaned[~df_cleaned["SMILES"].isin(conflicting_smiles)], rows_to_keep_df], ignore_index=True)

    # Step 4: Save the log file (all conflicting SMILES + removed SMILES)
    log_dir = os.getcwd()  # Change this if you want a custom folder
    log_file_path = os.path.join(log_dir, f"Conflicting_SMILES_Log_{sep_file_name}.csv")

    # Merge all conflicting SMILES (before filtering) and removed rows
    full_log_df = pd.concat([conflict_log_df, removed_df], ignore_index=True).drop_duplicates()

    if not full_log_df.empty:
        full_log_df.to_csv(log_file_path, index=False)
        print(f"Logged {len(full_log_df)} conflicting SMILES entries in {log_file_path}")
                
            '''
                
    for smiles in conflicting_smiles:
        subset = df_cleaned[df_cleaned["SMILES"] == smiles]  # Get all rows for this SMILES

        if subset["EASMS_ENRICHMENT"].max() <= 1:  
            # All values are < 1 → Keep the row with the lowest ENRICHMENT
            best_row = subset.loc[[subset["EASMS_ENRICHMENT"].idxmin()]]
            if not best_row.empty:
                rows_to_keep_df = pd.concat([rows_to_keep_df, best_row], ignore_index=True)
        elif subset["EASMS_ENRICHMENT"].min() > 1:  
            # All values are > 5 → Keep the row with the highest ENRICHMENT
            best_row = subset.loc[[subset["EASMS_ENRICHMENT"].idxmax()]]
            if not best_row.empty:
                rows_to_keep_df = pd.concat([rows_to_keep_df, best_row], ignore_index=True)

    # Add HAD_DUPLICATE_INTENSITY column
    df_cleaned.loc[~df_cleaned["SMILES"].isin(conflicting_smiles), "HAD_DUPLICATE_INTENSITY"] = "N"
    rows_to_keep_df["HAD_DUPLICATE_INTENSITY"] = "Y"

    # Step 3: Merge back with non-conflicting SMILES
    final_df = pd.concat([df_cleaned[~df_cleaned["SMILES"].isin(conflicting_smiles)], rows_to_keep_df], ignore_index=True)


    return final_df
