# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 12:14:55 2025

@author: shay
"""

import os
import pandas as pd

def handle_isomers(df, sep_file_name):
    """
    Splits rows containing multiple isomers (separated by ";") into separate rows and logs those with ENRICHMENT > 10.

    Args:
        df (pd.DataFrame): The input DataFrame.
        sep_file_name (str): The name of the separated CSV file being processed.

    Returns:
        pd.DataFrame: Expanded DataFrame with individual isomer rows and a new "ISOMERS" column.
    """

    # Ensure necessary columns exist
    required_columns = {"COMPOUND_ID", "COMPOUND_FORMULA", "SMILES", "ENRICHMENT"}
    if not required_columns.issubset(df.columns):
        raise ValueError(f"Missing required columns: {required_columns - set(df.columns)}")

    # Identify rows with isomers (rows where SMILES contains ";")
    isomer_rows = df[df["SMILES"].str.contains(";", na=False)].copy()

    # Storage for new expanded rows
    expanded_rows = []
    log_rows = []

    # Process each row with isomers
    for _, row in isomer_rows.iterrows():
        # Split isomer-related columns into lists
        compound_ids = row["COMPOUND_ID"].split(";")
        compound_formulas = row["COMPOUND_FORMULA"].split(";")
        smiles_list = row["SMILES"].split(";")

        # Ensure lists have the same length
        if not (len(compound_ids) == len(compound_formulas) == len(smiles_list)):
            raise ValueError(f"Inconsistent isomer data in row: {row}")

        # Create new rows for each isomer
        for i, (comp_id, comp_formula, smile) in enumerate(zip(compound_ids, compound_formulas, smiles_list)):
            new_row = row.copy()  # Copy original row
            new_row["COMPOUND_ID"] = comp_id
            new_row["COMPOUND_FORMULA"] = comp_formula
            new_row["SMILES"] = smile
            new_row["ISOMERS"] = ";".join([x for j, x in enumerate(compound_ids) if j != i])  # Store other isomers

            '''# Add to log if ENRICHMENT > 5
            if new_row["EASMS_ENRICHMENT"] > 5 and new_row["PVALUE"] < 0.05:
                #log_rows.append(new_row[["COMPOUND_ID", "COMPOUND_FORMULA", "SMILES", "ENRICHMENT", "ISOMERS"]])
                log_rows.append(new_row)'''

            expanded_rows.append(new_row)

    # Convert expanded rows into a DataFrame
    expanded_df = pd.DataFrame(expanded_rows)

    # Remove original isomer-containing rows and append the new ones
    df = df[~df["SMILES"].str.contains(";", na=False)]
    df = pd.concat([df, expanded_df], ignore_index=True)
    
    # Ensure ISOMERS column exists and is filled
    if "ISOMERS" not in df.columns:
        df["ISOMERS"] = ""
    else:
        df["ISOMERS"] = df["ISOMERS"].fillna("")

    '''# Save log file for ENRICHMENT > 10 cases
    if log_rows:
        log_df = pd.DataFrame(log_rows)
        log_file_path = os.path.join(os.getcwd(), f"IsomersLog_{sep_file_name}.csv")
        log_df.to_csv(log_file_path, index=False)
        print(f"Logged {len(log_df)} isomer entries with ENRICHMENT > 10 in {log_file_path}")'''

    return df
