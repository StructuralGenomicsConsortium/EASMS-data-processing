# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 10:10:24 2025

@author: shay
"""

import os
import pandas as pd

def split_protein_data(file_path, subfolder):
    """
    Splits a CSV file based on the 'PROTEIN_NUMBER' column and saves each part separately.
    
    Args:
        file_path (str): Path to the input CSV file.
        subfolder (str): Path to the folder where separated files should be stored.
        
    Returns:
        list: List of file paths for the separated CSV files.
    """
    # Load the CSV file
    print(file_path)

    df = pd.read_csv(file_path)

    # Ensure necessary columns exist
    required_columns = {"PROTEIN_NUMBER", "ASMS_BATCH_NUM", "TARGET_ID"}
    if not required_columns.issubset(df.columns):
        raise ValueError(f"Missing required columns: {required_columns - set(df.columns)}")

    # Dictionary to store output file paths
    separated_files = []

    # Group by PROTEIN_NUMBER and process each group
    for protein_number, group_df in df.groupby("TARGET_ID"):
        # Extract batch number (NUM) and protein name
        batch_number = group_df["ASMS_BATCH_NUM"].iloc[0]  # Take the first batch number
        protein_name = group_df["TARGET_ID"].iloc[0]  # Take the first protein name

        # Construct the filename
        file_name = f"{protein_name}_AsmBatchNumber{batch_number}.csv"

        # Define the file path
        output_path = os.path.join(subfolder, file_name)

        # Save the separated file
        group_df.to_csv(output_path, index=False) 
        separated_files.append(output_path)

        print(f"  Saved separated file: {output_path}")

    return separated_files
