# -*- coding: utf-8 -*-
"""
Created on Mon Mar  3 11:16:42 2025

@author: shagh
"""

import pandas as pd
from rdkit import Chem

# File paths
csv_file = r"D:\0000-UHN\03-DataAndCodes\AIRCHECK-workflow\SimpleML\Bootcamp\ExtractFingerPrints_Bootcamp\58Hits.csv"
excel_file = r"D:\0000-UHN\03-DataAndCodes\AIRCHECK-workflow\SimpleML\Bootcamp\Data\3March\250228_WDR91_70_ASMS_hits_SPR_Confirmed-v2.xlsx"
output_file = r"D:\0000-UHN\03-DataAndCodes\AIRCHECK-workflow\SimpleML\Bootcamp\Data\3March\70HitsMinus58Hits.csv"

# Load data
csv_df = pd.read_csv(csv_file, dtype=str)
excel_df = pd.read_excel(excel_file, dtype=str)

# Ensure consistent column name
smiles_col = "SMILES (Compounds)"

if smiles_col not in csv_df.columns or smiles_col not in excel_df.columns:
    raise ValueError(f"Column '{smiles_col}' not found in one of the files.")

# Function to convert SMILES to canonical form
def canonicalize_smiles(smiles):
    try:
        return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), canonical=True)
    except:
        return None

# Convert SMILES to canonical form
csv_df[smiles_col] = csv_df[smiles_col].dropna().apply(canonicalize_smiles)
excel_df[smiles_col] = excel_df[smiles_col].dropna().apply(canonicalize_smiles)

# Extract unique SMILES from both files
csv_smiles_set = set(csv_df[smiles_col].dropna().unique())
excel_filtered_df = excel_df[~excel_df[smiles_col].isin(csv_smiles_set)]

# Save the filtered rows to a new CSV file
excel_filtered_df.to_csv(output_file, index=False)

print(f"Filtered data saved to {output_file}")