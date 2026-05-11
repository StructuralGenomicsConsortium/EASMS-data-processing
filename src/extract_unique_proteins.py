#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
"""
Created on Wed Mar 19 14:58:58 2025

@author: shagh
"""
"""
import os
import pandas as pd

def extract_unique_protein_values(folder_path):
    # Get all CSV files in the folder
    csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]

    if not csv_files:
        print("No CSV files found in the folder.")
        return

    for file in csv_files:
        file_path = os.path.join(folder_path, file)
        
        try:
            # Read the CSV file
            df = pd.read_csv(file_path)
            
            # Extract unique values if the columns exist
            if 'PROTEIN_NUMBER' in df.columns:
                unique_protein_numbers = df['PROTEIN_NUMBER'].dropna().unique().tolist()
                # Count occurrences of each unique PROTEIN_NUMBER
                protein_number_counts = df['PROTEIN_NUMBER'].value_counts().to_dict()
                
            else:
                unique_protein_numbers = []
                protein_number_counts = {}

            #unique_protein_ids = df['PROTEIN_ID'].dropna().unique().tolist() if 'PROTEIN_ID' in df.columns else []
            unique_protein_ids = df['TARGET_ID'].dropna().unique().tolist() if 'TARGET_ID' in df.columns else []

            # Print results
            print(f"\nFile: {file}")
            print(f"Unique PROTEIN_NUMBER values: {unique_protein_numbers}")
            print(f"Number of rows per PROTEIN_NUMBER: {list(protein_number_counts.values())}")
            print(f"Unique PROTEIN_ID values: {unique_protein_ids}")
        
        except Exception as e:
            print(f"Error reading {file}: {e}")




# Example usage
folder_path = r"D:\0000-UHN\03-DataAndCodes\Data\ASMS\EASMS-7March\RawData"  # Replace with your actual folder path
extract_unique_protein_values(folder_path)
"""
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
import os
import pandas as pd

#### Chack data

def extract_unique_protein_values(folder_path):
    # Get all CSV files in the folder
    csv_files = [f for f in os.listdir(folder_path) if f.endswith('.csv')]

    if not csv_files:
        print("No CSV files found in the folder.")
        return

    for file in csv_files:
        file_path = os.path.join(folder_path, file)
        
        try:
            # Read the CSV file
            df = pd.read_csv(file_path)
            
            # Extract unique PROTEIN_NUMBER values
            '''if 'PROTEIN_NUMBER' in df.columns:
                unique_protein_numbers = df['PROTEIN_NUMBER'].dropna().unique().tolist()
            else:
                unique_protein_numbers = []'''

            # Extract unique TARGET_ID values and count rows per TARGET_ID
            if 'TARGET_ID' in df.columns:
                unique_protein_ids = df['TARGET_ID'].dropna().unique().tolist()
                target_id_counts = df['TARGET_ID'].value_counts().to_dict()
            else:
                unique_protein_ids = []
                target_id_counts = {}

            # Print results
            print(f"\nFile: {file}")
            # print(f"Unique PROTEIN_NUMBER values: {unique_protein_numbers}")
            print(f"Number of rows per TARGET_ID: {target_id_counts}")
            print(f"Unique TARGET_ID values: {unique_protein_ids}")
        
        except Exception as e:
            print(f"Error reading {file}: {e}")



# Example usage
folder_path = r"D:\0000-UHN\03-DataAndCodes\Data\ASMS\EASMS-7March\RawData"  # Replace with your actual folder path
extract_unique_protein_values(folder_path)




































