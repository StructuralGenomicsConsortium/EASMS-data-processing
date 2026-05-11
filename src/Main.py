# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 09:24:34 2025

@author: shay
test
"""



import os
import argparse
import pandas as pd
from separate_protein_files import split_protein_data
from add_scores import compute_and_add_scores
from anomaly_selection import filter_anomalous_data
from isomer_handling import handle_isomers
from produce_ml_labels import generate_ml_labels
from add_negatives import add_negative_samples_from_masterlist
from fingerprint_extraction import extract_fingerprints
from column_selection import select_final_columns

def process_csv_files(data_path, masterlist_path, separated_files_dir, output_dir1, output_dir2, output_dir3, MasterList_Information, DesiredColumns):
    """Processes all CSV files through data curation steps."""

    # Step 1: Separate protein-related data and store in subfolders
    for file_name in os.listdir(data_path):
        if file_name.endswith(".csv"):
            file_path = os.path.join(data_path, file_name)
            print(f"Processing: {file_name}")

            # Create a subfolder for the separated files
            subfolder = os.path.join(separated_files_dir, os.path.splitext(file_name)[0])
            os.makedirs(subfolder, exist_ok=True)

            # Step 1: Split protein files
            separated_files = split_protein_data(file_path, subfolder)
            
            # Step 2: Compute and Add Scores to all separated files together
            print("\nComputing and Adding Scores to All Separated Files...\n")
            compute_and_add_scores(separated_files)  # Process scores in a batch  


            # Step 3-8: Process each separated file after computing scores
            for sep_file in separated_files:
                sep_file_name = os.path.basename(sep_file)
                base_name = os.path.splitext(sep_file_name)[0]
                print(f"  Processing separated file: {sep_file_name}")
        
                # Load separated file
                df = pd.read_csv(sep_file)
        
                # Step 3: Identify and filter out anomalies
                df = filter_anomalous_data(df,sep_file_name)
        
                # Step 4: Handle isomer-specific corrections
                df = handle_isomers(df,sep_file_name)
        
                # Step 5: Add additional negative samples from master list
                df = add_negative_samples_from_masterlist(df, file_name, masterlist_path,MasterList_Information)
                
                # Step 6: Generate ML labels
                df = generate_ml_labels(df)
        
                # Save curated CSV file (MLReady)
                output_file1_csv = os.path.join(output_dir1, f"MLReady_{base_name}.csv")
                df.to_csv(output_file1_csv, index=False)

                output_file1_parquet = os.path.join(output_dir1, f"MLReady_{base_name}.parquet")
                df.to_parquet(output_file1_parquet, index=False)
                
                print(f"  Saved intermediate file: {output_file1_csv}")
        
                # Step 7: Extract chemical fingerprints
                df = extract_fingerprints(df) 
        
        
                # Step 8: 
                # RenameColumns
                df = df.rename(columns={"TARGET_VALUE": "TARGET_INTENSITY_VALUE"})
                df = df.rename(columns={"MEAN_NONTARGET_VALUES": "NONTARGET_INTENSITY_VALUE"})
                
                # Creates a new column LABEL with 1 if BINARY_LABEL is "Y", and 0 if it's "N":
                df["LABEL"] = (df["BINARY_LABEL"] == "Y").astype(int)
                
                # Step 9: Select final columns
                df = select_final_columns(df, DesiredColumns)                
       
                # Save as CSV
                output_file2_csv = os.path.join(output_dir2, f"MLReadyPlusFPs_{base_name}.csv")
                df.to_csv(output_file2_csv, index=False)
                
                # Save as Parquet
                output_file2_parquet = os.path.join(output_dir2, f"MLReadyPlusFPs_{base_name}.parquet")
                df.to_parquet(output_file2_parquet, index=False)
                
                print(f"Saved CSV: {output_file2_csv}")
                print(f"Saved Parquet: {output_file2_parquet}")
                #-----------------------
                df = select_final_columns(df, DesiredColumns2)                
       
                # Save as CSV
                output_file3_csv = os.path.join(output_dir3, f"MLReadyPlusFPs_{base_name}.csv")
                df.to_csv(output_file3_csv, index=False)
                
                # Save as Parquet
                output_file3_parquet = os.path.join(output_dir3, f"MLReadyPlusFPs_{base_name}.parquet")
                df.to_parquet(output_file3_parquet, index=False)
                print(f"Saved CSV: {output_file3_csv}")
                print(f"Saved Parquet: {output_file3_parquet}")
                


def main(data_path, masterlist_path, separated_files_dir, output_dir1, output_dir2, output_dir3, MasterList_Information, DesiredColumns):
    """Main function to execute the full data curation pipeline."""
    os.makedirs(separated_files_dir, exist_ok=True)
    os.makedirs(output_dir1, exist_ok=True)
    os.makedirs(output_dir2, exist_ok=True)
    os.makedirs(output_dir3, exist_ok=True)

    process_csv_files(data_path, masterlist_path, separated_files_dir, output_dir1, output_dir2, output_dir3, MasterList_Information, DesiredColumns)

if __name__ == "__main__":
    # Define paths (Modify as needed)
    parser = argparse.ArgumentParser(description="Run EASMS data processing pipeline.")
    parser.add_argument(
        "--path",
        default=os.path.abspath(os.path.join(os.getcwd(), "..")),
        help="Dataset root directory containing RawData/ and MasterLists/ (default: parent of cwd).",
    )
    args = parser.parse_args()
    path = args.path
    data_path = os.path.join(path, "RawData")
    masterlist_path = os.path.join(path, "MasterLists")  
    MasterList_Information = os.path.join(masterlist_path, "MasterList_Information.xlsx")   

    processed_data_dir = os.path.join(path, "ProcessedData")
    separated_files_dir = os.path.join(processed_data_dir, "Separated_Files")
    output_dir1 = os.path.join(processed_data_dir, "MLReady")
    output_dir2 = os.path.join(processed_data_dir, "MLReady_FullColumns")
    output_dir3 = os.path.join(processed_data_dir, "MLReady_KeyColumns")

    DesiredColumns = ['ASMS_BATCH_NUM',
     'COMPOUND_ID',
     'COMPOUND_FORMULA',
     'SMILES',
     'POOL_NAME',
     'POOL_ID',
     'POOL_SIZE',
     'PROTEIN_NUMBER',
     'TARGET_ID',
     'PROTEIN_ID',
     'PROTEIN_SEQ',
     'PROTEIN_TAG',
     'INCUBATION_VOLUME',
     'PROTEIN_CONC',
     'COMPOUND_CONC',
     'MS_REPRODUCABILITY',
     'POS_INT_REP1',
     'POS_INT_REP2',
     'POS_INT_REP3',
     'TARGET_INTENSITY_VALUE',
     'SELECTIVE_VALUE',
     'NTC_VALUE',
     'ENRICHMENT',
     'SELECTIVE_ENRICHMENT',
     'PVALUE',
     'BINARY_LABEL',
     'HAD_DUPLICATE_INTENSITY',
     'SPR',
     'KD',
     'ISOMERS',
     'MassSpec_Detected',
     'EASMS_ENRICHMENT',
     'NONTARGET_INTENSITY_VALUE',
     'LABEL',
     'AIRCHECK_LABEL',
     'MW',
     'ALOGP',
     'ECFP4',
     'ECFP6',
     'FCFP4',
     'FCFP6',
     'MACCS',
     'RDK',
     'AVALON',
     'TOPTOR',
     'ATOMPAIR']
    
    DesiredColumns2 = [
     'COMPOUND_ID',
     'SMILES',
     'TARGET_ID',
     'TARGET_INTENSITY_VALUE',
     'NONTARGET_INTENSITY_VALUE',
     'EASMS_ENRICHMENT',
     'PVALUE',
     'LABEL',
     'MW',
     'ALOGP',
     'ECFP4',
     'ECFP6',
     'FCFP4',
     'FCFP6',
     'MACCS',
     'RDK',
     'AVALON',
     'TOPTOR',
     'ATOMPAIR']

    main(data_path, masterlist_path, separated_files_dir, output_dir1, output_dir2, output_dir3, MasterList_Information, DesiredColumns)
