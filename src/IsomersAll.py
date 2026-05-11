import pandas as pd
import os

# Folder path
folder_path = r"D:\0000-UHN\03-DataAndCodes\Data\ASMS\LogsForSunny\A"

# Get full paths of all files in the folder
all_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]

# List to store successfully read DataFrames
dataframes = []

# Try reading each file with pandas
for file in all_files:
    try:
        df = pd.read_csv(file)
        dataframes.append(df)
        print(f"Loaded: {os.path.basename(file)}")
    except Exception as e:
        print(f"Skipped: {os.path.basename(file)} — {e}")

# Combine and save
if dataframes:
    combined_df = pd.concat(dataframes, ignore_index=True)
    output_path = os.path.join(folder_path, "Isomers_All.csv")
    combined_df.to_csv(output_path, index=False)
    print(f"\nCombined file saved as: {output_path}")
else:
    print("No valid CSV files found.")
