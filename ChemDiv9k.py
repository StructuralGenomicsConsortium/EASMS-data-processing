import pandas as pd

# Read the Excel file
input_file = r"D:\0000-UHN\03-DataAndCodes\Data\ASMS\EASMS-7March\MasterLists\Chemdiv+Chiral6k_15k.xlsx"  # Change to your actual file name
df = pd.read_excel(input_file)

# Select the first 9007 rows
df_selected = df.iloc[:9007]

# Write the new file
output_file = r"D:\0000-UHN\03-DataAndCodes\Data\ASMS\EASMS-7March\MasterLists\Chemdiv_9k.xlsx"
df_selected.to_excel(output_file, index=False)

print(f"Filtered file saved as {output_file}")