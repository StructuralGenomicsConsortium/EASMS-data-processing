# -*- coding: utf-8 -*-
"""
Created on Sun Mar  2 17:46:40 2025

@author: shagh
"""

import pandas as pd

def select_final_columns(df, DesiredColumns):
    """
    Selects only the specified columns from the DataFrame.

    Args:
        df (pd.DataFrame): The input DataFrame.
        DesiredColumns (list): A list of column names to keep in the final DataFrame.

    Returns:
        pd.DataFrame: A DataFrame containing only the selected columns.
    """

    # Get a list of existing columns in df
    existing_columns = df.columns.tolist()

    # Find columns that are available in the DataFrame
    available_columns = [col for col in DesiredColumns if col in existing_columns]

    # Warn if any DesiredColumns are missing
    missing_columns = set(DesiredColumns) - set(available_columns)
    if missing_columns:
        print(f"Warning: The following columns were not found in the DataFrame and will be ignored: {missing_columns}")

    # Select only the available columns
    df = df[available_columns]

    print(f"Final DataFrame contains {len(df.columns)} columns: {available_columns}")

    return df
