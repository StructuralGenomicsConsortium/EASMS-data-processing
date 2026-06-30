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


def load_column_tags(path):
    """Read a column/action spreadsheet and return (data_columns, metadata_columns).

    The spreadsheet has a column of raw/computed column names and a column of
    actions. Names tagged ``data`` go to the final data file, names tagged
    ``metadata`` go to the metadata file, and any other tag (e.g. ``-``) is
    ignored. Headers are matched case-insensitively as ``Column name`` and
    ``Action``, falling back to the first two columns. Accepts local paths and
    ``gs://`` URLs (via fsspec).

    Args:
        path (str): Path to the .xlsx (or .csv) tag file.

    Returns:
        tuple[list[str], list[str]]: (data_columns, metadata_columns), in the
        order they appear in the file.
    """
    if str(path).lower().endswith(".csv"):
        tags = pd.read_csv(path)
    else:
        tags = pd.read_excel(path)

    lower = {str(c).strip().lower(): c for c in tags.columns}
    name_col = lower.get("column name", tags.columns[0])
    action_col = lower.get("action", tags.columns[1])

    data_columns, metadata_columns = [], []
    for name, action in zip(tags[name_col], tags[action_col]):
        name = str(name).strip()
        action = str(action).strip().lower()
        if not name or name.lower() == "nan":
            continue
        if action == "data":
            data_columns.append(name)
        elif action == "metadata":
            metadata_columns.append(name)
    return data_columns, metadata_columns
