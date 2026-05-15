import os
import pandas as pd
import numpy as np
import scipy.stats as stats

def compute_and_add_scores(file_paths):
    """
    Computes TARGET_VALUE, ENRICHMENT, SELECTIVE_ENRICHMENT, EASMS_ENRICHMENT,
    MEAN_NONTARGET_VALUES, and PVALUE for a list of CSV files.
    Each file is processed individually, using all other files to compute comparison values.
    """
    if not file_paths:
        print("No files provided for score computation.")
        return

    # Load all CSV files
    dataframes = {f: pd.read_csv(f) for f in file_paths}

    # Ensure necessary columns exist
    required_columns = {"COMPOUND_ID", "POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"}
    for filename, df in dataframes.items():
        if not required_columns.issubset(df.columns):
            raise ValueError(f"Missing required columns in {filename}: {required_columns - set(df.columns)}")

    # Compute TARGET_VALUE for each df
    for df in dataframes.values():
        df["TARGET_VALUE"] = pd.to_numeric(
            df[["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]]
            .mean(axis=1, skipna=True), errors="coerce"
        )

    # Step 2: Process each file individually
    for current_file, df in dataframes.items():
        print(f"\n Processing: {os.path.basename(current_file)}")

        # Exclude current df
        other_dfs = [d for f, d in dataframes.items() if f != current_file]
        merged_other = pd.concat(other_dfs, ignore_index=True)

        # Ensure TARGET_VALUE is present in merged_other
        if "TARGET_VALUE" not in merged_other.columns:
            merged_other["TARGET_VALUE"] = pd.to_numeric(
                merged_other[["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]]
                .mean(axis=1, skipna=True), errors="coerce"
            )

        # SELECTIVE_VALUE & NTC_VALUE
        max_values = merged_other.groupby("COMPOUND_ID")["TARGET_VALUE"].max()
        min_values = merged_other.groupby("COMPOUND_ID")["TARGET_VALUE"].min()
        df["SELECTIVE_VALUE"] = df["COMPOUND_ID"].map(max_values)
        df["NTC_VALUE"] = df["COMPOUND_ID"].map(min_values)

        # ENRICHMENT calculations
        df["ENRICHMENT"] = df["TARGET_VALUE"] / df["NTC_VALUE"]
        df["SELECTIVE_ENRICHMENT"] = df["TARGET_VALUE"] / df["SELECTIVE_VALUE"]

        # EASMS_ENRICHMENT and MEAN_NONTARGET_VALUES
        def compute_easms_enrichment(row):
            cid = row["COMPOUND_ID"]
            other = merged_other[merged_other["COMPOUND_ID"] == cid]
            if other.empty:
                return pd.Series([None, None])
            other_mean = other[["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]]\
                .apply(pd.to_numeric, errors="coerce")\
                .mean(axis=1, skipna=True).mean()
            if pd.isna(other_mean) or other_mean == 0:
                return pd.Series([None, other_mean])
            enrichment = row["TARGET_VALUE"] / other_mean
            return pd.Series([enrichment, other_mean])

        df[["EASMS_ENRICHMENT", "MEAN_NONTARGET_VALUES"]] = df.apply(compute_easms_enrichment, axis=1)

        # PVALUE computation (robust)
        def calculate_p_value(row):
            try:
                compound_id = row["COMPOUND_ID"]
                if pd.isna(compound_id):
                    return None

                # Current compound's replicates
                protein_interest = pd.to_numeric(
                    row[["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]],
                    errors="coerce"
                ).dropna().values.astype(float)

                # Replicates from other files for same compound
                other = merged_other[merged_other["COMPOUND_ID"] == compound_id]
                protein_other_values = pd.to_numeric(
                    other[["POS_INT_REP1", "POS_INT_REP2", "POS_INT_REP3"]]
                    .stack(), errors="coerce"
                ).dropna().values.astype(float)

                if len(protein_interest) == 0 or len(protein_other_values) < 3:
                    return None

                if np.std(protein_interest) < 1e-8 or np.std(protein_other_values) < 1e-8:
                    return 1.0

                _, p_value = stats.ttest_ind(protein_interest, protein_other_values, equal_var=False)
                return p_value

            except Exception as e:
                print(f" Error calculating p-value for {row.get('COMPOUND_ID', 'UNKNOWN')}: {e}")
                return None

        df["PVALUE"] = df.apply(calculate_p_value, axis=1)

    # Step 3: Save the updated files
    for file_path, df in dataframes.items():
        df.to_csv(file_path, index=False)
        print(f"Updated and saved: {file_path}")

    print("\nAll files have been processed and saved with computed scores.")
