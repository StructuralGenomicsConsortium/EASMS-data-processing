import pandas as pd
import numpy as np

def generate_ml_labels(df):
    """
    Assigns AIRCHECK_LABEL based on EASMS_ENRICHMENT, PVALUE, ISOMERS, and HAD_DUPLICATE_INTENSITY:

    - AIRCHECK_LABEL = 3: if EASMS_ENRICHMENT ≥ 5 and PVALUE ≤ 0.05 and ISOMER is not empty
    - AIRCHECK_LABEL = 2: if 5 ≤ EASMS_ENRICHMENT < 10 and PVALUE ≤ 0.05
    - AIRCHECK_LABEL = 1: if EASMS_ENRICHMENT ≥ 10 and PVALUE ≤ 0.05
    - AIRCHECK_LABEL = 0: if 0 ≤ EASMS_ENRICHMENT ≤ 1 or PVALUE > 0.05
    - AIRCHECK_LABEL = -1: if 1 < EASMS_ENRICHMENT < 5 and PVALUE ≤ 0.05
    - AIRCHECK_LABEL = -2: if EASMS_ENRICHMENT is missing
    - AIRCHECK_LABEL = 4: if HAD_DUPLICATE_INTENSITY == "Y" and ENRICHMENT > 5
    """

    required_columns = {"EASMS_ENRICHMENT", "PVALUE", "ISOMERS"}
    if not required_columns.issubset(df.columns):
        raise ValueError(f"Missing required columns: {required_columns - set(df.columns)}")

    # Clean up values
    df["EASMS_ENRICHMENT"].replace("", np.nan, inplace=True)
    df["PVALUE"].replace("", np.nan, inplace=True)

    # Convert to numeric
    df["EASMS_ENRICHMENT"] = pd.to_numeric(df["EASMS_ENRICHMENT"], errors="coerce")
    df["PVALUE"] = pd.to_numeric(df["PVALUE"], errors="coerce")

    def assign_label(row):
        enrichment = row["EASMS_ENRICHMENT"]
        pvalue = row["PVALUE"]
        isomer = str(row["ISOMERS"]).strip()

        if pd.isna(enrichment):
            return -2
        elif enrichment >= 5 and pvalue <= 0.05 and isomer != 'nan' and isomer != "":
            return 3
        elif 5 <= enrichment < 10 and pvalue <= 0.05:
            return 2
        elif enrichment >= 10 and pvalue <= 0.05:
            return 1
        elif 0 <= enrichment <= 1 or pvalue > 0.05:
            return 0
        elif 1 < enrichment < 5 and pvalue <= 0.05:
            return -1
        else:
            return -2  # fallback

    # Assign labels
    df["AIRCHECK_LABEL"] = df.apply(assign_label, axis=1).astype("int8")

    # Apply the NA rule for high enrichment and duplicate intensity
    if "HAD_DUPLICATE_INTENSITY" in df.columns:
        mask = (df["HAD_DUPLICATE_INTENSITY"] == "Y") & (df["EASMS_ENRICHMENT"] > 5)
        df.loc[mask, "AIRCHECK_LABEL"] = 4

    return df
