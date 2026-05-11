# -*- coding: utf-8 -*-
"""
Created on Sun Mar  2 17:29:20 2025

@author: shagh
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from fingerprints import HitGenMACCS, HitGenECFP4, HitGenECFP6, HitGenFCFP4, HitGenFCFP6, HitGenRDK, HitGenAvalon, HitGenTopTor, HitGenAtomPair

def compute_molecular_properties(smiles):
    """
    Computes molecular properties (MW, ALOGP) for a given SMILES.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mw = Descriptors.MolWt(mol)
        alogp = Descriptors.MolLogP(mol)
    else:
        mw = np.nan
        alogp = np.nan
    return mw, alogp

def generate_fingerprints(smiles, fps_dict):
    """
    Generates fingerprints for a given SMILES string.

    Args:
        smiles (str): The input SMILES.
        fps_dict (dict): Dictionary of fingerprint classes.

    Returns:
        dict: Dictionary with fingerprint names as keys and fingerprint data as values.
    """
    fp_data = {}
    for fp_name, fp_class in fps_dict.items():
        try:
            fp_array = fp_class.generate_fps(smis=[smiles]).flatten()
            fp_data[fp_name] = ','.join(map(str, fp_array))
        except Exception:
            fp_data[fp_name] = ','.join(['nan'] * fp_class._dimension)  # Handle errors gracefully
    return fp_data

def extract_fingerprints(df):
    """
    Extracts molecular fingerprints and molecular properties for a given DataFrame.

    Args:
        df (pd.DataFrame): Input DataFrame containing a "SMILES" column.

    Returns:
        pd.DataFrame: Updated DataFrame with fingerprint features and molecular properties.
    """

    # Ensure the 'SMILES' column exists
    if "SMILES" not in df.columns:
        raise ValueError("Input DataFrame must contain a 'SMILES' column")

    # Define fingerprint classes
    fingerprint_classes = {
        'ECFP4': HitGenECFP4(),
        'ECFP6': HitGenECFP6(),
        'FCFP4': HitGenFCFP4(),
        'FCFP6': HitGenFCFP6(),
        'MACCS': HitGenMACCS(),
        'RDK': HitGenRDK(),
        'AVALON': HitGenAvalon(),
        'TOPTOR': HitGenTopTor(),
        'ATOMPAIR': HitGenAtomPair()
    }

    # Compute fingerprints and molecular properties
    fingerprint_data = []
    molecular_props = []

    for smiles in df["SMILES"]:
        fps = generate_fingerprints(smiles, fingerprint_classes)
        fingerprint_data.append(fps)
        mw, alogp = compute_molecular_properties(smiles)
        molecular_props.append({"MW": mw, "ALOGP": alogp})

    # Convert lists to DataFrames
    fingerprint_df = pd.DataFrame(fingerprint_data)
    molecular_props_df = pd.DataFrame(molecular_props)
 
    # Concatenate the original DataFrame with fingerprints and molecular properties
    df = pd.concat([df, molecular_props_df, fingerprint_df], axis=1)

    return df
