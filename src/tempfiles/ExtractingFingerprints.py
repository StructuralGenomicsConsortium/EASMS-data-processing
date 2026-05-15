import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from fingerprints import HitGenMACCS, HitGenECFP4, HitGenECFP6, HitGenFCFP4, HitGenFCFP6, HitGenRDK, HitGenAvalon, HitGenTopTor, HitGenAtomPair

def compute_molecular_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mw = Descriptors.MolWt(mol)
        alogp = Descriptors.MolLogP(mol)
    else:
        mw = np.nan
        alogp = np.nan
    return mw, alogp

def generate_fingerprints(smiles, fps_dict):
    fp_data = {}
    for fp_name, fp_class in fps_dict.items():
        try:
            fp_array = fp_class.generate_fps(smis=[smiles]).flatten()
            fp_data[fp_name] = ','.join(map(str, fp_array))
        except Exception:
            fp_data[fp_name] = ','.join(['nan'] * fp_class._dimension)
    return fp_data

def process_file(input_file, output_file, fingerprints, nrows=None):
    # Read the file
    df = pd.read_csv(input_file, nrows=nrows)

    '''# Add MW and ALOGP columns
    df[['MW', 'ALOGP']] = df['SMILES (Compounds)'].apply(
        lambda smi: pd.Series(compute_molecular_properties(smi))
    )'''

    # Generate fingerprint columns
    fingerprint_data = []
    for smiles in df['smiles']:
        fps = generate_fingerprints(smiles, fingerprints)
        fingerprint_data.append(fps)

    # Create a DataFrame from fingerprint data
    fingerprint_df = pd.DataFrame(fingerprint_data)

    # Concatenate fingerprint data with the main DataFrame
    df = pd.concat([df, fingerprint_df], axis=1)

    # Save the new DataFrame to a CSV file
    df.to_csv(output_file, index=False)

    print(f"The updated file with fingerprints has been saved as '{output_file}'")

def main():
    # Define fingerprint classes
    # You can also extract the binary versions (look at the fingerprints.py )
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

    nrows = None
    input_file = r"D:\0000-UHN\03-DataAndCodes\AIRCHECK-workflow\SimpleML\Bootcamp\Data\21Feb\ASMS_hits_clustered.csv"
    output_file = r"D:\0000-UHN\03-DataAndCodes\AIRCHECK-workflow\SimpleML\Bootcamp\Data\21Feb\ASMS_hits_clustered_with_fingerprints.csv"
    process_file(input_file, output_file, fingerprint_classes, nrows=nrows)

    '''input_file = "James_hits.csv"
    output_file = "James_hits_fingerprints.csv" 
    process_file(input_file, output_file, fingerprint_classes, nrows=nrows)
    input_file = "ASMS_hits.csv"
    output_file = "ASMS_hits_fingerprints.csv" 
    process_file(input_file, output_file, fingerprint_classes, nrows=nrows)'''
    
    '''
    fingerprint_classes = {
        'ECFP6': HitGenECFP4()
    }
    nrows = None # Use None to process the entire file
    # Example usage
    input_file = "ASMS_460K.csv"
    output_file = "ASMS_460K_with_ECFP6.csv"
    process_file(input_file, output_file, fingerprint_classes, nrows=nrows)
    '''


if __name__ == "__main__":
    main()
