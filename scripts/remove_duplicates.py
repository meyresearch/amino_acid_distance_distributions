"""
Remove duplicates from AlphaFold 2 PDB files
"""
import glob
import pandas as pd
import numpy as np
import argparse


def get_protein_names(filenames: np.ndarray) -> np.ndarray:
    """
    Loop through the filenames in the secondary structure dataframe and get the protein names
    @param filenames: array of full paths to pdb files in secondary structure dataframes
    @return: array of protein names
    """
    protein_names = []
    counter = 1
    for file in filenames:
        try:
            with open(file) as pdb:
                print(f"Progress: {counter}/{len(filenames)}") 
                header = [next(pdb) for i in range(2)]
                protein_name = header[-1].strip().split('FOR')[-1]
                protein_names.append(protein_name)
                counter += 1
        except TypeError as e:
            print(f"NaN in filename, excluding from list: {e}")
    return np.asarray(protein_names)
            

def filter_unique_secondary_structures(low_confidence: bool) -> None:
    """
    Open secondary structures csv and find unique protein structures and save as a csv.
    @param low_confidence: True only use low confidence structures, False use high confidence structures
    @return: None
    """
    chain_lengths = ["100", "200", "300"]
    for length in chain_lengths:
        if low_confidence:
            secondary_structures_df = pd.read_csv(f"../data/alphafold/low_secondary_structures_{length}.csv")
        elif not low_confidence:
            secondary_structures_df = pd.read_csv(f"../data/alphafold/secondary_structures_{length}.csv")

        filenames = secondary_structures_df["filename"]
        secondary_structures_df["protein_name"] = get_protein_names(filenames)
        unique_structures_df = secondary_structures_df.drop_duplicates(subset="protein_name")
        if low_confidence:
            unique_structures_df.to_csv(f"../data/alphafold/unique_low_secondary_structures_{length}.csv")
        elif not low_confidence:
            unique_structures_df.to_csv(f"../data/alphafold/unique_secondary_structures_{length}.csv")


def commandline_arguments() -> argparse.Namespace:
    """
    Parser for commandline arguments
    @return: commandline arguments from user
    """
    parser = argparse.ArgumentParser(description="filter AlphaFold 2 structures based on per-residue confidence scores")
    parser.add_argument("-l",
                        "--low",
                        help="get low-confidence AlphaFold 2 structures, default behaviour is to only get structures with confidence above 90",
                        action="store_true")
    return parser.parse_args()


def main():
    arguments = commandline_arguments()
    low_confidence = arguments.low
    filter_unique_secondary_structures(low_confidence)

    
if __name__ == "__main__":
    main()
    