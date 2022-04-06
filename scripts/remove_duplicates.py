"""
Remove duplicates from AlphaFold 2 PDB files
"""
import glob
import pandas as pd
import numpy as np


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
            

def filter_unique_secondary_structures() -> None:
    """
    Open secondary structures csv and find unique protein structures and save as a csv.
    @return: None
    """
    chain_lengths = ["100","200","300"]
    pdb_files = glob.glob("../data/alphafold/pdb_files/*.pdb")
    for length in chain_lengths:
        secondary_structures_df = pd.read_csv(f"../data/alphafold/secondary_structures_{length}.csv")
        filenames = secondary_structures_df["filename"]
        secondary_structures_df["protein_name"] = get_protein_names(filenames)
        unique_structures_df = secondary_structures_df.drop_duplicates(subset="protein_name")
        unique_structures_df.to_csv(f"../data/alphafold/unique_secondary_structures_{length}.csv")

        
def main():
    filter_unique_secondary_structures()

    
if __name__ == "__main__":
    main()
    