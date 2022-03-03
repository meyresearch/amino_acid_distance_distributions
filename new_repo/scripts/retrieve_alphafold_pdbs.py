"""
Get UniProt ids and loop through them to save AlphaFold PDBs in ../data/alphafold_data/pdb_files/
"""
import urllib.request
import numpy as np
import pandas as pd


def get_uniprot_ids(length_range: str) -> np.array:
    """
    Creates a numpy array containing uniprot ids for given length_range
    @param length_range:
    @return: uniprot_ids
    """
    uniprot_ids_df = pd.read_excel(f"../data/alphafold/uniprot_ids/uniprot_{length_range}s.xlsx")
    uniprot_ids = uniprot_ids_df["Entry"].to_numpy()
    return uniprot_ids


sequence_lengths = ["100", "200", "300"]

for sequence_length in sequence_lengths:
    counter = 0

    uniprot_IDs = get_uniprot_ids(sequence_length)
    for uniprot_ID in uniprot_IDs:
        print(f"At entry {counter}/{len(uniprot_IDs)}")
        download = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_ID}-F1-model_v2.pdb"

        try:
            file_name = f"../data/alphafold/pdb_files/AF-{uniprot_ID}-F1-model_v2.pdb"
            urllib.request.urlretrieve(download, file_name)
        except ValueError:
            print("Error downloading file.")

        counter += 1
         
