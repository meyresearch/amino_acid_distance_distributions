"""
Opens each chunked ID file and retrieves PDB files from RCSB to ../data/rcsb/pdb_files/.
"""
from Bio.PDB.PDBList import PDBList
import os
import pandas as pd


def concatenate_id_files():
    """
    Opens id_file_N.txt as a DF and concatenates them all into one DF
    @return:
    """
    frames = []
    for filename in os.listdir("../data/rcsb/ids/"):
        if filename.endswith(".txt"):
            temp_df = pd.read_csv(f"../data/rcsb/ids/{filename}", header=None)
            frames.append(temp_df)
    dataframe = pd.concat(frames).reset_index()
    dataframe = dataframe.drop(columns='index')
    pdb_array = dataframe[0].to_numpy()
    return pdb_array


pdblist = PDBList()
pdb_ids = concatenate_id_files()
counter = 1
for pdb in pdb_ids:
    print(f"At entry {counter}")
    pdblist.retrieve_pdb_file(pdb_code=pdb,
                              file_format="pdb",
                              pdir="../data/rcsb/pdb_files/")
    counter += 1

