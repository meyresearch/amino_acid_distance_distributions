"""
Exclusively to be run alongside the embarrassing_parallel.sh script!
Opens each chunked ID file and retrieves PDB files from RCSB to ../../pdbs/.
"""
import numpy as np
from Bio.PDB.PDBList import PDBList
import os
import pandas as pd


def concatenate_id_files() -> np.array:
    """
    Opens id_file_N.txt as a DF and concatenates them all into one DF
    @return: pdb_array
    """
    frames = []
    for filename in os.listdir("."):
        if filename.endswith(".txt"):
            temp_df = pd.read_csv(filename, header=None)
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
                              pdir="../../pdb_files/")
    counter += 1

