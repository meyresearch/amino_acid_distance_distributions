from Bio.PDB.PDBList import PDBList
import os
import pandas as pd


pdblist = PDBList()

def concat_ID_files():
    ''' Opens id_file_N.txt as a DF and concats them all into one DF '''
    frames = []
    for filename in os.listdir('../data/ids/'):
        if filename.endswith(".txt"):
            temp_df = pd.read_csv('../data/ids/'+str(filename), header=None)
            frames.append(temp_df)

    DF = pd.concat(frames).reset_index()
    DF = DF.drop(columns='index')
    PDB_list = DF[0].to_numpy()
    return PDB_list

PDBs = concat_ID_files()
counter = 1
for pdb in PDBs:
    print(f'At entry {counter}')
    pdblist.retrieve_pdb_file(pdb_code = pdb, 
                              file_format = "pdb", 
                              pdir = '/Volumes/Seagate_Extension_Plus/PDBs/')
    counter += 1

# pdb = '1AKI'
# pdblist.retrieve_pdb_file(pdb_code = pdb, 
#                           file_format = "bundle", 
#                           pdir = '../data/bioretrieve_PDBs/')