'''
Run through PDB IDs to select excluded IDs.

'''
import sys
sys.path.insert(1, '/Users/jasminguven/Documents/GitHub/sequence_distance_distribution/src/')
import PCN
import numpy as np
import pandas as pd 
import os

def concat_ID_files():
    '''
    Opens the chunked PDB id files as a DataFrame. Returns an array with of all the IDs.
    
    Parameters
    ----------
    None

    Return
    ------
    pdb_id_array: numpy array
        array containing all the PDB IDs
    '''
    frames = []
    for filename in os.listdir('/Users/jasminguven/Documents/GitHub/sequence_distance_distribution/data/ids/'):
        if filename.endswith('.txt'):
            temp_df = pd.read_csv('/Users/jasminguven/Documents/GitHub/sequence_distance_distribution/data/ids/'+str(filename), header = None)
            frames.append(temp_df)

    concatenated_DF = pd.concat(frames).reset_index()
    concatenated_DF = concatenated_DF.drop(columns = 'index')
    PDB_ID_array = concatenated_DF[0].to_numpy()
    
    return PDB_ID_array


PDB_IDs = concat_ID_files()
number_PDBs = len(PDB_IDs)
counter = 0
exclude_counter = 0
range_100_PDB = []
range_200_PDB = []
range_300_PDB = []
exclude_ids = []

for PDB_ID in PDB_IDs:
    print('--------------------------------------------------------------------------------------------------')
    print(f'At entry {counter+1} out of {number_PDBs+1} entries.')
    print(f'The PDB ID is {PDB_ID}')
    # This will be replaced with untarring of pdb files
    #-------------------------------------------------------------------
    try:
        # PDB_tarfile = tarfile.open(f'../data/PDBs.tar.gz', 'r:gz')
        PDB_ID_lower = str(PDB_ID).lower()
        # PDB_file = PDB_tarfile.extractfile(f'../data/PDBs/pdb{PDB_ID_lower}.ent')
        PDB_file = f'/Volumes/Seagate_Extension_Plus/PDBs/pdb{PDB_ID_lower}.ent'
        PDB_check = open(PDB_file, 'r')
        print('Successfully opened file.')
    except:
        print(f'Failed at PDB ID: {PDB_ID}. Could not find file \'pdb{PDB_ID_lower}.ent.\'')
        continue # Jump back to the top of the loop and try another PDB
    #-------------------------------------------------------------------

    if os.path.isfile(PDB_file):
        print('Creating PCN object.')
        # Create the PCN instance
        protein_contact_network = PCN.PCN(PDB_file)
        alpha_carbons = protein_contact_network.get_alpha_carbons()
        chain_length = protein_contact_network.get_chain_length(alpha_carbons)
        # Split into chain length ranges
        if chain_length in range(85,116):
            print('This PDB is in the length range 85-115. \nGetting link lengths.')
            range_100_PDB.append(PDB_ID)
        elif chain_length in range(185,216):
            print('This PDB is in the length range 185-215. \nGetting link lengths.')
            range_200_PDB.append(PDB_ID)    
        elif chain_length in range(285,316):
            print('This PDB is in the length range 285-315. \nGetting link lengths.')
            range_300_PDB.append(PDB_ID)
        else:
            print('This PDB is outside of the chosen ranges.')
            exclude_counter += 1
            exclude_ids.append(PDB_ID)

        counter += 1
        
    else:
        continue

print(f'{exclude_counter} out of {number_PDBs} PDB IDs were excluded because they were outside the range.')

print(f'100s: {len(range_100_PDB)}\n200s: {len(range_200_PDB)}\n300s: {len(range_300_PDB)}')
exclude_ids_file = 'exclude'
np.save(exclude_ids_file, exclude_ids)

