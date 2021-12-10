"""
Functions to be used in main.py. 

"""
import os
import pandas as pd
import Bio.PDB.PDBList as biopy
import PCN
import traceback
import tarfile
import random
import numpy as np

def concat_ID_files():
    """
    Opens the chunked PDB id files as a DataFrame. Returns an array with of all the IDs.
    
    Parameters
    ----------
    None

    Return
    ------
    PDB_ID_array: numpy array
        array containing all the PDB IDs
    """
    frames = []
    for filename in os.listdir('../data/ids/'):
        if filename.endswith('.txt'):
            temp_df = pd.read_csv('../data/ids/'+str(filename), header = None)
            frames.append(temp_df)

    concatenated_DF = pd.concat(frames).reset_index()
    concatenated_DF = concatenated_DF.drop(columns = 'index')
    PDB_ID_array = concatenated_DF[0].to_numpy()
    
    return PDB_ID_array

def bioretrieve_PDBs(PDB_ID_array, save_directory):
    """
    Take the array of PDB IDs and use BioPython to retrieve PDB files as .ent files.

    Parameters
    ----------
    PDB_ID_array: numpy array
        array containing all the PDB IDs
    save_directory: string 
        path to the directory where the file will be saved

    Return
    ------
    None
    """
    pdblist = biopy.PDBList()
    PDB_IDs = concat_ID_files()
    counter = 1
    for PDB_ID in PDB_IDs:
        print(f'At entry {counter} out of {len(PDB_IDs)} entries.')
        pdblist.retrieve_pdb_file(pdb_code = PDB_ID, file_format = 'pdb', pdir = save_directory)
        counter += 1

def bootstrap(length, link_lengths, length_range):
    """
    Take the link lengths and bootstrap the data. 

    Saves the mean and median of the bootstrapped data and the PDB ids belonging to
    that length range. Files are saved in .npy format. 

    Parameters
    ----------
    length: string
        defines the length range as '100', '200' or '300'

    link_lengths: list
        list containing lengths of links

    length_range: list
        list containing the PDB ids that fall in the length range given

    Return
    ------
    None
    """
    if len(link_lengths) > 0:
        sample_size = 1000
        histogram_list = []
        means = []

        for n in range(sample_size):
            # Randomly choose a list of observations from the sample
            bootstrap_sample = [random.choice(link_lengths) for _ in link_lengths]
            # Sample with replacement
            bootstrap_sample = [item for items in bootstrap_sample for item in items]
            histogram, _ = np.histogram(bootstrap_sample, bins = range(3, 1000), density = True)
            # Normalise over the number of PDBs in the set
            histogram = histogram / len(length_range)
            histogram_list.append(histogram)
        histogram_array = np.array(histogram_list)
        means = histogram_array.mean(axis = 0)
        medians = np.median(histogram_array, axis = 0)

        mean_save_file = f'data_for_plotting/means_{length}.npy'
        median_save_file = f'data_for_plotting/medians_{length}.npy'
        id_save_file = f'data_for_plotting/ids_{length}.npy'
        np.save(mean_save_file, means)
        np.save(median_save_file, medians)
        np.save(id_save_file, length_range)

def PDB_to_PCN(log_file):
    """
    Creates an instance of the PCN class from the PDB files. 

    Bootstraps the data, which in turn saves the data for plotting.

    Parameters
    ----------
    log_file: string
        name (and path) to the log.txt file where error messages are written
    Return
    ------
    None
    """
    # Get PDB IDs
    PDB_IDs = concat_ID_files()
    number_PDBs = len(PDB_IDs)

    range_100_PDB = []
    range_200_PDB = []
    range_300_PDB = []
    link_lengths_100 = []
    link_lengths_200 = []
    link_lengths_300 = []

    counter = 0
    exclude_counter = 0

    exclude_IDs = np.load('../data/exclude.npy')

    # Write errors into a log file
    with open(log_file, 'w') as log:     

        for i in range(number_PDBs):

            if PDB_IDs[i] not in exclude_IDs:
                print('--------------------------------------------------------------------------------------------------')
                print(f'At entry {counter+1} out of {number_PDBs+1} entries.')
                print(f'The PDB ID is {PDB_IDs[i]}')

                try:
                    # PDB_tarfile = tarfile.open(f'../data/PDBs.tar.gz', 'r:gz')
                    PDB_IDs_lower = str(PDB_IDs[i]).lower()
                    # PDB_file = PDB_tarfile.extractfile(f'../data/PDBs/pdb{PDB_IDs_lower}.ent')
                    PDB_file = f'/Volumes/Seagate_Extension_Plus/PDBs/pdb{PDB_IDs_lower}.ent'
                    PDB_check = open(PDB_file, 'r')
                    print('Successfully opened file.')
                except:
                    traceback.print_exc(file=log) # Log the exception
                    continue # Jump back to the top of the loop and try another PDB


                if os.path.isfile(PDB_file):
                    print('Creating PCN object.')
                    # Create the PCN instance
                    protein_contact_network = PCN.PCN(PDB_file)
                    C_alphas = protein_contact_network.get_C_alphas()
                    chain_length = protein_contact_network.get_chain_length(C_alphas)
                    # Split into chain length ranges
                    if chain_length in range(85,116):
                        print('This PDB is in the length range 85-115. \nGetting link lengths.')
                        range_100_PDB.append(PDB_IDs[i])
                        link_lengths = protein_contact_network.get_link_lengths(C_alphas)
                        link_lengths_100.append(link_lengths)
                    elif chain_length in range(185,216):
                        print('This PDB is in the length range 185-215. \nGetting link lengths.')
                        range_200_PDB.append(PDB_IDs[i])
                        link_lengths = protein_contact_network.get_link_lengths(C_alphas)
                        link_lengths_200.append(link_lengths)    
                    elif chain_length in range(285,316):
                        print('This PDB is in the length range 285-315. \nGetting link lengths.')
                        range_300_PDB.append(PDB_IDs[i])
                        link_lengths = protein_contact_network.get_link_lengths(C_alphas)
                        link_lengths_300.append(link_lengths)
                
                else:
                    continue
           
            else:
                print('This PDB is outside of the chosen ranges.')

            counter += 1

        print('--------------------------------------------------------------------------------------------------')
        print(f'100s: {len(range_100_PDB)}\n200s: {len(range_200_PDB)}\n300s: {len(range_100_PDB)}')
        print('Bootstrapping...')
        # Bootstrap the data
        bootstrap('100', link_lengths_100, range_100_PDB)
        bootstrap('200', link_lengths_200, range_200_PDB)
        bootstrap('300', link_lengths_300, range_300_PDB)
        print('Done.')
        print('--------------------------------------------------------------------------------------------------')


def main():
    PDB_to_PCN('log.txt')











