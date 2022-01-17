"""
Functions to be used in main.py. 

"""
import os
import pandas as pd
import Bio.PDB.PDBList as biopy
import PCN
import traceback
import tarfile
import numpy as np

import warnings
warnings.filterwarnings("ignore")

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

def is_int(input_string, input_var): # NOT WORKING
    """
    Check if user input is an integer.

    Parameters
    ----------
    input_string: str
        string to be printed out to get user input
    input_var: int
        variable to hold user input, should be integer otherwise ask for new input
    
    Return
    ------ 
    integer_input: int
        user input converted to an integer
    
    """
    while True:

        try:
            integer_input = int(input_var)
        
        except ValueError:

            print('Did not understand the command. Try again.')
            continue
        
        else:
            break
    
    return integer_input


def is_file(filename):
    """
    Check if given filename exists. 
    
    Parameters
    ----------
    filename: str
        filename given by user

    Return
    ------
    checked_filename: str
        correct filename

    """
    is_file = os.path.isfile(str(filename))
    while is_file == False:
        checked_filename = input(f'Invalid filename. Please enter again: ')
        is_file = os.path.isfile(str(checked_filename))
    return checked_filename

def get_uniprot_IDs(sequence_length):
    """
    Take the sequence length and return an array of Uniprot IDs in that range.
    Parameters
    ----------
    sequence_length: str
        sequence length of PDBs, either 100, 200 or 300
    
    Return
    ------
    uniprot_IDs: numpy.array
        array containing all the uniprot IDs of the given sequence length
    """
    uniprot_IDs_df = pd.read_excel(f'../data/alphafold_data/uniprot_{sequence_length}s.xlsx')
    uniprot_IDs = uniprot_IDs_df['Entry'].to_numpy()

    return uniprot_IDs

def get_alphafold_PDBs(uniprot_IDs):
    """
    Parameters
    ----------
    uniprot_IDs: numpy.array
        uniprot IDs of proteins in a given sequence length range
    
    Return
    ------

    """
    sequence_length = '100'
    uniprot_IDs = get_uniprot_IDs(sequence_length)
    

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

def bootstrap():
    """
    Take the link lengths and bootstrap the data to get means and 33% confidence intervals.

    Saves the mean and lower and upper bound C.L. of the bootstrapped data and the PDB ids belonging to
    that length range as a csv file.

    Parameters
    ----------
    None

    Return
    ------
    None
    """
    input_file = input('Please enter the input file including the full path to the file: ')
  #  link_length_file = is_file(link_length_file)

    # Get the link length file and concatenate arrays in the file into one
    link_lengths_file = np.load(input_file, allow_pickle = True)
    number_of_PBDs_in_range = link_lengths_file.shape[0]
    link_lengths = np.concatenate(link_lengths_file)
    number_of_link_lengths = link_lengths.shape[0]
    max_link_length = max(link_lengths)

    bootstrap_sample_size = 10
    histogram_list = []
    edges_list = []
    means = []

    rng = np.random.default_rng(seed = 1234)


    bootstrapped_samples = {}

    if len(link_lengths) > 0:
        print('Bootstrapping...')
        for n in range(bootstrap_sample_size):
            print(f'At step {n+1}')
            # Randomly choose a list of observations from the sample            
            random_choice = rng.choice(np.arange(number_of_link_lengths), replace = True, size = number_of_link_lengths)
            bootstrap_sample = [link_lengths[i] for i in random_choice]
            # Get number of counts in each link length bin
            values, counts = np.unique(bootstrap_sample, return_counts = True)    
            # Put all the bootstrapped samples into a dict
            for value, count in zip(values, counts):
                # Try to add a count to a link length bin
                try: 
                    bootstrapped_samples[value].append(count)
                except KeyError:
                    # Create the link length bin in the dict
                    bootstrapped_samples[value] = []
                    # Append the count to the newly created bin
                    bootstrapped_samples[value].append(count)

        # Create a DataFrame from the dictionary
        bootstrap_df_T = pd.DataFrame.from_dict(bootstrapped_samples, orient = 'index')
        bootstrap_df = bootstrap_df_T.transpose() 
        bootstrap_df.to_csv(f'../data/data_for_plotting/bs_300_df.csv', index = False)

        # Melt the columns into the link length bin and their corresponding values
        melted_bootstrap_df = bootstrap_df.melt()

        # Normalise by the number of link lengths and get the mean and upper and lower bounds of the CL
        #melted_bootstrap_df.loc[:, 'value'] /= number_of_PBDs_in_range
        
        # Group all the values into the corresponding bins, i.e. by the default 'variable' column
        bootstrap_df_stats = melted_bootstrap_df.groupby('variable', as_index = False).agg(mean = ('value', np.mean),
                                                                                           lower_bound = ('value', lambda val: np.quantile(val, q = 0.05)),
                                                                                           upper_bound = ('value', lambda val: np.quantile(val, q = 0.95)))
        print(bootstrap_df_stats.head())
        output_file = input('Please enter filename for output file: ')
        print(f'The file {output_file} will be saved in \'../data/exclude_subgraphs/\'')
        bootstrap_df_stats.to_csv(f'../data/exclude_subgraphs/{output_file}', index = False)

        #     # Sample with replacement
        #     bootstrap_sample = [item for items in bootstrap_sample for item in items]
        #     histogram, edges = np.histogram(bootstrap_sample, bins = range(0, max_link_length), density = True)
        #     # Normalise over the number of PDBs in the set
        #     histogram = histogram / len(length_range)
        #     histogram_list.append(histogram)
        #     edges_list.append(edges)
      
        # histogram_array = np.array(histogram_list)
        # means = histogram_array.mean(axis = 0)
        # edges_array = np.array(edges_list)
      
        # mean_save_file = f'../data/data_for_plotting/means_{length}.npy'
        # edges_save_file = f'../data/data_for_plotting/edges_{length}.npy'
        # id_save_file = f'../data/data_for_plotting/ids_{length}.npy'
     
        # np.save(mean_save_file, means)
        # np.save(edges_save_file, edges_array)
        # np.save(id_save_file, length_range)

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
            print('--------------------------------------------------------------------------------------------------')
            print(f'At entry {counter+1} out of {number_PDBs+1} entries.')
            if PDB_IDs[i] not in exclude_IDs:
                
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
                    #print('Creating PCN object.')
                    # Create the PCN instance
                    protein_contact_network = PCN.PCN(PDB_file)
                    alpha_carbons = protein_contact_network.get_alpha_carbons()
                    chain_length = protein_contact_network.get_chain_length(alpha_carbons)
                    # Split into chain length ranges
                    if chain_length in range(85,116):
                        print('This PDB is in the length range 85-115. \nGetting link lengths.')
                        range_100_PDB.append(PDB_IDs[i])
                        link_lengths = protein_contact_network.get_link_lengths(alpha_carbons, PDB_IDs[i])
                        link_lengths_100.append(link_lengths)
                        
                    if chain_length in range(185,216):
                        print('This PDB is in the length range 185-215. \nGetting link lengths.')
                        range_200_PDB.append(PDB_IDs[i])
                        link_lengths = protein_contact_network.get_link_lengths(alpha_carbons, PDB_IDs[i])
                        link_lengths_200.append(link_lengths)    
                    if chain_length in range(285,316):
                        print('This PDB is in the length range 285-315. \nGetting link lengths.')
                        range_300_PDB.append(PDB_IDs[i])
                        link_lengths = protein_contact_network.get_link_lengths(alpha_carbons, PDB_IDs[i])
                        link_lengths_300.append(link_lengths)
                
                else:
                    continue
           
            else:
                print('This PDB is outside of the chosen ranges.')

            counter += 1

        print('--------------------------------------------------------------------------------------------------')
        print(f'100s: {len(range_100_PDB)}\n200s: {len(range_200_PDB)}\n300s: {len(range_300_PDB)}')

        np.save(f'../data/exclude_subgraphs/lls_100.npy', link_lengths_100)
        np.save(f'../data/exclude_subgraphs/lls_200.npy', link_lengths_200)
        np.save(f'../data/exclude_subgraphs/lls_300.npy', link_lengths_300)
        np.save(f'../data/exclude_subgraphs/ids_100.npy', range_100_PDB)
        np.save(f'../data/exclude_subgraphs/ids_200.npy', range_200_PDB)
        np.save(f'../data/exclude_subgraphs/ids_300.npy', range_300_PDB)
        # print('Bootstrapping...')
        # # Bootstrap the data
        # bootstrap('100', link_lengths_100, range_100_PDB)
        # bootstrap('200', link_lengths_200, range_200_PDB)
        # bootstrap('300', link_lengths_300, range_300_PDB)
        print('Done.')
        print('--------------------------------------------------------------------------------------------------')


def user_options():
    """
    Print out options for user when main() is run.

    Parameters
    ----------
    None

    Return
    ------
    choice: int
        user input option to either 1 get link lengths or 2 bootstrap
    """

    print('--------------------------------------------------------------------------------------------------')
    print('What would you like to do?')
    print('1) Get link lengths\n2) Bootstrap link lengths')
    choice = input('Please choose: ')
    # choice = is_int('Please choose: ', choice)

    return choice

def run():
    """
    Run the options again as many times as user wants. 

    Parameters
    ----------
    None

    Return
    ------
    None

    """
    # repeat = True

    # while repeat == True:
        
    user_choice = user_options()
    user_choice = int(user_choice)
    # User choice: get link lengths
    if user_choice == 1: 
        PDB_to_PCN('log.txt')
    
    # User choice: bootstrap
    elif user_choice == 2: 
        bootstrap()

    # print('Do you want to do another run or quit?')
    # repeat_option = input('1. Run again \n2. Quit')
    # repeat_option = is_int(repeat_option)

        # if repeat_option == 2:
        #     repeat = False

def main():

    run()

    











