"""
functions.py

"""
import os
import pandas as pd
import Bio.PDB.PDBList as BioPy
import PCN
import traceback
import numpy as np
import argparse
from warnings import filterwarnings

filterwarnings("ignore")


def concat_id_files() -> np.array:
    """
    @return: concatenated array of pdb ids
    """
    frames = []
    for filename in os.listdir("../data/ids/"):
        if filename.endswith(".txt"):
            temp_df = pd.read_csv("../data/ids/" + str(filename), header=None)
            frames.append(temp_df)

    concatenated_df = pd.concat(frames).reset_index()
    concatenated_df = concatenated_df.drop(columns="index")
    pdb_id_array = concatenated_df[0].to_numpy()

    return pdb_id_array




def is_file(filename):
    is_file = os.path.isfile(str(filename))
    while is_file == False:
        checked_filename = input(f"Invalid filename. Please enter again: ")
        is_file = os.path.isfile(str(checked_filename))
    return checked_filename


def get_uniprot_ids(sequence_length: str) -> np.array:
    """

    @param sequence_length: 
    @return: 
    """
    uniprot_ids_df = pd.read_excel(f"../data/alphafold_data/uniprot_{sequence_length}s.xlsx")
    uniprot_ids = uniprot_ids_df["Entry"].to_numpy()

    return uniprot_ids


def get_alphafold_pdbs(uniprot_ids):
    sequence_length = "100"
    uniprot_ids = get_uniprot_ids(sequence_length)


def bioretrieve_pdbs(pdb_id_array: np.array, save_directory: str) -> None:
    """

    @param pdb_id_array: 
    @param save_directory: 
    """
    pdblist = BioPy.PDBList()
    pdb_ids = concat_id_files()
    counter = 1
    for pdb_id in pdb_ids:
        print(f"At entry {counter} out of {len(pdb_ids)} entries.")
        pdblist.retrieve_pdb_file(pdb_code=pdb_id, file_format="pdb", pdir=save_directory)
        counter += 1


def bootstrap(inputfile: str, sample_replacement: bool, command_line_link_length: str) -> None:
    """
    @param command_line_link_length:
    @param inputfile: name and full path to input file
    @param sample_replacement: sample with or without replacement
    @return: None
    """

    # Get the link length file and concatenate arrays in the file into one
    link_lengths_file = np.load(inputfile, allow_pickle=True)
    # number_of_pdbs_in_range = link_lengths_file.shape[0]
    link_lengths = np.concatenate(link_lengths_file)
    number_of_link_lengths = link_lengths.shape[0]

    bootstrap_sample_size = 1000
    random_number_generator = np.random.default_rng(seed=1234)
    bootstrapped_samples = {}
    if len(link_lengths) > 0:
        print("Bootstrapping...")
        for n in range(bootstrap_sample_size):
            print(f"At step {n + 1}")
            # Randomly choose a list of observations from the sample            
            random_choice = random_number_generator.choice(np.arange(number_of_link_lengths),
                                                           replace=sample_replacement, size=number_of_link_lengths)
            bootstrap_sample = [link_lengths[i] for i in random_choice]
            # Get number of counts in each link length bin
            values, counts = np.unique(bootstrap_sample, return_counts=True)
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
        bootstrap_df_t = pd.DataFrame.from_dict(bootstrapped_samples, orient="index")
        bootstrap_df = bootstrap_df_t.transpose()
        bootstrap_df_nona = bootstrap_df.fillna(0)
        bootstrap_df_nona.to_csv(f"../data/rcsb_data/bs_{command_line_link_length}_df_no_stats.csv", index=False)

        # Melt the columns into the link length bin and their corresponding values
        melted_bootstrap_df = bootstrap_df_nona.melt()

        # Group all the values into the corresponding bins, i.e. by the default "variable" column
        bootstrap_df_stats = melted_bootstrap_df.groupby("variable", as_index=False).agg(mean=("value", np.mean),
                                                                                         lower_bound=("value", lambda
                                                                                             val: np.quantile(val,
                                                                                                              q=0.05)),
                                                                                         upper_bound=("value", lambda
                                                                                             val: np.quantile(val,
                                                                                                              q=0.95)))
        print(bootstrap_df_stats.head())
        bootstrap_df_stats.to_csv(f"../data/rcsb_data/bs_{command_line_link_length}_stats.csv", index=False)
        # if outputfile is None:
        #     bootstrap_df_stats.to_csv(f"bootstrapped_data", index=False)
        # else:


def pdb_to_pcn(log_file: str) -> None:
    """

    @param log_file:
    @return:
    """
    # Get PDB IDs
    pdb_ids = concat_id_files()
    # pdb_ids = np.load("../data/pdb_data_NaNs/ids_200.npy")
    number_PDBs = len(pdb_ids)

    range_100_PDB = []
    range_200_PDB = []
    range_300_PDB = []
    link_lengths_100 = []
    link_lengths_200 = []
    link_lengths_300 = []

    counter = 0
    exclude_counter = 0

    exclude_IDs = np.load("../data/exclude.npy")

    # Write errors into a log file
    with open(log_file, "w") as log:

        for i in range(number_PDBs):
            print("--------------------------------------------------------------------------------------------------")
            print(f"At entry {counter + 1} out of {number_PDBs + 1} entries.")
            if pdb_ids[i] not in exclude_IDs:

                print(f"The PDB ID is {pdb_ids[i]}")

                try:
                    # PDB_tarfile = tarfile.open(f"../data/PDBs.tar.gz", "r:gz")
                    pdb_ids_lower = str(pdb_ids[i]).lower()
                    # PDB_file = PDB_tarfile.extractfile(f"../data/PDBs/pdb{pdb_ids_lower}.ent")
                    PDB_file = f"/Volumes/Seagate_Extension_Plus/PDBs/pdb{pdb_ids_lower}.ent" # replace with get_pdbs.py
                    PDB_check = open(PDB_file, "r")
                    print("Successfully opened file.")
                except:
                    traceback.print_exc(file=log)  # Log the exception
                    continue  # Jump back to the top of the loop and try another PDB

                if os.path.isfile(PDB_file):
                    # print("Creating PCN object.")
                    # Create the PCN instance
                    protein_contact_network = PCN.PCN(PDB_file)
                    alpha_carbons = protein_contact_network.get_alpha_carbons()
                    chain_length = protein_contact_network.get_chain_length(alpha_carbons)
                    # Split into chain length ranges
                    if chain_length in range(85, 116):
                        print("This PDB is in the length range 85-115. \nGetting link lengths.")
                        range_100_PDB.append(pdb_ids[i])
                        link_lengths = protein_contact_network.get_link_lengths(alpha_carbons, pdb_ids[i])
                        link_lengths_100.append(link_lengths)
                        # link_lengths_100 = np.array(link_lengths_100)
                    if chain_length in range(185, 216):
                        print("This PDB is in the length range 185-215. \nGetting link lengths.")
                        range_200_PDB.append(pdb_ids[i])
                        link_lengths = protein_contact_network.get_link_lengths(alpha_carbons, pdb_ids[i])
                        link_lengths_200.append(link_lengths)
                        # link_lengths_200 = np.array(link_lengths_200)
                    if chain_length in range(285, 316):
                        print("This PDB is in the length range 285-315. \nGetting link lengths.")
                        range_300_PDB.append(pdb_ids[i])
                        link_lengths = protein_contact_network.get_link_lengths(alpha_carbons, pdb_ids[i])
                        link_lengths_300.append(link_lengths)
                        # link_lengths_300 = np.array(link_lengths_300)

                else:
                    continue

                counter += 1

            else:
                print("This PDB is outside of the chosen ranges.")
                counter += 1



        print("--------------------------------------------------------------------------------------------------")
        print(f"100s: {len(range_100_PDB)}\n200s: {len(range_200_PDB)}\n300s: {len(range_300_PDB)}")
        link_lengths_100 = np.array(link_lengths_100)
        link_lengths_200 = np.array(link_lengths_200)
        link_lengths_300 = np.array(link_lengths_300)

        np.save(f"../data/rcsb_data/lls_100.npy", link_lengths_100)
        np.save(f"../data/rcsb_data/lls_200.npy", link_lengths_200)
        np.save(f"../data/rcsb_data/lls_300.npy", link_lengths_300)
        np.save(f"../data/rcsb_data/ids_100.npy", range_100_PDB)
        np.save(f"../data/rcsb_data/ids_200.npy", range_200_PDB)
        np.save(f"../data/rcsb_data/ids_300.npy", range_300_PDB)
        # print("Bootstrapping...")
        # # Bootstrap the data
        # bootstrap("100", link_lengths_100, range_100_PDB)
        # bootstrap("200", link_lengths_200, range_200_PDB)
        # bootstrap("300", link_lengths_300, range_300_PDB)
        print("Done.")
        print("--------------------------------------------------------------------------------------------------")


def run(arguments: argparse.Namespace) -> None:
    """
    @param arguments: command line arguments given by user
    @return: None
    """
    algorithm = arguments.algorithm
    inputfile = arguments.inputfile
    # outputfile = arguments.outputfile
    link_length = arguments.link_length
    if algorithm == "LL":
        pdb_to_pcn("log.txt")
    elif algorithm == "BS":
        bootstrap(inputfile, True, link_length)
    elif algorithm == "C":
        bootstrap(inputfile, False, link_length)


def command_line_arguments() -> argparse.Namespace:
    """
    @return: command line arguments passed by user
    """
    parser = argparse.ArgumentParser(description="Link lengths code.")
    parser.add_argument("algorithm",
                        type=str,
                        choices=["LL", "BS", "C"],
                        help="get link lengths (LL), bootstrap (BS) or chunk (C)")
    parser.add_argument("link_length", type=str, choices=["100", "200", "300"], help="link length range")
    parser.add_argument("-i",
                        dest="inputfile",
                        type=str,
                        help="full path and name of input file")
    parser.add_argument("-o",
                        dest="outputfile",
                        type=str,
                        help="name of output file")

    arguments: argparse.Namespace = parser.parse_args()

    if arguments.algorithm == "BS" and arguments.inputfile is None:
        parser.error("BS requires -i")
    elif arguments.algorithm == "C" and arguments.inputfile is None:
        parser.error("C requires -i")

    return arguments


def main():
    arguments = command_line_arguments()
    run(arguments)


if __name__ == "__main__":
    main()
