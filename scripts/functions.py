"""
Functions for getting link lengths and bootstrapping.
"""
import argparse
import glob
import os
import traceback
import networkx as nx
import numpy as np
import pandas as pd
import protein_contact_map
import subprocess as sp
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)


def concatenate_rcsb_id_files() -> np.array:
    """
    Opens id_file_N.txt as a DF and concatenates them all into one DF for RCSB PDBs.
    @return: pdb_array
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


def get_uniprot_ids(sequence_length: str) -> np.ndarray:
    """
    Return Uniprot IDs in the given sequence length range.
    @param sequence_length
    @return: np.ndarray of Uniprot IDs
    """
    id_dataframe = pd.read_excel(f"../data/alphafold/uniprot_ids/uniprot_{sequence_length}s.xlsx")
    return id_dataframe["Entry"].to_numpy()


def commandline_arguments() -> argparse.Namespace:
    """
    Parser for commandline arguments
    @return: commandline arguments from user
    """

    parser = argparse.ArgumentParser(description="Get amino acid distance distributions.")
    parser.add_argument("algorithm", type=str, choices=["rcsb", "alpha", "3d-sim"],
                        help="Get amino acid distances from rcsb, [alpha]fold or 3D simulations")
    parser.add_argument("-r", "--range", dest="length_range", type=str, choices=[None, "100", "200", "300"],
                        help="Chain length range")
    parser.add_argument("-p", "--path", dest="path_to_pdbs", type=str,
                        help="Full path to file containing PDB ID csv file")
    cl_arguments = parser.parse_args()
    check_arguments(cl_arguments, parser)
    return cl_arguments


def check_arguments(arguments: argparse.Namespace, argument_parser: argparse.ArgumentParser) -> None:
    """
    Check the requirements for command line arguments
    @param arguments: command line arguments
    @param argument_parser: parser for command line arguments
    @return: None
    """
    if arguments.algorithm == "rcsb" and arguments.path_to_pdbs is None:
        argument_parser.error("rcsb requires --path")
    elif arguments.algorithm == "alpha" and arguments.length_range is None:
        argument_parser.error("alpha requires --range")
    elif arguments.algorithm == "alpha" and arguments.path_to_pdbs is None:
        argument_parser.error("alpha requires --path")
    elif arguments.algorithm == "2d-sim" and arguments.length_range is not None:
        argument_parser.error("2d-sim requires --range")
    elif arguments.algorithm == "2d-sim" and arguments.path_to_pdbs is not None:
        argument_parser.error("2d-sim requires --path=None")
    elif arguments.algorithm == "3d-sim" and arguments.length_range is None:
        argument_parser.error("3d-sim requires --range")
    elif arguments.algorithm == "3d-sim" and arguments.path_to_pdbs is not None:
        argument_parser.error("3d-sim requires --path=None")
    elif arguments.algorithm == "3d-sim" and arguments.inputfile is not None:
        argument_parser.error("3d-sim requires --file=None")
        
        
def pdb_to_adjacency(pdb_file: str) -> tuple:
    """
    Convert given PDB file to an adjacency matrix
    @param pdb_file: PDB file from RCSB or AlphaFold in each length range
    @return: tuple of the adjacency and distance matrices as numpy arrays
    """
    pcm = protein_contact_map.ProteinContactMap(pdb_file)
    alpha_carbons = pcm.get_alpha_carbons
    distance_array = protein_contact_map.get_distance_array(alpha_carbons)
    distance_matrix = protein_contact_map.get_distance_matrix(alpha_carbons, distance_array)
    adjacency_matrix = pcm.get_adjacency_matrix(alpha_carbons, distance_array)
    return distance_matrix, adjacency_matrix


def get_distances(adjacency_matrix: np.ndarray) -> np.ndarray:
    """
    Use adjacency array to get the amino acid distances
    @param adjacency_matrix: adjacency matrix from PDB file
    @return: array of distances in each range
    """
    distances_list = []
    for row_value in range(len(adjacency_matrix)):
        for col_value in range(len(adjacency_matrix)):
            if adjacency_matrix[row_value][col_value] == 1:
                distance = np.abs(col_value - row_value)
                distances_list.append(distance)
    return np.asarray(distances_list)


def get_shadow_distances(path: str):
    """
    Read in .csv file containing paths to pdb files and return amino acid distances for Shadow map
    @param path: full path to the csv file containing paths to pdb files
    @return: ???
    """
    dataframe = pd.read_csv(path)
    pdb_files = dataframe["filename"].tolist()
    pdb_files = pdb_files[:10] # uncomment for debugging
    for pdb in pdb_files:
        


def return_distance_histogram(log_file: str, given_algorithm: str, length_range: str, path_to_csvs: str) -> np.ndarray:
    """
    Compute the amino acid distance distribution for PDB files in given range from adjacency matrix
    and save in a numpy file.
    @param log_file: file to save exceptions in
    @param given_algorithm: alpha or rcsb
    @param length_range: 100, 200 or 300
    @param path_to_csvs: full path to csv files
    @return: None
    """
    dataframe = pd.read_csv(path_to_csvs)
    pdb_files = dataframe["filename"].to_numpy()
    histogram_list = []
    counter = 1
    with open(log_file, "w") as log_file:
        for pdb_file in pdb_files:
            print(f"Progress: {counter}/{len(pdb_files)}")
            if given_algorithm == "alpha":
                clean_pdb_filename = pdb_file.replace("/home/jguven/Projects/sequence_distance_distribution", "..")
            else:
                clean_pdb_filename = pdb_file
            try:
                adjacency_matrix = pdb_to_adjacency(clean_pdb_filename)[1]
                distances = get_distances(adjacency_matrix)
                # bins = np.linspace(start=1, stop=200, num=100)
                bins = np.linspace(start=1, stop=350, num=350)
                histogram = np.histogram(distances, bins=bins, density=False)[0]
                histogram_list.append(histogram)
                counter += 1
            except FileNotFoundError:
                traceback.print_exc(file=log_file)
    histogram_array = np.asarray(histogram_list)
    if not histogram_list:
        print("Warning: Histogram list is empty. Check log file.")
    if given_algorithm == "alpha" and "low" in path_to_csvs:
        np.save(f"../data/alphafold/histogram_low_conf_{length_range}_not_normed.npy", histogram_array)
    elif given_algorithm == "alpha" and "low" not in path_to_csvs:
        np.save(f"../data/alphafold/histogram_{length_range}_not_normed.npy", histogram_array)
    elif given_algorithm == "rcsb":
        np.save(f"../data/rcsb/histogram_{length_range}_not_normed.npy", histogram_array)


def pdb_to_pcm(log_file: str, given_algorithm: str, length_range: str, path_to_pdbs: str) -> None:
    """
    Convert given PDB files to ProteinContactNetworks.
    @rtype: object
    @param path_to_pdbs: path to directory containing PDB files
    @param length_range: chain length range
    @param given_algorithm: PDB or aF
    @param log_file: log.txt to save errors
    @return: None
    """
    ids_100, ids_200, ids_300 = [], [], []
    distances_100, distances_200, distances_300 = [], [], []
    pdb_ids = None
    n_pdbs = 0
    if given_algorithm == "PDB":
        print("Getting amino acid distances for RCSB data.")
        pdb_ids = concatenate_rcsb_id_files()
        n_pdbs = len(pdb_ids)
    elif given_algorithm == "aF":
        print("Getting amino acid distances for AlphaFold data.")
        pdb_ids = get_uniprot_ids(length_range)
        n_pdbs = len(pdb_ids)
    counter = 0

    with open(log_file, "w") as log:
        for i in range(n_pdbs):
            print(f"At entry {counter + 1} out of {n_pdbs + 1}.")
            try:
                pdb_id = str(pdb_ids[i]).lower()
                print(f"PDB ID is: {pdb_id}")
                pdb_file = None
                if given_algorithm == "PDB":
                    pdb_file = f"{path_to_pdbs}/pdb{pdb_id}.ent"
                elif given_algorithm == "aF":
                    pdb_file = f"{path_to_pdbs}/AF-{pdb_id}-F1-model_v2.pdb"
                print("Successfully opened file.")
            except FileNotFoundError:
                traceback.print_exc(file=log)
                continue

            if os.path.isfile(pdb_file):
                pcm = protein_contact_map.ProteinContactMap(pdb_file)
                alpha_carbons = pcm.get_alpha_carbons
                chain_length = protein_contact_map.get_chain_length(alpha_carbons)

                if chain_length in range(85, 116):
                    print("This PDB is in the length range 85-115. \nGetting amino acid distances.")
                    ids_100.append(pdb_id)
                    distances_100.append(pcm.get_link_lengths(alpha_carbons))
                elif chain_length in range(185, 216):
                    print("This PDB is in the length range 185-215. \nGetting link lengths.")
                    ids_200.append(pdb_id)
                    distances_200.append(pcm.get_link_lengths(alpha_carbons))
                elif chain_length in range(285, 316):
                    print("This PDB is in the length range 285-315. \nGetting link lengths.")
                    ids_300.append(pdb_id)
                    distances_300.append(pcm.get_link_lengths(alpha_carbons))
                else:
                    print("PDB is outside chosen length ranges.")
            else:
                continue

            counter += 1

        print(f"100s: {len(ids_100)}\n200s: {len(ids_200)}\n300s: {len(ids_300)}")
        if given_algorithm == "PDB":
            np.save(f"../data/rcsb/lls_100.npy", distances_100)
            np.save(f"../data/rcsb/lls_200.npy", distances_200)
            np.save(f"../data/rcsb/lls_300.npy", distances_300)
        elif given_algorithm == "aF":
            np.save(f"../data/alphafold/lls_100.npy", distances_100)
            np.save(f"../data/alphafold/lls_200.npy", distances_200)
            np.save(f"../data/alphafold/lls_300.npy", distances_300)
        print("Done.")


def bootstrap(inputfile: str, sample_replacement: bool, length_range: str) -> None:
    """
    Bootstrap data either with or without replacement.
    @param inputfile: name and full path
    @param sample_replacement
    @param length_range
    @return: None
    """
    distances_file = np.load(inputfile, allow_pickle=True)
    distances = np.concatenate(distances_file)
    n_distances = distances.shape[0]
    n_bootstrapping_samples = 1000
    random_number_generator = np.random.default_rng(seed=1234)
    bootstrapped_samples = {}
    if n_distances > 0:
        print("Bootstrapping.")
        for n in range(n_bootstrapping_samples):
            print(f"At step {n + 1}.")
            random_observations_from_sample = random_number_generator.choice(np.arange(n_distances),
                                                                             replace=sample_replacement,
                                                                             size=n_distances)
            bootstrapping_sample = [distances[i] for i in random_observations_from_sample]
            values, counts = np.unique(bootstrapping_sample, return_counts=True)
            for value, count in zip(values, counts):
                try:
                    bootstrapped_samples[value].append(count)
                except KeyError:
                    bootstrapped_samples[value] = []
                    bootstrapped_samples[value].append(count)

            bootstrap_dataframe_tr = pd.DataFrame.from_dict(bootstrapped_samples, orient="index")
            bootstrap_dataframe = bootstrap_dataframe_tr.transpose()
            bootstrap_dataframe_no_nan = bootstrap_dataframe.fillna(0)
            if sample_replacement:
                bootstrap_dataframe_no_nan.to_csv(f"../data/rcsb/bootstrap_{length_range}_raw.csv", index=False)
            elif not sample_replacement:
                bootstrap_dataframe_no_nan.to_csv(f"../data/alphafold/chunk_{length_range}_raw.csv", index=False)

            melted_bootstrap_dataframe = bootstrap_dataframe_no_nan.melt()
            bootstrap_dataframe_stats = melted_bootstrap_dataframe.groupby("variable", as_index=False).agg(
                mean=("value", np.mean), lower_bound=("value", lambda val: np.quantile(val, q=0.05)),
                upper_bound=("value", lambda val: np.quantile(val, q=0.95)))
            if sample_replacement:
                bootstrap_dataframe_stats.to_csv(f"../data/rcsb/bootstrap_{length_range}_stats.csv", index=False)
            elif not sample_replacement:
                bootstrap_dataframe_stats.to_csv(f"../data/alphafold/chunk_{length_range}_stats.csv", index=False)


def get_3d_simulation_files(length_range: str) -> list:
    """
    Get all the 3D simulation adjacency matrix files belonging to the length range
    @param length_range: chain length range
    @return: list containing all simulation files
    """
    return glob.glob(f"../data/simulations/3d/matrices/matrix_{length_range}_*")


def get_3d_simulation_adjacency_matrix(simulation_file: str) -> np.ndarray:
    """
    Opens a 3d simulation file and returns a binary adjacency matrix
    @param simulation_file: 3d simulation file of an adjacency matrix
    @return: numpy array of a binary adjacency matrix
    """
    adjacency_matrix = np.loadtxt(simulation_file)
    adjacency_matrix[adjacency_matrix > 1] = 0
    return adjacency_matrix


def get_3d_simulation_distances(adjacency_matrix: np.ndarray) -> np.ndarray:
    """
    Open 3D simulation adjacency matrix files and compute amino acid distances with statistics
    @param adjacency_matrix: adjacency matrix from 3d simulation files
    @return: numpy array of distances
    """
    distances_list = []
    for row_value in range(len(adjacency_matrix)):
        for col_value in range(len(adjacency_matrix)):
            if adjacency_matrix[row_value][col_value] == 1:
                distance = np.abs(col_value - row_value)
                distances_list.append(distance)
    return np.asarray(distances_list)


def return_3d_simulation_distance_histogram(length_range: str) -> None:
    """
    Open 3D simulation files and return frequencies of amino acid distances
    @param length_range: chain length range
    @return: None
    """
    simulation_files = get_3d_simulation_files(length_range)
    histogram_list = []
    counter = 1
    for simulation_file in simulation_files:
        print(f"Progress: {counter}/{len(simulation_files)}")
        adjacency_matrix = get_3d_simulation_adjacency_matrix(simulation_file)
        distances = get_3d_simulation_distances(adjacency_matrix)
        bins = np.linspace(start=1, stop=350, num=350)
        histogram = np.histogram(distances, bins=bins, density=False)[0]
        histogram_list.append(histogram)
        counter += 1
    histogram_array = np.asarray(histogram_list)
    np.save(f"../data/simulations/3d/histogram_{length_range}_not_normed.npy", histogram_array)
