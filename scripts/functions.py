import pandas as pd
import numpy as np
import protein_contact_map
import traceback
import glob
import os


def pdb_to_adjacency(pdb_file: str, cutoff=8.0) -> tuple:
    """
    Convert given PDB file to an adjacency matrix
    @param pdb_file: PDB file from RCSB or AlphaFold in each length range
    @return: tuple of the adjacency and distance matrices as numpy arrays
    """
    pcm = None
    if cutoff != 8.0:
        pcm = protein_contact_map.ProteinContactMap(pdb_file, cutoff)
    else:
        pcm = protein_contact_map.ProteinContactMap(pdb_file)
    alpha_carbons = pcm.get_alpha_carbons
    distance_array = protein_contact_map.get_distance_array(alpha_carbons)
    return pcm.get_adjacency_matrix(alpha_carbons, distance_array)


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


def return_distance_histogram(log_file: str, given_algorithm: str, length_range: str, path_to_csvs: str) -> None:
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
                adjacency_matrix = pdb_to_adjacency(clean_pdb_filename)
                distances = get_distances(adjacency_matrix)
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