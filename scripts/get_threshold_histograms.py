import pandas as pd
import functions
import numpy as np
from tqdm import tqdm
import protein_contact_map
import traceback

def pdb_to_adjacency(pdb_file: str, threshold: float) -> tuple:
    """
    Convert given PDB file to an adjacency matrix
    @param pdb_file: PDB file from RCSB or AlphaFold in each length range
    @param threshold: threshold value to use for counting contacts
    @return: tuple of the adjacency and distance matrices as numpy arrays
    """
    pcm = protein_contact_map.ProteinContactMap(pdb_file, default_threshold=threshold)
    alpha_carbons = pcm.get_alpha_carbons
    distance_array = protein_contact_map.get_distance_array(alpha_carbons)
    adjacency_matrix = pcm.get_adjacency_matrix(alpha_carbons, distance_array)
    return adjacency_matrix


def return_distance_histogram(log_file: str, given_algorithm: str, length_range: str,
                              path_to_csvs: str, threshold: float) -> np.ndarray:
    """
    Compute the amino acid distance distribution for PDB files in given range from adjacency matrix
    and save in a numpy file.
    @param log_file: file to save exceptions in
    @param given_algorithm: alphafold or rcsb
    @param length_range: 100, 200 or 300
    @param path_to_csvs: full path to csv files
    @param threshold: threshold value to use for counting contacts
    @return: None
    """
    dataframe = pd.read_csv(path_to_csvs)
    pdb_files = dataframe["filename"].to_numpy()
    histogram_list = []
    counter = 1
    with open(log_file, "w") as log_file:

        for pdb_file in tqdm(pdb_files):
            # print(f"Progress: {counter}/{len(pdb_files)}")
            if given_algorithm == "alphafold":
                clean_pdb_filename = pdb_file.replace("/home/jguven/Projects/sequence_distance_distribution", "..")
            else:
                clean_pdb_filename = pdb_file
            try:
                adjacency_matrix = pdb_to_adjacency(clean_pdb_filename, threshold)
                distances = functions.get_distances(adjacency_matrix)
                bins = np.linspace(start=1, stop=350, num=350)
                histogram = np.histogram(distances, bins=bins, density=False)[0]
                histogram_list.append(histogram)
                counter += 1
            except FileNotFoundError:
                traceback.print_exc(file=log_file)
    histogram_array = np.asarray(histogram_list)
    if not histogram_list:
        print("Warning: Histogram list is empty. Check log file.")
    print("Saving histograms.")
    if given_algorithm == "alphafold":
        np.save(f"../data/alphafold/threshold_histograms/histogram_{length_range}_{threshold}.npy", histogram_array)
    elif given_algorithm == "rcsb":
        np.save(f"../data/rcsb/threshold_histograms/histogram_{length_range}_{threshold}.npy", histogram_array)


lengths = [100, 200, 300]
algorithms = ["rcsb", "alphafold"]
thresholds = [10, 15, 21]
filename = ""
for algorithm in algorithms:
    print(f"Algorithm: {algorithm}")
    for length in lengths:
        print(f"Chain length range: {length}")
        for dc in thresholds:
            print(f"Threshold: {dc}")
            if algorithm == "alphafold":
                filename = f"unique_secondary_structures_{length}.csv"
            elif algorithm == "rcsb":
                filename = f"secondary_structures_{length}.csv"
            return_distance_histogram(log_file="log.txt", given_algorithm=algorithm, length_range=str(length),
                                      path_to_csvs=f"../data/{algorithm}/{filename}", threshold=dc)
