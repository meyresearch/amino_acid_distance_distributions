import pandas as pd
import functions
import numpy as np


def return_distance_histogram(log_file: str, given_algorithm: str, length_range: str, path_to_csvs: str, threshold: float) -> np.ndarray:
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
        for pdb_file in pdb_files:
            # print(f"Progress: {counter}/{len(pdb_files)}")
            if given_algorithm == "alphafold":
                clean_pdb_filename = pdb_file.replace("/home/jguven/Projects/sequence_distance_distribution", "..")
            else:
                clean_pdb_filename = pdb_file
            try:
                adjacency_matrix = pdb_to_adjacency(clean_pdb_filename)[1]
                distances = functions.get_distances(adjacency_matrix)
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
    if given_algorithm == "alphafold":
        np.save(f"histogram_alphafold_{length_range}_{threshold}.npy", histogram_array)
    elif given_algorithm == "rcsb":
        np.save(f"histogram_{length_range}_{threshold}.npy", histogram_array)


lengths = [100, 200, 300]
algorithms = ["rcsb", "alphafold"]
filename = ""
for algorithm in algorithms:
    if algorithm == "alphafold":
        filename = f"unique_secondary_structures_{length}.csv"
    elif algorithm == "rcsb":
        filename = f"secondary_structures_{length}.csv"
    for length in lengths:
    return_distance_histogram(log_file="log.txt", given_algorithm=algorithm, length_range=str{length}, path_to_csvs=f"../data/{algorithm}/{filename}")
