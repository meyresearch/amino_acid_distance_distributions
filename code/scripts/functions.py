"""
Functions for getting link lengths and bootstrapping.
"""
import glob
import traceback

import networkx as nx

import pcm
import numpy as np
import os
import pandas as pd
import argparse


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
    id_dataframe = pd.read_csv(f"../data/alphafold/uniprot_ids/uniprot_{sequence_length}s.xlsx")
    return id_dataframe["Entry"].to_numpy()


def commandline_arguments() -> argparse.Namespace:
    """
    Parser for commandline arguments
    @return: commandline arguments from user
    """

    parser = argparse.ArgumentParser(description="Get amino acid distance distributions.")
    parser.add_argument("algorithm", type=str, choices=["PDB", "aF", "SIM", "BS", "C"],
                        help="Get amino acid distances from RCSB (PDB), AlphaFold (aF) or simulations (SIM), bootstrap "
                             "(BS) or chunk (C)")
    parser.add_argument("-r", dest="length_range", type=str, choices=[None, "100", "200", "300"], help="Chain length "
                                                                                                       "range")
    parser.add_argument("-p", dest="path_to_pdbs", type=str, help="Full path to directory containing PDB files")
    parser.add_argument("-i", dest="inputfile", type=str, help="Full path and name of input file")
    cl_arguments = parser.parse_args()

    if cl_arguments.algorithm == "BS" and cl_arguments.inputfile is None:
        parser.error("BS requires -i")
    elif cl_arguments.algorithm == "C" and cl_arguments.inputfile is None:
        parser.error("C requires -i")
    elif cl_arguments.algorithm == "PDB" and cl_arguments.length_range is not None:
        parser.error("PDB requires -r=None")
    elif cl_arguments.algorithm == "PDB" and cl_arguments.path_to_pdbs is None:
        parser.error("PDB requires -p")
    elif cl_arguments.algorithm == "aF" and cl_arguments.length_range is None:
        parser.error("aF requires -r")
    elif cl_arguments.algorithm == "aF" and cl_arguments.path_to_pdbs is None:
        parser.error("aF requires -p")
    elif cl_arguments.algorithm == "SIM" and cl_arguments.length_range is None:
        parser.error("SIM requires -r")
    return cl_arguments


def pdb_to_pcm(log_file: str, given_algorithm: str, length_range: str, path_to_pdbs: str) -> None:
    """
    Convert given PDB files to ProteinContactNetworks.
    @rtype: object
    @param path_to_pdbs: path to directory contiaining PDB files
    @param length_range: chain length range
    @param given_algorithm: PDB or aF
    @param log_file: log.txt to save errors
    @return:
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
                protein_contact_map = pcm.ProteinContactMap(pdb_file)
                alpha_carbons = protein_contact_map.get_alpha_carbons
                chain_length = pcm.get_chain_length(alpha_carbons)

                if chain_length in range(85, 116):
                    print("This PDB is in the length range 85-115. \nGetting amino acid distances.")
                    ids_100.append(pdb_id)
                    distances_100.append(protein_contact_map.get_link_lengths(alpha_carbons))
                elif chain_length in range(185, 216):
                    print("This PDB is in the length range 185-215. \nGetting link lengths.")
                    ids_200.append(pdb_id)
                    distances_200.append(protein_contact_map.get_link_lengths(alpha_carbons))
                elif chain_length in range(285, 316):
                    print("This PDB is in the length range 285-315. \nGetting link lengths.")
                    ids_300.append(pdb_id)
                    distances_300.append(protein_contact_map.get_link_lengths(alpha_carbons))
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


def get_3d_sim_files(length_range: str) -> list:
    """
    Get all the 3D simulation adjacency matrix files belonging to the length range
    @param length_range: chain length range
    @return: list containing all simulation files
    """
    return glob.glob(f"../data/simulations/3d/matrices/matrix_{length_range}_*")


def sim_to_pcm(length_range: str) -> None:
    """
    Open 3D simulation adjacency matrix files and compute amino acid distances with statistics
    @param length_range: chain length range
    @return: None
    """
    files = get_3d_sim_files(length_range)
    all_distances = []
    for file in files:
        distances_list = []
        adjacency_matrix = np.loadtxt(file)
        adjacency_matrix[adjacency_matrix > 1] = 0
        protein_graph = nx.from_numpy_matrix(adjacency_matrix)
        distances = list(protein_graph.edges())
        for distance in distances:
            distances_list.append(abs(distance[0] - distance[1]))
        all_distances.append(np.asarray(distances_list))
        np.save(f"../data/simulations/3d/lls_{length_range}.npy", all_distances)


def get_simulation_stats(length_range: str):
    """
    Open 3D simulation files and return frequencies of amino acid distances
    @param length_range: chain length range
    @return: pd.DataFrame
    """
    simulation_samples = np.load(f"../data/simulations/3d/lls_{length_range}.npy", allow_pickle=True)
    distances_dict = {}
    for sample in simulation_samples:
        values, counts = np.unique(sample, return_counts=True)
        for value, count in zip(values, counts):
            try:
                distances_dict[value].append(count)
            except KeyError:
                distances_dict[value] = []
                distances_dict[value].append(count)
    dataframe_from_dict = pd.DataFrame.from_dict(distances_dict, orient="index").reset_index()
    dataframe_no_nans = dataframe_from_dict.fillna(0)
    dataframe_no_nans.to_csv(f"../data/simulations/3d/simulation_{length_range}_raw.csv")
    melted_dataframe = dataframe_no_nans.melt(id_vars="index")
    stats_dataframe = melted_dataframe.groupby("index", as_index=False).agg(mean=("value", np.mean),
                                                                            lower_bound=("value", lambda val:
                                                                            np.quantile(val, q=0.05)),
                                                                            upper_bound=("value", lambda val:
                                                                            np.quantile(val, q=0.95)))
    stats_dataframe.to_csv(f"../data/simulations/3d/simulation_{length_range}_stats.csv")


def compute_3d_simulation_distribution(length_range: str) -> None:
    """

    @param length_range:
    @return:
    """
    sim_to_pcm(length_range=length_range)
    get_simulation_stats(length_range=length_range)

