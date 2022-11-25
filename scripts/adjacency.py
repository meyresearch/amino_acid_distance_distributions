""" Functions for plotting an adjacency matrix from 3D simulations or PDBs"""
import matplotlib
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import protein_contact_map
import functions
from colour_palette import _COLOUR_PALETTE


def get_pdb_matrix(pdb_file: str) -> np.ndarray:
    """
    Open a PDB file and return its adjacency matrix
    @param pdb_file: PDB file from either RCSB or AlphaFold (default)
    @return: adjacency matrix as a Numpy array
    """
    return functions.pdb_to_adjacency(pdb_file)[1]


def set_adjacency_matrix_ticks(plot: matplotlib.axes.Axes) -> None:
    """
    Format tick labels for adjacency matrix plots
    @param plot: seaborn heatmap
    @return:
    """
    for index, label in enumerate(plot.get_xticklabels()):
        if index % 4 == 0:
            label.set_visible(True)
            label.set_rotation(360)
            label.set_font("Helvetica")
        else:
            label.set_visible(False)
    for index, label in enumerate(plot.get_yticklabels()):
        if index % 4 == 0:
            label.set_visible(True)
            label.set_font("Helvetica")
        else:
            label.set_visible(False)


def remove_diagonal_elements(matrix: np.ndarray):
    """
    Remove diagonal elements and add contacts to first off-diagonals
    @param matrix: adjacency matrix from simulations
    @return: simulation matrix with diagonals set to 1
    """
    np.fill_diagonal(matrix, 0)
    off_diagonal = np.ones(len(np.diag(matrix, 1)))
    np.fill_diagonal(matrix[1:], off_diagonal)
    np.fill_diagonal(matrix[:, 1:], off_diagonal)
    return matrix


def plot_adjacency_matrix(file: str, data_type: str) -> None:
    """
    Plot adjacency matrix for given protein file (or simulation matrix)
    @param file: PDB or matrix file
    @param data_type: PDB or SIM
    @return: None
    """
    adjacency_matrix = []
    if data_type == "sim":
        simulation_matrix = functions.get_3d_simulation_adjacency_matrix(file)
        adjacency_matrix = remove_diagonal_elements(simulation_matrix)
    elif data_type == "pdb":
        pdb_matrix = get_pdb_matrix(file)
        adjacency_matrix = remove_diagonal_elements(pdb_matrix)
    plt.figure(figsize=(6, 6))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=1.8, font="Helvetica")
    colormap = [_COLOUR_PALETTE["NO_CONTACT"], _COLOUR_PALETTE["CONTACT"]]
    heatmap = sns.heatmap(adjacency_matrix, cmap=colormap, cbar=False)
    set_adjacency_matrix_ticks(heatmap)
    heatmap.set_xlabel("Amino acid")
    heatmap.set_ylabel("Amino acid")
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"../plots/adjacency_matrices/{data_type}_matrix.jpeg", dpi=900)
    plt.show()


# skeleton code for average matrix

# read in secondary_structure_{range}.csv
# dataframe = pd.read_csv(path_to_csv) 
# pdb_filenames = dataframe["filename"].tolist()

# Getting adjacency matrix data:
# with open("log.txt", "w") as log_file:
#   for pdb_filename in pdb_filenames:
        # given_algortihm = [rcsb/alpha]
#       if given_algorithm == "alpha":
#               clean_pdb_filename = pdb_file.replace("/home/jguven/Projects/sequence_distance_distribution", "..")
#       else:
#               clean_pdb_filename = pdb_filename
#       try:
#           adjacency_matrix = functions.adjacency_matrix(clean_pdb_filename)[1]
#           all_adjacency_matrices_list.append(adjacency_matrix)
#       except FileNotFoundError:
#           traceback.print_exc(file=log_file)
#       all_adjacency_matrices = np.array(all_adjacency_matrices_list)
#    if not histogram_list:
#        print("Warning: Histogram list is empty. Check log file.")
#   if given_algorithm == "alpha":
#        np.save(f"../data/alphafold/adjacency_matrix_{length_range}.npy", all_adjacency_matrices)
#   elif given_algorithm == "rcsb":
#        np.save(f"../data/rcsb/adjacency_matrix_{length_range}.npy", all_adjacency_matrices)
       
# Get average adjacency matrix:

# total_matrix = np.load(f"../data/rcsb/adjacency_matrix_{length_range}.npy")
# n_total_counts = np.sum(total_matrix)
# average_adjacency_matrix = np.sum(total_matrix, axis=0) / n_total_counts

# plot