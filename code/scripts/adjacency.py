""" Functions for plotting an adjacency matrix from 3D simulations or PDBs"""
import matplotlib
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import protein_contact_map
import theory_functions
from plot_functions import _COLOUR_PALETTE


def get_simulation_matrix(matrix_file: str) -> np.ndarray:
    """
    Open a 3D simulation file and return the adjacency matrix
    @param matrix_file: matrix_<length_range>_<file_no>.txt
    @return: adjacency matrix as a Numpy array
    """
    adjacency_matrix = np.loadtxt(matrix_file)
    adjacency_matrix[adjacency_matrix > 1] = 0
    return adjacency_matrix


def get_pdb_matrix(pdb_file: str) -> np.ndarray:
    """
    Open a PDB file, create a graph and return its adjacency matrix
    @param pdb_file: PDB file from either RCSB or AlphaFold (default)
    @return: adjacency matrix as a Numpy array
    """
    pcm = protein_contact_map.ProteinContactMap(pdb_file)
    alpha_carbons = pcm.get_alpha_carbons
    protein_graph = pcm.get_protein_graph(alpha_carbons)
    return theory_functions.get_adjacency_matrix(protein_graph)


def set_adjacency_matrix_ticks(plot: matplotlib.axes.Axes) -> None:
    """
    Format tick labels for adjacency matrix plots
    @param plot: seaborn heatmap
    @return:
    """
    for index, label in enumerate(plot.get_xticklabels()):
        if index % 5 == 0:
            label.set_visible(True)
            label.set_rotation(360)
            label.set_font("Helvetica")
            label.set_fontsize(14)
        else:
            label.set_visible(False)
    for index, label in enumerate(plot.get_yticklabels()):
        if index % 5 == 0:
            label.set_visible(True)
            label.set_font("Helvetica")
            label.set_fontsize(14)
        else:
            label.set_visible(False)


def plot_adjacency_matrix(file: str, data_type: str) -> None:
    """
    Plot adjacency matrix for given protein file (or simulation matrix)
    @param file: PDB or matrix file
    @param data_type: PDB or SIM
    @return: None
    """
    adjacency_matrix = []
    if data_type == "SIM":
        adjacency_matrix = get_simulation_matrix(file)
    elif data_type == "PDB":
        adjacency_matrix = get_pdb_matrix(file)

    plt.figure(figsize=(8, 8))
    colormap = [_COLOUR_PALETTE["NO_CONTACT"], _COLOUR_PALETTE["CONTACT"]]
    heatmap = sns.heatmap(adjacency_matrix, cmap=colormap, cbar=False)
    set_adjacency_matrix_ticks(heatmap)
    heatmap.set_xlabel("Amino acid", fontsize=20, font="Helvetica")
    heatmap.set_ylabel("Amino acid", fontsize=20, font="Helvetica")
    sns.despine()
    plt.savefig(f"../plots/adjacency_matrices/{data_type}_matrix.jpeg", dpi=900)
    plt.show()
