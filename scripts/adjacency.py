""" Functions for plotting an adjacency matrix from 3D simulations or PDBs"""
import matplotlib
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import protein_contact_map
import functions
from colour_palette import _COLOUR_PALETTE
import pandas as pd
import traceback
import argparse
import os
import sys


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


def get_all_adjacency_matrices(path_to_csv: str) -> str:
    """
    Open secondary structure csv file and use pdb files to store all adjacency matrices in a numpy file
    @param path_to_csv: full path to secondary structure csv file
    @return: chain length range as a string 
    """
    dataframe = pd.read_csv(path_to_csv) 
    length_range = path_to_csv.split("_")[-1].replace(".csv", "")
    pdb_filenames = dataframe["filename"].tolist()
    all_adjacency_matrices_list = []
    counter = 1
    for pdb_filename in pdb_filenames:
        print(f"Progress: {counter}/{len(pdb_filenames)}")
        try:
            adjacency_matrix = get_pdb_matrix(pdb_filename)
            rows, columns = adjacency_matrix.shape
            set_correct_rows = np.vstack([adjacency_matrix, np.zeros((350 - rows, columns), dtype=adjacency_matrix.dtype)])
            new_rows, _ = set_correct_rows.shape
            reshaping_columns = np.zeros((new_rows, 350 - columns))
            square_matrix = np.hstack((set_correct_rows, reshaping_columns))
            all_adjacency_matrices_list.append(square_matrix)

        except FileNotFoundError as err:
            print(err)
        all_adjacency_matrices = np.array(all_adjacency_matrices_list)
        counter += 1
    if not all_adjacency_matrices:           
        print("Warning: Histogram list is empty. Check log file.")
    np.save(f"../data/rcsb/adjacency_matrix_{length_range}.npy", all_adjacency_matrices)
    return length_range
    

# Get average adjacency matrix:
def get_average_adjacency_matrix(length_range: str):
    """
    compute the average adjacency matrix for a given chain length range
    @param length_range: chain length range as a string
    @return: average adjacency matrix
    """
    total_matrix = np.load(f"../data/rcsb/adjacency_matrix_{length_range}.npy", allow_pickle=True)
    n_total_counts = np.sum(total_matrix)
    return np.sum(total_matrix, axis=0) / n_total_counts


def plot_average_matrix(length_range: str, average_matrix: np.ndarray):
    """
    Plot average adjacency matrix for given chain length range
    @param length_range: chain length range as as tring
    @return: None
    """
    average_adjacency_matrix = remove_diagonal_elements(average_matrix)
    plt.figure(figsize=(6, 6))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=1.8, font="Helvetica")
    colormap = [_COLOUR_PALETTE["NO_CONTACT"], _COLOUR_PALETTE["CONTACT"]]
    heatmap = sns.heatmap(average_adjacency_matrix, cmap=colormap, cbar=False)
    set_adjacency_matrix_ticks(heatmap)
    heatmap.set_xlabel("Amino acid")
    heatmap.set_ylabel("Amino acid")
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"../plots/adjacency_matrices/average_matrix_{length_range}.jpeg", dpi=900)


def command_line_arguments() -> argparse.Namespace:
    """
    Parser for command line arguments
    @return: command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_file",
                        type=str,
                        help="secondary structure csv file")                 
    return parser.parse_args()

def main():
    arguments = command_line_arguments()
    csv_file = arguments.csv_file
    if not os.path.isfile(csv_file):
        print(f"The input file {csv_file} does not exist")
        sys.exit()
    length_range = get_all_adjacency_matrices(csv_file)
    average_matrix = get_average_adjacency_matrix(length_range)
    plot_average_matrix(length_range, average_matrix)


if __name__ == "__main__":
    main()
    