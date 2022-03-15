"""
Functions used for plotting amino acid distance distributions.
"""
import argparse
import glob

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pcm
import theory_functions
import seaborn as sns

_COLOUR_PALETTE = {"PDB_SCATTER": "#006374",
                   "SIM_SCATTER": "#fbafe4",
                   "ALPHA_SCATTER": "#1ac938",
                   "THEORY": "#006374",
                   "RESIDUALS": "#fbafe4",
                   "DATABANK": "#006374",
                   "USED": "#fbafe4",
                   "CONTACT": "#fbafe4",
                   "NO_CONTACT": "#006374"}


def create_plot_label(length_range: str, algorithm: str):
    """
    Create custom labels for Alphafold PDB data and RCSB PDB data with the given link length
    @param algorithm: data source from RCSB (bootstrapped) or AlphaFold (chunked)
    @param length_range: chain length range
    @return:
    """
    plot_label = " "
    if algorithm == "BS" or algorithm == "BOTH":
        plot_label = f"PDB {length_range}s"
    elif algorithm == "C":
        plot_label = f"AF PDB {length_range}"
    return plot_label


def grid_plot_distances(dimensionality_range: np.ndarray, exponent_range: np.ndarray, n_points: int,
                        half_n_harmonic: float, plotting_sumrange: np.ndarray, normalised_means: np.ndarray,
                        distances: np.ndarray, normalised_sim_means: np.ndarray, sim_distances: np.ndarray,
                        lower_cl: np.ndarray, upper_cl: np.ndarray, length_range: str, algorithm: str) -> None:
    """
    Plot amino acid distances and theory plot with different values of exponent and dimensionality scaling constant
    @param algorithm: data source from RCSB (bootstrapped) or AlphaFold (chunked)
    @param length_range: chain length range
    @param sim_distances: amino acid distances for simulated data (x-axis)
    @param upper_cl: upper confidence level
    @param lower_cl: lower confidence level
    @param normalised_sim_means: normalised simulated means of amino acid distance frequencies
    @param normalised_means: normalised means of PDB amino acid distance frequencies
    @param distances: amino acid distances for PDB data (x-axis)
    @param plotting_sumrange: Range to be used in sum for the theory plot
    @param half_n_harmonic: Harmonic number for N/2
    @param n_points: number of data points
    @param exponent_range: range of exponent values (constant a)
    @type dimensionality_range: range of dimensionality scaling constant (constant A)
    @return: None
    """
    plot_label = create_plot_label(length_range, algorithm)
    fig = plt.figure(figsize=(16, 8))
    sns.set(context="notebook", palette="colorblind", style='ticks', font='Helvetica')
    ax = fig.subplots(len(exponent_range), len(dimensionality_range), sharex=True, sharey=True)
    for row in range(len(exponent_range)):
        for col in range(len(dimensionality_range)):
            theory = [
                theory_functions.amino_acid_distance_distribution(s, n_points, half_n_harmonic, exponent_range[row],
                                                                  dimensionality_range[col]) for s in plotting_sumrange]
            ax[row][col].scatter(distances, normalised_means, s=10, c=_COLOUR_PALETTE["PDB_SCATTER"], label=plot_label)
            ax[row][col].fill_between(distances, upper_cl, lower_cl, color=_COLOUR_PALETTE["PDB_SCATTER"], alpha=0.4,
                                      label="95% CL", zorder=-100)
            ax[row][col].scatter(sim_distances, normalised_sim_means, s=10, marker="^",
                                 c=_COLOUR_PALETTE["SIM_SCATTER"], label=f"SIM {length_range}s", zorder=-50)
            ax[row][col].plot(plotting_sumrange, theory, c=_COLOUR_PALETTE["THEORY"], label="Theory")
            ax[row][col].set_yscale("log")
            ax[row][col].set_xscale("log")
            ax[row][col].set_ylim(0.0001, 1)
            ax[row][col].set_xlim(2, n_points / 2)
            ax[row][-1].set_ylabel(f"a = {exponent_range[row]}", fontsize=13, rotation=0, labelpad=21)
            ax[row][-1].yaxis.set_label_position("right")
            ax[0][col].set_title(f"A = {dimensionality_range[col]:2f}")
    plt.legend(bbox_to_anchor=[1.57, 4.65], fontsize=9)
    fig.text(0.5, 0.025, "s", ha="center", fontsize=15.5)
    fig.text(0.005, 0.5, "P(s)", va="center", rotation="vertical", fontsize=15.5)
    plt.subplots_adjust(left=0.06, bottom=0.08, top=0.95, wspace=0.05, right=0.909)
    plt.savefig(f"../plots/supplementary_information/{algorithm}_{length_range}.pdf")


def grid_plot_residuals(algorithm: str, length_range: str, dimensionality_range, exponent_range: np.ndarray,
                        n_points: int, half_n_harmonic_number: float, plotting_sumrange: np.ndarray,
                        normalised_means: np.ndarray) -> None:
    """
    Plot residuals grid and print out statistics.
    @param algorithm:
    @param length_range:
    @param dimensionality_range:
    @param exponent_range:
    @param n_points:
    @param half_n_harmonic_number:
    @param plotting_sumrange:
    @param normalised_means:
    @return:
    """
    fig = plt.figure(figsize=(16, 8))
    sns.set(context="notebook", palette="colorblind", style='ticks', font='Helvetica')
    ax = fig.subplots(len(exponent_range), len(dimensionality_range), sharex=True, sharey=True)
    for row in range(len(exponent_range)):
        for col in range(len(dimensionality_range)):
            residuals = theory_functions.plotting_statistics(dimensionality_range[col], exponent_range[row],
                                                             n_points, half_n_harmonic_number, plotting_sumrange,
                                                             normalised_means)[0]

            ax[row][col].scatter(plotting_sumrange, residuals, marker=".", c=_COLOUR_PALETTE["RESIDUALS"],
                                 label="Residuals", zorder=10)
            ax[row][col].hlines(0, plotting_sumrange[0], plotting_sumrange[-1], color=_COLOUR_PALETTE["THEORY"],
                                label="Theory")
            ax[row][-1].set_ylabel(f"a = {exponent_range[row]}", fontsize=13, rotation=0, labelpad=21)
            ax[row][-1].yaxis.set_label_position("right")
            ax[row][col].set_ylim(-0.025, 0.025)
            ax[0][col].set_title(f"A = {dimensionality_range[col]:2f}")
            ax[row][col].legend(fontsize=9)

    fig.text(0.5, 0.025, "s", ha="center", fontsize=15.5)
    fig.text(0.005, 0.5, "Residuals", va="center", rotation="vertical", fontsize=15.5)
    plt.subplots_adjust(left=0.06, bottom=0.08, top=0.95, wspace=0.1, right=0.95)
    plt.savefig(f"../plots/supplementary_information/{algorithm}_{length_range}_r.pdf")


def parse_command_line_arguments() -> argparse.Namespace:
    """
    Parse CL arguments for plotting.
    @return: Namespace containing CL arguments
    """
    parser = argparse.ArgumentParser(description="Plot amino acid distances and residuals.")
    parser.add_argument("-r", dest="length_range", type=str, choices=["100", "200", "300", None],
                        help="chain length range to be plotted")
    parser.add_argument("algorithm", type=str, choices=["BS", "C", "BOTH", "2D-SIM", "A", "B"],
                        help="get distances from bootstrapping (BS), chunking (C), compare both, distances from 2D "
                             "simulations, plot adjacency matrix (A) or plot bar plots (B)")
    parser.add_argument("-t", dest="data_type", type=str, choices=["PDB", "SIM"], help="data type for adjacency matrix")
    parser.add_argument("-f", dest="file", type=str, help="PDB file or simulation matrix for adjacency matrix")
    parser.add_argument("--d-begin", dest="start_dimensionality", type=float, help="starting value for dimensionality "
                                                                                   "constant (A)")
    parser.add_argument("--d-end", dest="end_dimensionality", type=float, help="last value for dimensionality "
                                                                               "constant (A)")
    parser.add_argument("--e-begin", dest="start_exponent", type=int, help="starting value for exponent (constant a)")
    parser.add_argument("--e-end", dest="end_exponent", type=int, help="starting value for exponent (constant a)")
    parser.add_argument("--e-step", dest="step_exponent", type=int, nargs="?", const=1, default=1,
                        help="step size for constant a range, default=1")
    parser.add_argument("--d-step", dest="step_dimensionality", type=float, nargs="?", const=0.0001, default=0.0001,
                        help="step size for constant A range, default=0.0001")
    parser.add_argument("--d-BS", dest="d_bs", type=float, help="constant A for BS algorithm in comparison plot")
    parser.add_argument("--d-C", dest="d_c", type=float, help="constant A for C algorithm in comparison plot")
    parser.add_argument("--e-BS", dest="e_bs", type=float, help="constant a for BS algorithm in comparison plot")
    parser.add_argument("--e-C", dest="e_c", type=float, help="constant a for C algorithm in comparison plot")
    arguments = parser.parse_args()
    check_required_arguments(parser, arguments.algorithm, arguments.length_range, arguments.start_dimensionality,
                             arguments.end_dimensionality, arguments.start_exponent, arguments.end_exponent,
                             arguments.d_bs, arguments.e_bs, arguments.d_c, arguments.e_c, arguments.data_type,
                             arguments.file)
    return arguments


def check_required_arguments(argument_parser: argparse.ArgumentParser, given_algorithm: str, given_length: str,
                             starting_d: str, ending_d: str, starting_e: str, ending_e: str, bootstrap_d: str,
                             bootstrap_e: str, chunk_d: str, chunk_e: str,
                             data_type: str, file: str) -> None:
    """
    Check given CL arguments so that all requirements are met
    @param file: PDB file or simulation matrix
    @param data_type: SIM or PDB
    @param given_length: given chain length range
    @param argument_parser: parser for CL arguments
    @param given_algorithm: either BS, C or BOTH
    @param starting_d: starting value for dimensionality (constant A)
    @param ending_d: ending value for dimensionality (constant A)
    @param starting_e: starting value for exponent (constant a)
    @param ending_e: ending value for exponent (constant a)
    @param bootstrap_d: dimensionality constant (A) for bootstrapping algorithm
    @param bootstrap_e: exponent constant (a) for bootstrapping algorithm
    @param chunk_d: dimensionality constant (A) for chunking algorithm
    @param chunk_e: exponent constant (a) for chunking algorithm
    @return: None
    """
    if (given_algorithm == "BS" and starting_d is None) or (given_algorithm == "BS" and ending_d is None) or \
            (given_algorithm == "BS" and starting_e is None) or (given_algorithm == "BS" and ending_e is None) \
            or (given_algorithm == "BS" and given_length is None):
        argument_parser.error("BS requires -r, --d-begin, --d-end, --e-begin, --e-end")
    elif (given_algorithm == "C" and starting_d is None) or (given_algorithm == "C" and ending_d is None) or \
            (given_algorithm == "C" and starting_e is None) or (given_algorithm == "C" and ending_e is None) \
            or (given_algorithm == "C" and given_length is None):
        argument_parser.error("BS requires -r, --d-begin, --d-end, --e-begin, --e-end")
    elif (given_algorithm == "BOTH" and bootstrap_d is None) or (given_algorithm == "BOTH" and chunk_d is None) \
            or (given_algorithm == "BOTH" and bootstrap_e is None) or (given_algorithm == "BOTH" and chunk_e is None) \
            or (given_algorithm == "BOTH" and given_length is None):
        argument_parser.error("BOTH requires -r, --d-BS, --d-C, --e-BS, --e-C")
    elif given_algorithm == "B" and given_length is not None:
        argument_parser.error("B requires -r=None")
    elif (given_algorithm == "A" and data_type is None) or (given_algorithm == "A" and file is None) or \
            (given_algorithm == "A" and given_length is not None):
        argument_parser.error("A requires -t, -f")


def get_dataframe(arguments: argparse.Namespace) -> pd.DataFrame:
    """
    Get distribution of amino acid distances from csv.
    @param arguments: command line arguments
    @return: dataframe of bootstrapped/chunked data
    """
    dataframe = None
    if arguments.algorithm == "BS":
        dataframe = pd.read_csv(f"../data/rcsb/bootstrap_{arguments.length_range}_stats.csv")
    elif arguments.algorithm == "C":
        dataframe = pd.read_csv(f"../data/alphafold/chunk_{arguments.length_range}_stats.csv")
    return dataframe


def get_data_for_plotting(dataframe: pd.DataFrame) -> tuple:
    """
    Get protein data for plotting from dataframe
    @param dataframe: contains amino acid distances and confidence levels
    @return: tuple of number of datapoints, normalised means and confidence level bounds
    """
    try:
        distances = dataframe["variable"].to_numpy()
    except KeyError:
        distances = dataframe["index"].to_numpy()
    means = dataframe["mean"].to_numpy()
    lower_bound = dataframe["lower_bound"].to_numpy()
    upper_bound = dataframe["upper_bound"].to_numpy()
    n_datapoints = int(distances[-1] + 1)
    return n_datapoints, means / np.sum(means), lower_bound / np.sum(means), upper_bound / np.sum(means)


def create_grid_plots(arguments: argparse.Namespace, pdb_dataframe: pd.DataFrame, sim_dataframe: pd.DataFrame) -> None:
    """
    Function to bring together everything needed to plot the grid plots of amino acid distance distributions
    and residuals
    @param arguments: command line arguments
    @param pdb_dataframe: dataframe containing data for PDBs
    @param sim_dataframe: 3D simulation data for plotting
    @return:
    """
    pdb_plotting_tuple = get_data_for_plotting(pdb_dataframe)
    pdb_distances = pdb_dataframe["variable"].to_numpy()
    half_n_harmonic = theory_functions.harmonic_number(n_numbers=(pdb_plotting_tuple[0] // 2))
    pdb_plotting_sum_range = np.array(range(int(pdb_distances[0]), pdb_plotting_tuple[0]))
    sim_distances = sim_dataframe["index"].to_numpy()
    sim_plotting_tuple = get_data_for_plotting(sim_dataframe)

    dimensionality_range = np.arange(arguments.start_dimensionality,
                                     arguments.end_dimensionality,
                                     arguments.step_dimensionality)
    exponent_range = np.arange(arguments.start_exponent,
                               arguments.end_exponent,
                               arguments.step_exponent)
    grid_plot_distances(dimensionality_range=dimensionality_range,
                        exponent_range=exponent_range,
                        n_points=pdb_plotting_tuple[0],
                        half_n_harmonic=half_n_harmonic,
                        plotting_sumrange=pdb_plotting_sum_range,
                        normalised_means=pdb_plotting_tuple[1],
                        distances=pdb_distances,
                        normalised_sim_means=sim_plotting_tuple[1],
                        sim_distances=sim_distances,
                        lower_cl=pdb_plotting_tuple[2],
                        upper_cl=pdb_plotting_tuple[3],
                        length_range=arguments.length_range,
                        algorithm=arguments.algorithm)
    grid_plot_residuals(algorithm=arguments.algorithm,
                        length_range=arguments.length_range,
                        dimensionality_range=dimensionality_range,
                        exponent_range=exponent_range,
                        n_points=pdb_plotting_tuple[0],
                        half_n_harmonic_number=half_n_harmonic,
                        plotting_sumrange=pdb_plotting_sum_range,
                        normalised_means=pdb_plotting_tuple[1])
    plt.show()


def get_data_for_bars(path_to_data: str) -> tuple:
    """
    Read csv file from path and return numpy arrays as tuple
    @param path_to_data:
    @return:
    """
    dataframe = pd.read_csv(path_to_data)
    bins = dataframe["Bins"].to_numpy()
    frequencies = dataframe["Number"].to_numpy()
    return bins, frequencies


def calculate_adjusted_frequency(frequency_to_adjust: np.ndarray, bottom_frequency: np.ndarray) -> np.ndarray:
    """

    @param frequency_to_adjust: "top" bars i.e. either all SwissProt or all RCSB frequencies
    @param bottom_frequency: "bottom" bars i.e. the used PDBs' frequencies
    @return: difference between top and bottom bars
    """
    return frequency_to_adjust - bottom_frequency


def create_bar_plots() -> None:
    """
    Bring all stats together and plot subplots of PDB statistics
    @return: None
    """
    path = "../data/pdb_statistics/"
    alphafold_file = "AlphaFold_used.csv"
    pdb_file = "ll_used.csv"
    swissprot_file = "UniProt_Swiss_Prot.csv"
    rcsb_file = "RCSB_by_length.csv"

    alpha_bins, alpha_frequencies = get_data_for_bars(path + alphafold_file)
    swiss_bins, swiss_frequencies = get_data_for_bars(path + swissprot_file)
    pdb_bins, pdb_frequencies = get_data_for_bars(path + pdb_file)
    bins, rcsb_frequencies = get_data_for_bars(path + rcsb_file)

    adjusted_swiss = calculate_adjusted_frequency(swiss_frequencies, alpha_frequencies)
    adjusted_rcsb = calculate_adjusted_frequency(rcsb_frequencies, pdb_frequencies)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.subplots(2, 1, sharex=True)

    ax[0].bar(bins, adjusted_swiss,
              color=_COLOUR_PALETTE["DATABANK"],
              label="Swiss-Prot Frequencies",
              bottom=alpha_frequencies)
    ax[0].bar(bins, alpha_frequencies,
              color=_COLOUR_PALETTE["USED"],
              label="Used AlphaFold Frequencies")
    ax[1].bar(bins, adjusted_rcsb,
              color=_COLOUR_PALETTE["DATABANK"],
              label="RCSB Frequencies",
              bottom=pdb_frequencies)
    ax[1].bar(bins, pdb_frequencies,
              color=_COLOUR_PALETTE["USED"],
              label="Used RCSB Frequencies")

    ax[0].tick_params(axis="x", labelrotation=90)
    ax[1].tick_params(axis="x", labelrotation=90)

    ax[0].legend(fontsize=14)
    ax[1].legend(fontsize=14)

    plt.subplots_adjust(left=0.09, bottom=0.08, top=0.99, wspace=0.05, right=1)
    plt.tight_layout()
    sns.despine()
    plt.show()


def get_comparison_data(path_to_file: str, is_simulation: bool) -> tuple:
    """
    Get csv data from path
    @param is_simulation: boolean to determine if file contains simulation data
    @param path_to_file: RCSB, AlphaFold or Simulation folder
    @return: tuple containing data required for comparison plot
    """
    dataframe = pd.read_csv(path_to_file)
    if is_simulation:
        distances = dataframe["index"].to_numpy()
    else:
        distances = dataframe["variable"].to_numpy()
    means = dataframe["mean"].to_numpy()
    normalised_means = means / np.sum(means)
    lower_bound = dataframe["lower_bound"].to_numpy()
    normalised_lower_bound = lower_bound / np.sum(means)
    upper_bound = dataframe["upper_bound"].to_numpy()
    normalised_upper_bound = upper_bound / np.sum(means)
    return distances, normalised_means, normalised_lower_bound, normalised_upper_bound


def create_comparison_plot(arguments: argparse.Namespace) -> None:
    """
    Create amino acid distance distribution plot (Fig.4)
    @param arguments: command line arguments from user
    @return: None
    """
    rcsb_tuple = get_comparison_data(f"../data/rcsb/bootstrap_{arguments.length_range}_stats.csv",
                                     is_simulation=False)
    alphafold_tuple = get_comparison_data(f"../data/alphafold/chunk_{arguments.length_range}_stats.csv",
                                          is_simulation=False)
    sim_tuple = get_comparison_data(f"../data/simulations/3d/lls_{arguments.length_range}.csv",
                                    is_simulation=True)
    rcsb_means = rcsb_tuple[1]
    alphafold_means = alphafold_tuple[1]
    rcsb_n_points = len(rcsb_means)
    alphafold_n_points = len(alphafold_means)

    rcsb_half_n_harmonic = theory_functions.harmonic_number(rcsb_n_points // 2)
    alphafold_half_n_harmonic = theory_functions.harmonic_number(alphafold_n_points // 2)
    rcsb_starting_distance = int(rcsb_means[0])
    alphafold_starting_distance = int(alphafold_means[0])
    rcsb_sum_range = np.array(range(rcsb_starting_distance, rcsb_n_points))
    alphafold_sum_range = np.array(range(alphafold_starting_distance, alphafold_n_points))
    rcsb_theory = [theory_functions.amino_acid_distance_distribution(s,
                                                                     rcsb_n_points,
                                                                     rcsb_half_n_harmonic,
                                                                     arguments.e_bs,
                                                                     arguments.d_bs) for s in rcsb_sum_range]
    alphafold_theory = [theory_functions.amino_acid_distance_distribution(s,
                                                                          alphafold_n_points,
                                                                          alphafold_half_n_harmonic,
                                                                          arguments.e_c,
                                                                          arguments.d_c) for s in alphafold_sum_range]
    plt.figure(figsize=(8, 8))
    sns.set(context="notebook", palette="colorblind", style='ticks', font_scale=1.8, font='Helvetica')
    plt.scatter(sim_tuple[0], sim_tuple[1], label="Simulation", color=_COLOUR_PALETTE["SIM_SCATTER"], marker="^")
    plt.fill_between(sim_tuple[0], sim_tuple[3], sim_tuple[2], color=_COLOUR_PALETTE["SIM_SCATTER"], alpha=0.4,
                     label="Simulation 95% C.L.", zorder=-100)
    plt.scatter(rcsb_tuple[0], rcsb_means, label=f"RCSB {arguments.length_range}", color=_COLOUR_PALETTE["PDB_SCATTER"])
    plt.scatter(alphafold_tuple[0], alphafold_means, label=f"AlphaFold {arguments.length_range}",
                color=_COLOUR_PALETTE["ALPHA_SCATTER"])
    plt.plot(rcsb_sum_range, rcsb_theory, label="Theory RCSB", color=_COLOUR_PALETTE["THEORY"], lw=1.5)
    plt.plot(alphafold_sum_range, alphafold_theory, label="Theory AlphaFold", color=_COLOUR_PALETTE["ALPHA_SCATTER"],
             lw=1.5)
    plt.yscale("log")
    plt.xscale("log")
    plt.xlim(2.5, int(alphafold_n_points / 2))
    plt.ylim(0.001, 0.2)
    plt.xlabel("s")
    plt.ylabel("P(s)")
    plt.legend()
    sns.despine()
    plt.tight_layout()
    plt.show()


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
    protein_contact_map = pcm.ProteinContactMap(pdb_file)
    alpha_carbons = protein_contact_map.get_alpha_carbons
    protein_graph = protein_contact_map.get_protein_graph(alpha_carbons)
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


def get_2d_sim_files():
    """
    Get all the files for 2D simulations
    @return: list containing all simulation files
    """
    return glob.glob(f"../data/simulations/2d/2d_sim_no_circle_*")


def get_2d_distance_frequencies(files: list) -> np.ndarray:
    """
    Get the frequencies from 2D simulation files
    @param files: 2D simulation files
    @return: array containing frequencies for amino acid distances
    """
    distances_list = []
    for file in files:
        distances = np.loadtxt(file)
        distances_list.append(distances)
    return np.asarray(distances_list)


def get_2d_simulation_stats() -> pd.DataFrame:
    """
    Compute the mean and confidence interval for 2D simulations
    @return: None
    """
    files = get_2d_sim_files()
    frequencies = get_2d_distance_frequencies(files)
    simulation_dataframe = pd.DataFrame(frequencies)
    melted_dataframe = simulation_dataframe.melt()
    simulation_stats_dataframe = melted_dataframe.groupby("variable",
                                                          as_index=False).agg(mean=("value", np.mean),
                                                                              lower_bound=("value", lambda val:
                                                                              np.quantile(val, q=0.05)),
                                                                              upper_bound=("value", lambda val:
                                                                              np.quantile(val, q=0.95)))
    simulation_stats_dataframe.to_csv(f"../data/simulations/2d/simulation_stats.csv")
    return simulation_stats_dataframe


def get_2d_plotting_data(simulation_dataframe: pd.DataFrame) -> tuple:
    """
    Get 2D simulation data from dataframe
    @param simulation_dataframe: dataframe containing 2D simulation stats
    @return: tuple containing the mean frequencies of amino acid distances, distances and the CI bounds
    """
    mean_frequencies = simulation_dataframe["mean"].to_numpy()[2:]
    distances = simulation_dataframe["variable"].to_numpy()[2:]
    lower_bound = simulation_dataframe["lower_bound"].to_numpy()[2:]
    upper_bound = simulation_dataframe["upper_bound"].to_numpy()[2:]
    return mean_frequencies, distances, lower_bound, upper_bound


