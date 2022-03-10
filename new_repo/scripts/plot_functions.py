"""
Functions used for plotting amino acid distance distributions.
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import theory_functions
import seaborn as sns

_COLOUR_PALETTE = {"PDB_SCATTER": "#006374",
                   "SIM_SCATTER": "#fbafe4",
                   "ALPHA_SCATTER": "#1ac938",
                   "THEORY": "#006374",
                   "RESIDUALS": "#fbafe4",
                   "DATABANK": "#006374",
                   "USED": "#fbafe4"}


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
    parser.add_argument("algorithm", type=str, choices=["BS", "C", "BOTH", "B"],
                        help="get distances from bootstrapping (BS), chunking (C), compare both or plot bars (B)")
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
                             arguments.d_bs, arguments.e_bs, arguments.d_c, arguments.e_c)
    return arguments


def check_required_arguments(argument_parser: argparse.ArgumentParser, given_algorithm: str, given_length: str,
                             starting_d: str, ending_d: str, starting_e: str, ending_e: str, bootstrap_d: str,
                             bootstrap_e: str, chunk_d: str, chunk_e: str) -> None:
    """
    Check given CL arguments so that all requirements are met
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
    @param arguments:
    @param pdb_dataframe:
    @param sim_dataframe:
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

    alpha_bins, alpha_frequencies = get_data_for_bars(path+alphafold_file)
    swiss_bins, swiss_frequencies = get_data_for_bars(path+swissprot_file)
    pdb_bins, pdb_frequencies = get_data_for_bars(path+pdb_file)
    rcsb_bins, rcsb_frequencies = get_data_for_bars(path+rcsb_file)

    adjusted_swiss = calculate_adjusted_frequency(swiss_frequencies, alpha_frequencies)
    adjusted_rcsb = calculate_adjusted_frequency(rcsb_frequencies, pdb_frequencies)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.subplots(2, 1, sharex=True)

    ax[0].bar(swiss_bins, adjusted_swiss,
              color=_COLOUR_PALETTE["DATABANK"],
              label="Swiss-Prot Frequencies",
              bottom=alpha_frequencies)
    ax[0].bar(alpha_bins, alpha_frequencies,
              color=_COLOUR_PALETTE["USED"],
              label="Used AlphaFold Frequencies")
    ax[1].bar(rcsb_bins, adjusted_swiss,
              color=_COLOUR_PALETTE["DATABANK"],
              label="RCSB Frequencies",
              bottom=pdb_frequencies)
    ax[1].bar(pdb_bins, pdb_frequencies,
              color=_COLOUR_PALETTE["USED"],
              label="Used RCSB Frequencies")

    ax[0].tick_params(axis="x", labelrotation=90)
    ax[1].tick_params(axis="x", labelrotation=90)

    ax[0].legend()
    ax[1].legend()

    fig.text(0.5, 0.027, "Chain length", ha="center")
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















