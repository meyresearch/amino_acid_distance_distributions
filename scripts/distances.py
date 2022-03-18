"""Functions for plotting RCSB and AlphaFold amino acid distance distributions"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import theory_functions
from colour_palette import _COLOUR_PALETTE
import pandas as pd
import argparse


def get_dataframe(arguments: argparse.Namespace) -> pd.DataFrame:
    """
    Get distribution of amino acid distances from csv.
    @param arguments: command line arguments
    @return: dataframe of bootstrapped/chunked data
    """
    dataframe = None
    if arguments.algorithm == "rcsb":
        dataframe = pd.read_csv(f"../data/rcsb/bootstrap_{arguments.length_range}_stats.csv")
    elif arguments.algorithm == "alpha":
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


def create_plot_label(length_range: str, algorithm: str):
    """
    Create custom labels for Alphafold PDB data and RCSB PDB data with the given link length
    @param algorithm: data source from RCSB (bootstrapped) or AlphaFold (chunked)
    @param length_range: chain length range
    @return:
    """
    plot_label = " "
    if algorithm == "rcsb":
        plot_label = f"PDB {length_range}s"
    elif algorithm == "alpha":
        plot_label = f"AF PDB {length_range}"
    return plot_label


def plot_distances(pdb_distances: np.ndarray, normalised_means: np.ndarray,
                   sim_distances: np.ndarray, normalised_sim_means: np.ndarray,
                   lower_cl: np.ndarray, upper_cl: np.ndarray,
                   plot_sum_range: np.ndarray, n_points: int, half_n_harmonic: float, exponent: int,
                   dimensionality: float, length_range: str, algorithm: str) -> None:
    """
    Plot RCSB or AlphaFold amino acid distance distributions
    @param pdb_distances: amino acid distances
    @param normalised_means: normalised frequencies of amino acid distances
    @param sim_distances: amino acid distances from 3D simulations
    @param normalised_sim_means: normalised frequencies of amino acid distances from 3D simulations
    @param lower_cl: lower bound for confidence level
    @param upper_cl: upper bound for confidence level
    @param plot_sum_range: range to sum over
    @param n_points: number of points
    @param half_n_harmonic: the "N/2"-th harmonic number
    @param exponent: exponent constant a
    @param dimensionality: dimensionality constant A
    @param length_range: chain length range
    @param algorithm: either rcsb or alpha
    @return: None
    """
    theory = [theory_functions.amino_acid_distance_distribution(s, n_points, half_n_harmonic, exponent, dimensionality)
              for s in plot_sum_range]
    plt.figure(figsize=(8, 8))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2.4, font="Helvetica")
    plt.scatter(pdb_distances, normalised_means, label=create_plot_label(length_range, algorithm),
                color=_COLOUR_PALETTE["PDB_SCATTER"])
    plt.scatter(sim_distances, normalised_sim_means, s=10, marker="^", c=_COLOUR_PALETTE["3D_SIM_SCATTER"],
                label=f"SIM {length_range}s", zorder=-50)
    plt.plot(plot_sum_range, theory, linestyle="--", label="Theory", color=_COLOUR_PALETTE["THEORY"], lw=1.5)
    plt.fill_between(pdb_distances, upper_cl, lower_cl, color=_COLOUR_PALETTE["CL"], zorder=-1, label="95% C.L.",
                     alpha=0.4)
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("s")
    plt.ylabel("P(s)")
    plt.xlim(2, n_points/2)
    plt.legend(loc="lower right")
    sns.despine()
    plt.tight_layout()


def plot_residuals(normalised_means: np.ndarray, n_points: int, half_n_harmonic: float, plot_sum_range: np.ndarray,
                   exponent: int, dimensionality: float) -> None:
    """
    Plot residuals for RCSB or AlphaFold amino acid distributions around theory
    @param normalised_means: normalised frequencies of amino acid distances
    @param n_points: number of points
    @param half_n_harmonic: the "N/2"-th harmonic number
    @param plot_sum_range: range to sum over
    @param exponent: exponent constant a
    @param dimensionality: dimensionality constant A
    @return: NoNE
    """
    plt.figure(figsize=(8, 8))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2, font="Helvetica")
    residuals = theory_functions.plotting_statistics(dimensionality, exponent, n_points, half_n_harmonic,
                                                     plot_sum_range, normalised_means)[0]
    plt.scatter(plot_sum_range, residuals, color=_COLOUR_PALETTE["RESIDUALS"], label="Residuals")
    plt.hlines(0, plot_sum_range[0], plot_sum_range[-1], color=_COLOUR_PALETTE["THEORY"], label="Theory")
    plt.ylim(-1, 1)
    plt.legend()
    plt.xlabel("s")
    plt.ylabel("Residuals")
    plt.tight_layout()


def create_plots(arguments: argparse.Namespace, exponent: int, dimensionality: float,
                 pdb_dataframe: pd.DataFrame, sim_dataframe: pd.DataFrame) -> None:
    """
    Plot amino acid distances and residuals in single plots
    @param arguments: command line arguments
    @param exponent: exponent constant a
    @param dimensionality: dimensionality constant A
    @param pdb_dataframe: dataframe containing amino acid distances and stats from RCSB or AlphaFold
    @param sim_dataframe: dataframe containing amino acid distances and stats from 3D simulations
    @return: None
    """
    pdb_plotting_tuple = get_data_for_plotting(pdb_dataframe)
    pdb_distances = pdb_dataframe["variable"].to_numpy()
    half_n_harmonic = theory_functions.harmonic_number(n_numbers=(pdb_plotting_tuple[0] // 2))
    pdb_plotting_sum_range = np.array(range(int(pdb_distances[0]), pdb_plotting_tuple[0]))
    sim_distances = sim_dataframe["index"].to_numpy()
    sim_plotting_tuple = get_data_for_plotting(sim_dataframe)
    plot_distances(pdb_distances=pdb_distances,
                   normalised_means=pdb_plotting_tuple[1],
                   sim_distances=sim_distances,
                   normalised_sim_means=sim_plotting_tuple[1],
                   lower_cl=pdb_plotting_tuple[2],
                   upper_cl=pdb_plotting_tuple[3],
                   plot_sum_range=pdb_plotting_sum_range,
                   n_points=pdb_plotting_tuple[0],
                   half_n_harmonic=half_n_harmonic,
                   exponent=exponent,
                   dimensionality=dimensionality,
                   length_range=arguments.length_range,
                   algorithm=arguments.algorithm)
    plot_residuals(normalised_means=pdb_plotting_tuple[1],
                   n_points=pdb_plotting_tuple[0],
                   half_n_harmonic=half_n_harmonic,
                   plot_sum_range=pdb_plotting_sum_range,
                   exponent=exponent,
                   dimensionality=dimensionality)
    plt.show()



def grid_plot_distances(dimensionality_range: np.ndarray, exponent_range: np.ndarray, n_points: int,
                        half_n_harmonic: float, plotting_sumrange: np.ndarray, normalised_means: np.ndarray,
                        distances: np.ndarray, normalised_sim_means: np.ndarray, sim_distances: np.ndarray,
                        lower_cl: np.ndarray, upper_cl: np.ndarray, length_range: str, algorithm: str) -> None:
    """
    Plot amino acid distances and theory plot with different values of exponent and dimensionality scaling constant
    @param algorithm: data source from RCSB (bootstrapped) or AlphaFold (chunked)
    @param length_range: chain length range
    @param sim_distances: amino acid distances for 3D simulated data (x-axis)
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
                                 c=_COLOUR_PALETTE["3D_SIM_SCATTER"], label=f"SIM {length_range}s", zorder=-50)
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
