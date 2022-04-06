"""Functions for plotting RCSB and AlphaFold amino acid distance distributions"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import seaborn as sns
import theory_functions
from colour_palette import _COLOUR_PALETTE
import pandas as pd
import argparse


def get_histogram(arguments: argparse.Namespace, algorithm: str) -> np.ndarray:
    """
    Get distribution of amino acid distances from histogram.
    @param arguments: command line arguments
    @return: histogram of RCSB/AlphaFold 2 data
    """
    histogram = None
    if algorithm == "rcsb":
        histogram = np.load(f"../data/rcsb/histogram_{arguments.length_range}.npy", allow_pickle=True)
    elif algorithm == "alpha":
        histogram = np.load(f"../data/alphafold/histogram_{arguments.length_range}.npy", allow_pickle=True)
    return histogram


def get_simulation_histogram(length_range: str) -> np.ndarray:
    """
    Get distribution of amino acid distances from histogram.
    @param length_range: chain length
    @return: histogram of 3d simulation data
    """
    return np.load(f"../data/simulations/3d/histogram_{length_range}.npy", allow_pickle=True)


def get_confidence_interval(histogram: np.ndarray, quantile: int) -> tuple:
    """
    Get the confidence interval chosen by the user
    @param histogram: numpy array that contains amino acid distances
    @param quantile: 1, 2 or 3 standard deviations, given by user
    @return: tuple containing lower and upper bounds of the confidence interval
    """
    lower_bound, upper_bound = [], []
    if quantile == 1:
        lower_bound, upper_bound = np.quantile(histogram, q=[1 / 3, 2 / 3], axis=0)
    elif quantile == 2:
        lower_bound, upper_bound = np.quantile(histogram, q=[0.05, 0.95], axis=0)
    elif quantile == 3:
        lower_bound, upper_bound = np.quantile(histogram, q=[0.003, 0.997], axis=0)
    return lower_bound, upper_bound


def get_measure_of_central_tendency(histogram: np.ndarray, user_measure: str) -> np.ndarray:
    """
    Return the desired measure of central tendency, median or mean
    @param histogram: numpy array that contains amino acid distances
    @param user_measure: user specified measure of central tendency
    @return: desired measure
    """
    measure = []
    if user_measure == "median":
        measure = np.quantile(histogram, q=0.5, axis=0)
    elif user_measure == "mean":
        measure = np.mean(histogram, axis=0)
    return measure


def filter_histogram(histogram: np.ndarray, start_point: int, end_point: int) -> np.ndarray:
    """
    Use
    @param histogram: numpy array that contains amino acid distances
    @param start_point: point from which to start plotting
    @param end_point: point at which to stop plotting
    @return: numpy array of the filtered histogram range
    """
    rows = histogram.shape[0]
    cols = histogram.shape[1]
    histogram_list = []
    for i in range(rows):
        for j in range(cols):
            filtered_histogram = histogram[i][start_point:end_point]
            histogram_list.append(filtered_histogram)
    return np.asarray(histogram_list)


def get_data_for_plotting(filtered_histogram: np.ndarray, arguments: argparse.Namespace) -> tuple:
    """
    Get protein data for plotting from dataframe
    @param arguments: command line arguments
    @param filtered_histogram: numpy array that contains amino acid distances
    @return: tuple of number of datapoints, normalised means and confidence level bounds
    """
    distances = np.linspace(start=arguments.start_point, stop=arguments.end_point, num=len(filtered_histogram[0]))
    lower_bound, upper_bound = get_confidence_interval(filtered_histogram, arguments.quantile)
    measure = get_measure_of_central_tendency(filtered_histogram, arguments.measure)
    normalised_measure = measure / np.sum(measure)
    normalised_lower_bound = lower_bound / np.sum(measure)
    normalised_upper_bound = upper_bound / np.sum(measure)
    # n_datapoints = len(distances)
    n_datapoints = int(arguments.length_range)
    return n_datapoints, distances, normalised_measure, normalised_lower_bound, normalised_upper_bound


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


def plot_distances(pdb_distances: np.ndarray, normalised_measure: np.ndarray,
                   sim_distances: np.ndarray, normalised_sim_measure: np.ndarray,
                   lower_cl: np.ndarray, upper_cl: np.ndarray,
                   plot_sum_range: np.ndarray, n_points: int, half_n_harmonic: float, exponent: int,
                   dimensionality: float, length_range: str, algorithm: str) -> None:
    """
    Plot RCSB or AlphaFold amino acid distance distributions
    @param pdb_distances: amino acid distances
    @param normalised_measure: normalised frequencies of amino acid distances
    @param sim_distances: amino acid distances from 3D simulations
    @param normalised_sim_measure: normalised frequencies of amino acid distances from 3D simulations
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
    plt.scatter(pdb_distances, normalised_measure, label=create_plot_label(length_range, algorithm),
                color=_COLOUR_PALETTE["PDB_SCATTER"])
    plt.scatter(sim_distances, normalised_sim_measure, s=10, marker="^", c=_COLOUR_PALETTE["3D_SIM_SCATTER"],
                label=f"SIM {length_range}s", zorder=-50)
    plt.plot(plot_sum_range, theory, linestyle="--", label="Theory", color=_COLOUR_PALETTE["THEORY"], lw=1.5)
    plt.fill_between(pdb_distances, upper_cl, lower_cl, color=_COLOUR_PALETTE["CL"], zorder=-1, label="95% C.L.",
                     alpha=0.4)
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("s")
    plt.ylabel("P(s)")
    plt.xlim(2, n_points / 2)
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
                 pdb_histogram: np.ndarray, sim_histogram: np.ndarray) -> None:
    """
    Plot amino acid distances and residuals in single plots
    @param arguments: command line arguments
    @param exponent: exponent constant a
    @param dimensionality: dimensionality constant A
    @param pdb_histogram: histogram containing amino acid distances from RCSB or AlphaFold
    @param sim_histogram: dataframe containing amino acid distances and stats from 3D simulations
    @return: None
    """
    pdb_plotting_tuple = get_data_for_plotting(pdb_histogram, arguments)
    pdb_distances = pdb_plotting_tuple[1]
    half_n_harmonic = theory_functions.harmonic_number(n_numbers=(pdb_plotting_tuple[0] // 2))
    pdb_plotting_sum_range = np.array(range(int(pdb_distances[0]), pdb_plotting_tuple[0]))
    sim_plotting_tuple = get_data_for_plotting(sim_histogram, arguments)
    sim_distances = sim_plotting_tuple[1]
    plot_distances(pdb_distances=pdb_distances, normalised_measure=pdb_plotting_tuple[2], sim_distances=sim_distances,
                   normalised_sim_measure=sim_plotting_tuple[2], lower_cl=pdb_plotting_tuple[3],
                   upper_cl=pdb_plotting_tuple[4], plot_sum_range=pdb_plotting_sum_range,
                   n_points=pdb_plotting_tuple[0], half_n_harmonic=half_n_harmonic, exponent=exponent,
                   dimensionality=dimensionality, length_range=arguments.length_range, algorithm=arguments.algorithm)
    plot_residuals(normalised_means=pdb_plotting_tuple[2],
                   n_points=pdb_plotting_tuple[0],
                   half_n_harmonic=half_n_harmonic,
                   plot_sum_range=pdb_plotting_sum_range,
                   exponent=exponent,
                   dimensionality=dimensionality)
    plt.show()


def grid_plot_distances(dimensionality_range: np.ndarray, exponent_range: np.ndarray, n_points: int,
                        half_n_harmonic: float, plotting_sumrange: np.ndarray, normalised_measure: np.ndarray,
                        distances: np.ndarray, normalised_sim_measure: np.ndarray, sim_distances: np.ndarray,
                        lower_cl: np.ndarray, upper_cl: np.ndarray, length_range: str, algorithm: str,
                        start_point: int, end_point: int) -> None:
    """
    Plot amino acid distances and theory plot with different values of exponent and dimensionality scaling constant
    @param normalised_measure:
    @param normalised_sim_measure:
    @param start_point:
    @param end_point:
    @param algorithm: data source from RCSB (bootstrapped) or AlphaFold (chunked)
    @param length_range: chain length range
    @param sim_distances: amino acid distances for 3D simulated data (x-axis)
    @param upper_cl: upper confidence level
    @param lower_cl: lower confidence level
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
            theory = theory_functions.amino_acid_distance_distribution(plotting_sumrange,
                                                                       n_points,
                                                                       half_n_harmonic,
                                                                       exponent_range[row],
                                                                       dimensionality_range[col])
            ax[row][col].scatter(distances, normalised_measure, s=10, c=_COLOUR_PALETTE["PDB_SCATTER"],
                                 label=plot_label)
            ax[row][col].fill_between(distances, upper_cl, lower_cl, color=_COLOUR_PALETTE["PDB_SCATTER"], alpha=0.4,
                                      label="95% CL", zorder=-100)
            ax[row][col].scatter(sim_distances, normalised_sim_measure, s=10, marker="^",
                                 c=_COLOUR_PALETTE["3D_SIM_SCATTER"], label=f"SIM {length_range}s", zorder=-50)
            ax[row][col].plot(plotting_sumrange, theory, c=_COLOUR_PALETTE["THEORY"], label="Theory")
            ax[row][col].set_yscale("log")
            ax[row][col].set_xscale("log")
            ax[row][col].set_ylim(0.0001, 1)
            ax[row][col].set_xlim(start_point, end_point)
            ax[row][-1].set_ylabel(f"a = {exponent_range[row]}", fontsize=13, rotation=0, labelpad=21)
            ax[row][-1].yaxis.set_label_position("right")
            ax[0][col].set_title(f"A = {dimensionality_range[col]:2f}")
    plt.legend(bbox_to_anchor=[1.57, 4.65], fontsize=9)
    fig.text(0.5, 0.025, "s", ha="center", fontsize=15.5)
    fig.text(0.005, 0.5, "P(s)", va="center", rotation="vertical", fontsize=15.5)
    plt.subplots_adjust(left=0.06, bottom=0.08, top=0.95, wspace=0.05, right=0.909)
    # plt.savefig(f"../plots/supplementary_information/{algorithm}_{length_range}.pdf")


def fit_distances(dimensionality_range: np.ndarray, exponent_range: np.ndarray, n_points: int,
                  half_n_harmonic: float, plotting_sumrange: np.ndarray, normalised_measure: np.ndarray,
                  distances: np.ndarray, normalised_sim_measure: np.ndarray, sim_distances: np.ndarray,
                  lower_cl: np.ndarray, upper_cl: np.ndarray, length_range: str, algorithm: str,
                  start_point: int, end_point: int) -> None:
    """
    Plot amino acid distances and theory plot with different values of exponent and dimensionality scaling constant
    @param algorithm: data source from RCSB (bootstrapped) or AlphaFold (chunked)
    @param length_range: chain length range
    @param sim_distances: amino acid distances for 3D simulated data (x-axis)
    @param upper_cl: upper confidence level
    @param lower_cl: lower confidence level
    @param normalised_sim_measure: normalised simulated means of amino acid distance frequencies
    @param normalised_measure: normalised means or medians of PDB amino acid distance frequencies
    @param distances: amino acid distances for PDB data (x-axis)
    @param plotting_sumrange: Range to be used in sum for the theory plot
    @param half_n_harmonic: Harmonic number for N/2
    @param n_points: number of data points
    @param exponent_range: range of exponent values (constant a)
    @type dimensionality_range: range of dimensionality scaling constant (constant A)
    @param start_point: point from which to start plotting
    @param end_point: point at which to end plotting
    @return: None
    """
    plot_label = create_plot_label(length_range, algorithm)
    plt.figure(figsize=(8, 8))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=1.8, font="Helvetica")
    theory_parameters, covariance = scipy.optimize.curve_fit(
        f=theory_functions.amino_acid_distance_distribution, xdata=plotting_sumrange, ydata=normalised_measure,
        bounds=([n_points - 1, half_n_harmonic - 0.00001, exponent_range[0], dimensionality_range[0]],
                [n_points, half_n_harmonic, exponent_range[-1], dimensionality_range[-1]]))
    print("-----Theory parameters-----------")
    print(f"N: {theory_parameters[0]}, H-N/2: {theory_parameters[1]}, a: {theory_parameters[2]}, A: {theory_parameters[3]}")
    power_law_parameters, covariance = scipy.optimize.curve_fit(f=theory_functions.power_law,
                                                                xdata=plotting_sumrange,
                                                                ydata=normalised_measure)
    print("-----Power law parameters--------")
    print(f"gamma: {power_law_parameters[0]}, C: {power_law_parameters[1]}")
    plt.plot(plotting_sumrange, theory_functions.power_law(plotting_sumrange, *power_law_parameters),
             c=_COLOUR_PALETTE["POWER"], label="Power law", lw=1.5)
    plt.scatter(distances, normalised_measure, c=_COLOUR_PALETTE["PDB_SCATTER"],
                label=plot_label)
    plt.scatter(sim_distances, normalised_sim_measure, marker="^", c=_COLOUR_PALETTE["3D_SIM_SCATTER"],
                label=f"SIM {length_range}s", zorder=-50)
    plt.fill_between(distances, upper_cl, lower_cl, color=_COLOUR_PALETTE["PDB_SCATTER"], alpha=0.4, label="95% CL",
                     zorder=-100)
    plt.plot(plotting_sumrange,
             theory_functions.amino_acid_distance_distribution(plotting_sumrange, *theory_parameters),
             c=_COLOUR_PALETTE["THEORY"], label="Theory", lw=1.5)
    plt.loglog()
    plt.ylim(0.001, 1)
    plt.xlim(start_point, int(length_range) // 2)
    plt.legend()
    plt.xlabel("s")
    plt.ylabel("P(s)")
    sns.despine()


def grid_plot_residuals(algorithm: str, length_range: str, dimensionality_range, exponent_range: np.ndarray,
                        n_points: int, half_n_harmonic_number: float, plotting_sumrange: np.ndarray,
                        normalised_measure: np.ndarray, start_point: int, end_point: int) -> None:
    """
    Plot residuals grid and print out statistics.
    @param algorithm: rcsb or alpha
    @param length_range: chain length
    @param dimensionality_range: range of values for dimensionality constant A
    @param exponent_range: range of values for exponent constant a
    @param n_points: number of datapoints
    @param half_n_harmonic_number: the "N/2-th" harmonic number
    @param plotting_sumrange: range to sum over
    @param normalised_measure: normalised median or mean
    @param start_point: point from which to start plotting
    @param end_point: point at which to end plotting
    @return:
    """
    fig = plt.figure(figsize=(16, 8))
    sns.set(context="notebook", palette="colorblind", style='ticks', font='Helvetica')
    ax = fig.subplots(len(exponent_range), len(dimensionality_range), sharex=True, sharey=True)
    for row in range(len(exponent_range)):
        for col in range(len(dimensionality_range)):
            residuals = theory_functions.plotting_statistics(dimensionality_range[col], exponent_range[row], n_points,
                                                             half_n_harmonic_number, plotting_sumrange,
                                                             normalised_measure)[0]

            ax[row][col].scatter(plotting_sumrange, residuals, marker=".", c=_COLOUR_PALETTE["RESIDUALS"],
                                 label="Residuals", zorder=10)
            ax[row][col].hlines(0, plotting_sumrange[0], plotting_sumrange[-1], color=_COLOUR_PALETTE["THEORY"],
                                label="Theory")
            ax[row][-1].set_ylabel(f"a = {exponent_range[row]}", fontsize=13, rotation=0, labelpad=21)
            ax[row][-1].yaxis.set_label_position("right")
            ax[row][col].set_ylim(-0.025, 0.025)
            ax[row][col].set_xlim(start_point, end_point)
            ax[0][col].set_title(f"A = {dimensionality_range[col]:2f}")
            ax[row][col].legend(fontsize=9)

    fig.text(0.5, 0.025, "s", ha="center", fontsize=15.5)
    fig.text(0.005, 0.5, "Residuals", va="center", rotation="vertical", fontsize=15.5)
    plt.subplots_adjust(left=0.06, bottom=0.08, top=0.95, wspace=0.1, right=0.95)
    plt.savefig(f"../plots/supplementary_information/{algorithm}_{length_range}_r.pdf")


def create_grid_plots(arguments: argparse.Namespace, pdb_histogram: np.ndarray, sim_histogram: np.ndarray) -> None:
    """
    Function to bring together everything needed to plot the grid plots of amino acid distance distributions
    and residuals
    @param arguments: command line arguments
    @param pdb_histogram: histogram containing data for PDBs
    @param sim_histogram: 3D simulation data for plotting
    @return:
    """
    filtered_pdb_histogram = filter_histogram(pdb_histogram, arguments.start_point, arguments.end_point)
    pdb_plotting_tuple = get_data_for_plotting(filtered_pdb_histogram, arguments)
    pdb_distances = pdb_plotting_tuple[1]
    half_n_harmonic = theory_functions.harmonic_number(n_numbers=(pdb_plotting_tuple[0] // 2))
    # pdb_plotting_sum_range = np.linspace(start=1, stop=int(arguments.length_range), num=100)[4:50]
    pdb_plotting_sum_range = pdb_distances
    filtered_sim_histogram = filter_histogram(sim_histogram, arguments.start_point, arguments.end_point)
    sim_plotting_tuple = get_data_for_plotting(filtered_sim_histogram, arguments)
    sim_distances = sim_plotting_tuple[1]
    dimensionality_range = np.arange(arguments.start_dimensionality,
                                     arguments.end_dimensionality,
                                     arguments.step_dimensionality)
    exponent_range = np.arange(arguments.start_exponent,
                               arguments.end_exponent,
                               arguments.step_exponent)
    grid_plot_distances(dimensionality_range=dimensionality_range, exponent_range=exponent_range,
                        n_points=pdb_plotting_tuple[0], half_n_harmonic=half_n_harmonic,
                        plotting_sumrange=pdb_plotting_sum_range, normalised_measure=pdb_plotting_tuple[2],
                        distances=pdb_distances, normalised_sim_measure=sim_plotting_tuple[2],
                        sim_distances=sim_distances, lower_cl=pdb_plotting_tuple[3], upper_cl=pdb_plotting_tuple[4],
                        length_range=arguments.length_range, algorithm=arguments.algorithm,
                        start_point=arguments.start_point,
                        end_point=arguments.end_point)
    fit_distances(dimensionality_range=dimensionality_range, exponent_range=exponent_range,
                  n_points=pdb_plotting_tuple[0], half_n_harmonic=half_n_harmonic,
                  plotting_sumrange=pdb_plotting_sum_range, normalised_measure=pdb_plotting_tuple[2],
                  distances=pdb_distances, normalised_sim_measure=sim_plotting_tuple[2],
                  sim_distances=sim_distances,
                  lower_cl=pdb_plotting_tuple[3], upper_cl=pdb_plotting_tuple[4],
                  length_range=arguments.length_range, algorithm=arguments.algorithm,
                  start_point=arguments.start_point, end_point=arguments.end_point)
    # grid_plot_residuals(algorithm=arguments.algorithm, length_range=arguments.length_range,
    #                     dimensionality_range=dimensionality_range, exponent_range=exponent_range,
    #                     n_points=pdb_plotting_tuple[0], half_n_harmonic_number=half_n_harmonic,
    #                     plotting_sumrange=pdb_plotting_sum_range, normalised_measure=pdb_plotting_tuple[2],
    #                     start_point=arguments.start_point, end_point=arguments.end_point)
    plt.show()
