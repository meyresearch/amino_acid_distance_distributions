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


def filter_histogram(histogram: np.ndarray, arguments: argparse.Namespace) -> np.ndarray:
    """
    Use

    @param histogram: numpy array that contains amino acid distances
    @param arguments: command line arguments
    @return: numpy array of the filtered histogram range
    """
    rows = histogram.shape[0]
    cols = histogram.shape[1]
    half = len(histogram[0]) // 2
    all_bins = np.linspace(start=1, stop=int(arguments.length_range), num=len(histogram[0]))
    start_index = 0
    for i in range(len(all_bins)):
        if int(all_bins[i]) == arguments.start_point:
            start_index = i
    histogram_list = []
    for i in range(rows):
        for j in range(cols):
            filtered_histogram = histogram[i][start_index:half]
            histogram_list.append(filtered_histogram)
    return np.asarray(histogram_list)


def get_data_for_plotting(filtered_histogram: np.ndarray, arguments: argparse.Namespace) -> tuple:
    """
    Get protein data for plotting from dataframe
    @param arguments: command line arguments
    @param filtered_histogram: numpy array that contains amino acid distances
    @return: tuple of number of datapoints, normalised means and confidence level bounds
    """
    distances = np.linspace(start=arguments.start_point, stop=arguments.end_point,
                            num=len(filtered_histogram[0]))
    lower_bound, upper_bound = get_confidence_interval(filtered_histogram, arguments.quantile)
    measure = get_measure_of_central_tendency(filtered_histogram, arguments.measure)
    normalised_measure = measure / np.sum(measure)
    normalised_lower_bound = lower_bound / np.sum(measure)
    normalised_upper_bound = upper_bound / np.sum(measure)
    # chain_length = len(distances)
    chain_length = int(arguments.length_range)
    return chain_length, distances, normalised_measure, normalised_lower_bound, normalised_upper_bound


def create_comparison_plot(arguments: argparse.Namespace, rcsb_histogram: np.ndarray,
                           alpha_histogram: np.ndarray, sim_histogram: np.ndarray) -> None:
    """
    Create amino acid distance distribution plot (Fig.4)
    @param arguments: command line arguments from user
    @param rcsb_histogram: numpy array containing RCSB PDB data
    @param alpha_histogram: numpy array containing AlphaFold 2 PDB data
    @param sim_histogram: numpy array containing 3d simulation data
    @return: None
    """
    filtered_rcsb_histogram = filter_histogram(rcsb_histogram, arguments)
    rcsb_plotting_tuple = get_data_for_plotting(filtered_rcsb_histogram, arguments)
    rcsb_chain_length = rcsb_plotting_tuple[0]
    rcsb_distances = rcsb_plotting_tuple[1]
    rcsb_half_n_harmonic = theory_functions.harmonic_number(n_numbers=(rcsb_chain_length // 2))
    # rcsb_plotting_sum_range = np.linspace(start=11,
    #                                       stop=70,
    #                                       num=len(filtered_rcsb_histogram[0]))
    rcsb_plotting_sum_range = np.linspace(start=arguments.start_point,
                                          stop=int(arguments.length_range) ,
                                          num=len(filtered_rcsb_histogram[0]))

    filtered_alpha_histogram = filter_histogram(alpha_histogram, arguments)
    alpha_plotting_tuple = get_data_for_plotting(filtered_alpha_histogram, arguments)
    alpha_chain_length = alpha_plotting_tuple[0]
    alpha_distances = alpha_plotting_tuple[1]
    alpha_half_n_harmonic = theory_functions.harmonic_number(n_numbers=(alpha_chain_length // 2))
    # alpha_plotting_sum_range = np.linspace(start=arguments.start_point,
    #                                        stop=int(arguments.length_range)//2,
    #                                        num=len(filtered_rcsb_histogram[0]))

    filtered_sim_histogram = filter_histogram(sim_histogram, arguments)
    sim_plotting_tuple = get_data_for_plotting(filtered_sim_histogram, arguments)
    sim_distances = sim_plotting_tuple[1]

    rcsb_dimensionality_range = np.arange(arguments.rcsb_startd,
                                          arguments.rcsb_endd,
                                          arguments.step_dimensionality)
    rcsb_exponent_range = np.arange(arguments.rcsb_starte,
                                    arguments.rcsb_ende,
                                    arguments.step_exponent)
    alpha_dimensionality_range = np.arange(arguments.alpha_startd,
                                           arguments.alpha_endd,
                                           arguments.step_dimensionality)
    alpha_exponent_range = np.arange(arguments.alpha_starte,
                                     arguments.alpha_ende,
                                     arguments.step_exponent)

    rcsb_parameters, rcsb_cov = scipy.optimize.curve_fit(f=theory_functions.amino_acid_distance_distribution,
                                                         xdata=rcsb_plotting_sum_range, ydata=rcsb_plotting_tuple[2],
                                                         sigma=rcsb_plotting_tuple[4] - rcsb_plotting_tuple[3],
                                                         bounds=([rcsb_chain_length - 1, rcsb_half_n_harmonic - 0.00001,
                                                                  rcsb_exponent_range[0], rcsb_dimensionality_range[0]],
                                                                 [rcsb_chain_length, rcsb_half_n_harmonic,
                                                                  rcsb_exponent_range[-1],
                                                                  rcsb_dimensionality_range[-1]]))

    rcsb_power, rcsb_power_cov = scipy.optimize.curve_fit(f=theory_functions.power_law,
                                                          xdata=rcsb_plotting_sum_range, ydata=rcsb_plotting_tuple[2])
    rcsb_theory = theory_functions.amino_acid_distance_distribution(rcsb_plotting_sum_range, *rcsb_parameters)
    rcsb_powerlaw = theory_functions.power_law(rcsb_plotting_sum_range, *rcsb_power)

    fig = plt.figure()
    fig.set_size_inches((8, 8))
    sns.set(context="notebook", palette="colorblind", style='ticks', font_scale=1.8, font="Helvetica")

    plt.scatter(sim_distances, sim_plotting_tuple[2], label="Simulation", color=_COLOUR_PALETTE["3D_SIM_SCATTER"],
                marker="^")
    # plt.fill_between(sim_distances, sim_plotting_tuple[4], sim_plotting_tuple[3],
    #                  color=_COLOUR_PALETTE["3D_SIM_SCATTER"], alpha=0.4, label="Simulation 95% C.L.", zorder=-100)

    plt.scatter(rcsb_distances, rcsb_plotting_tuple[2], label=f"RCSB {arguments.length_range}",
                color=_COLOUR_PALETTE["PDB_SCATTER"], marker="o")
    plt.fill_between(rcsb_distances, rcsb_plotting_tuple[4], rcsb_plotting_tuple[3],
                     color=_COLOUR_PALETTE["PDB_SCATTER"], alpha=0.4, label="RCSB 95% C.L.", zorder=-99)
    plt.scatter(alpha_distances, alpha_plotting_tuple[2], label=f"AlphaFold 2 {arguments.length_range}",
                color=_COLOUR_PALETTE["ALPHA_SCATTER"], marker="s")

    plt.plot(rcsb_plotting_sum_range, rcsb_theory, label="Theory RCSB", color="k", lw=1.5)
    plt.plot(rcsb_plotting_sum_range, rcsb_powerlaw, label="Power law RCSB", color="k", lw=1.5, ls="--")

    print("-----------------------RCSB THEORY-----------------------")
    print(f"N: {rcsb_parameters[0]}")
    print(f"H-N/2: {rcsb_parameters[1]}")
    print(f"a: {rcsb_parameters[2]}")
    print(f"A: {rcsb_parameters[3]}")
    print("----------------------RCSB POWER LAW----------------------")
    print(f"gamma: {rcsb_power[0]}")
    print(f"constant: {rcsb_power[1]}")

    plt.yscale("log")
    plt.xscale("log")
    # plt.xlim(arguments.start_point, int(arguments.length_range)//2)
    plt.xlim(1, int(arguments.length_range))
    # plt.ylim(0.001, 0.2)
    plt.xlabel("s")
    plt.ylabel("P(s)")
    plt.legend()
    sns.despine()
    plt.tight_layout()
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
    # plot_label = create_plot_label(length_range, algorithm)
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
            ax[row][col].scatter(distances, normalised_measure, s=10, c=_COLOUR_PALETTE["PDB_SCATTER"])
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


def create_grid_plots(arguments: argparse.Namespace, pdb_histogram: np.ndarray, sim_histogram: np.ndarray) -> None:
    """
    Function to bring together everything needed to plot the grid plots of amino acid distance distributions
    and residuals
    @param arguments: command line arguments
    @param pdb_histogram: histogram containing data for PDBs
    @param sim_histogram: 3D simulation data for plotting
    @return:
    """
    filtered_pdb_histogram = filter_histogram(pdb_histogram, arguments)
    pdb_plotting_tuple = get_data_for_plotting(filtered_pdb_histogram, arguments)
    pdb_distances = pdb_plotting_tuple[1]
    half_n_harmonic = theory_functions.harmonic_number(n_numbers=(pdb_plotting_tuple[0] // 2))
    # pdb_plotting_sum_range = np.linspace(start=1, stop=int(arguments.length_range), num=100)[4:50]
    pdb_plotting_sum_range = np.linspace(start=arguments.start_point,
                                         stop=int(arguments.length_range) // 2,
                                         num=len(filtered_rcsb_histogram[0]))
    filtered_sim_histogram = filter_histogram(sim_histogram, arguments)
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
    plt.show()
