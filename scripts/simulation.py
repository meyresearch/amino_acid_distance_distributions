"""Functions for plotting 2D simulation amino acid distances."""
import argparse
import glob
import numpy as np
import pandas as pd
import scipy.optimize
import random
import theory_functions
import matplotlib.pyplot as plt
import seaborn as sns
from colour_palette import _COLOUR_PALETTE


def get_2d_files():
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


def get_2d_stats() -> pd.DataFrame:
    """
    Compute the mean and confidence interval for 2D simulations
    @return: None
    """
    files = get_2d_files()
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


def plot_2d_distances(distances: np.ndarray, simulation_means: np.ndarray, lower_cl: np.ndarray,
                      upper_cl: np.ndarray, chain_length: int,
                      half_n_harmonic: float, exponent_range: np.ndarray, dimensionality_range: np.ndarray,
                      start_point: int, end_point: int) -> None:
    """
    Plot 2D simulation distance distribution and compare to theory
    @param distances: amino acid distances in chain length range
    @param simulation_means: mean frequencies of amino acid distances
    @param lower_cl: lower bound for confidence interval
    @param upper_cl: upper bound for confidence interval
    @param chain_length: length of chain
    @param half_n_harmonic: the "N/2"-th harmonic number
    @param exponent_range: constant a
    @param dimensionality_range: constant A
    @param start_point: point from which to start plotting
    @param end_point: point at which to stop plotting
    @return: None
    """
    random.seed(1234)
    exponent_index = random.randint(0, len(exponent_range) - 1)
    dimensionality_index = random.randint(0, len(dimensionality_range) - 1)
    print(f"exponent index: {exponent_index}, dimensionality_index: {dimensionality_index}")
    # theory_x = np.linspace(2, chain_length//2, chain_length//2)
    parameters = [chain_length//2, half_n_harmonic,
                  exponent_range[exponent_index], dimensionality_range[dimensionality_index]]
    least_squares_result = scipy.optimize.least_squares(theory_functions.vector_of_residuals, parameters,
                                                        args=(distances, simulation_means))
    optimised_parameters = least_squares_result.x
    jacobian = least_squares_result.jac
    covariance_matrix = theory_functions.covariance_matrix(jacobian)
    sigma = theory_functions.error_on_least_squares(covariance_matrix)
    plt.plot(distances, theory_functions.amino_acid_distance_distribution(distances, *optimised_parameters))
    plt.scatter(distances, simulation_means)
    plt.loglog()
    plt.show()

    # print(distances.size)
    # parameters, covariance = scipy.optimize.curve_fit(f=theory_functions.amino_acid_distance_distribution,
    #                                                   xdata=distances, ydata=simulation_means,
    #                                                   p0=[chain_length, half_n_harmonic,
    #                                                       exponent_range[exponent_index],
    #                                                       dimensionality_range[dimensionality_index]],
    #                                                   sigma=upper_cl-lower_cl, absolute_sigma=True)
    # bounds=([0, 0,
    #          exponent_range[0], dimensionality_range[0]],
    #         [chain_length, half_n_harmonic,
    #          exponent_range[-1], dimensionality_range[-1]]))
    # theory = theory_functions.amino_acid_distance_distribution(distances, *parameters)
    # sigma = np.sqrt(np.diag(covariance))
    # plt.figure(figsize=(8, 8))
    # sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2.4, font="Helvetica")
    # plt.scatter(distances, simulation_means, label="2D Simulation", color=_COLOUR_PALETTE["2D_SIM_SCATTER"])
    # plt.plot(distances, theory, linestyle="--", label="Theory", color=_COLOUR_PALETTE["THEORY"], lw=1.5)
    # plt.fill_between(distances, upper_cl, lower_cl, color=_COLOUR_PALETTE["CL"], zorder=-10, label="95% C.L.",
    #                  alpha=0.4)
    #
    # print("----------------Parameters----------------")
    # print(f"Chain length: {parameters[0]}")
    # print(f"N/2 Harmonic: {parameters[1]}")
    # print(f"Exponent: {parameters[2]}")
    # print(f"Dimensionality: {parameters[3]}")
    # print("----------------Covariance----------------")
    # print(covariance)
    # print("-------------------Sigma------------------")
    # print(f"Chain length: {sigma[0]}")
    # print(f"N/2 Harmonic: {sigma[1]}")
    # print(f"Exponent: {sigma[2]}")
    # print(f"Dimensionality: {sigma[3]}")
    #
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.xlabel("s")
    # plt.ylabel("P(s)")
    # plt.xlim(start_point, end_point)
    # plt.legend(loc="upper right")
    # sns.despine()
    # plt.tight_layout()


def plot_2d_residuals(simulation_means: np.ndarray, n_points: int, half_n_harmonic: float,
                      plot_sum_range: np.ndarray, exponent: int, dimensionality: float) -> None:
    """
    Plot residuals for 2D simulations about the theory
    @param simulation_means: mean frequencies of amino acid distances for 2D simulation
    @param n_points: number of points
    @param half_n_harmonic: the "N/2"-th harmonic number
    @param plot_sum_range: range to sum over
    @param exponent: constant a
    @param dimensionality: constant A
    @return: None
    """
    plt.figure(figsize=(8, 8))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2, font="Helvetica")
    residuals = theory_functions.plotting_statistics(dimensionality, exponent, n_points, half_n_harmonic,
                                                     plot_sum_range, simulation_means)[0]
    plt.scatter(plot_sum_range, residuals, color=_COLOUR_PALETTE["RESIDUALS"], label="Residuals")
    plt.hlines(0, plot_sum_range[0], plot_sum_range[-1], color=_COLOUR_PALETTE["THEORY"], label="Theory")
    plt.ylim(-5, 5)
    plt.legend()
    plt.xlabel("s")
    plt.ylabel("Residuals")
    plt.tight_layout()


def create_2d_plots(arguments: argparse.Namespace) -> None:
    """
    Plot amino acid distances and residuals
    @param arguments: command line arguments
    @return: None
    """
    dataframe = pd.read_csv(f"../data/simulations/2d/simulation_stats.csv")
    data_tuple = get_2d_plotting_data(dataframe)
    distances = data_tuple[1].astype(np.float)
    chain_length = np.float64(distances[-1] + 1)
    half_n_harmonic = np.float64(theory_functions.harmonic_number(chain_length // 2))
    dimensionality_range = np.arange(arguments.start_dimensionality,
                                     arguments.end_dimensionality,
                                     arguments.step_dimensionality)
    exponent_range = np.arange(arguments.start_exponent,
                               arguments.end_exponent,
                               arguments.step_exponent).astype(np.float)

    plot_2d_distances(distances=distances,
                      simulation_means=data_tuple[0],
                      lower_cl=data_tuple[2],
                      upper_cl=data_tuple[3],
                      chain_length=chain_length,
                      half_n_harmonic=half_n_harmonic,
                      exponent_range=exponent_range,
                      dimensionality_range=dimensionality_range,
                      start_point=np.float64(arguments.start_point),
                      end_point=np.float64(arguments.end_point))
    # plot_2d_residuals(simulation_means=data_tuple[0], chain_length=chain_length, half_n_harmonic=half_n_harmonic,
    # plot_sum_range=plotting_sum_range, exponent=exponent_range, dimensionality=dimensionality_range)
    plt.show()


def grid_plot_2d_distances(dimensionality_range: np.ndarray, exponent_range: np.ndarray, chain_length: int,
                           half_n_harmonic: float, plotting_sum_range: np.ndarray,
                           distances: np.ndarray, means: np.ndarray,
                           lower_cl: np.ndarray, upper_cl: np.ndarray, start_point: int, end_point: int) -> None:
    """
    Plot amino acid distances and residuals grids for different values of exponent and dimensionality
    @param dimensionality_range: range for dimensionality constant A
    @param exponent_range: range for exponent constant a
    @param chain_length: number of datapoints
    @param half_n_harmonic: the "N/2"-th harmonic number
    @param plotting_sum_range: range to sum over
    @param distances: amino acid distances
    @param means: mean frequencies of amino acid distances
    @param lower_cl: lower bound of confidence level
    @param upper_cl: upper bound for confidence level
    @param start_point: point from which to start plotting
    @param end_point: point at which to stop plotting
    @return: None
    """
    fig = plt.figure(figsize=(16, 10))
    ax = fig.subplots(len(exponent_range), len(dimensionality_range), sharex=True, sharey=True)
    sns.set(context="notebook", palette="colorblind", style="ticks", font="Helvetica")
    for row in range(len(exponent_range)):
        for col in range(len(dimensionality_range)):
            theory = theory_functions.amino_acid_distance_distribution(distances, chain_length, half_n_harmonic,
                                                                       exponent_range[row],
                                                                       dimensionality_range[col])

            ax[row][col].scatter(distances, means, s=10, label="2D SIM", c=_COLOUR_PALETTE["2D_SIM_SCATTER"])
            ax[row][col].plot(plotting_sum_range, theory, linestyle="--", label="Theory", c=_COLOUR_PALETTE["THEORY"],
                              lw=1.5)
            ax[row][col].fill_between(distances, upper_cl, lower_cl, color=_COLOUR_PALETTE["CL"],
                                      alpha=0.4, label="95% C.L.", zorder=-10)
            ax[row][col].set_yscale("log")
            ax[row][col].set_xscale("log")
            ax[row][col].set_ylim(0.1, 100)
            ax[row][col].set_xlim(start_point, end_point)
            ax[row][-1].set_ylabel(f"a = {exponent_range[row]}", fontsize=13, rotation=0, labelpad=21)
            ax[row][-1].yaxis.set_label_position("right")
            ax[0][col].set_title(f"A = {dimensionality_range[col]:2f}")
    fig.text(0.5, 0.025, "s", ha="center", fontsize=15.5)
    fig.text(0.005, 0.5, "P(s)", va="center", rotation="vertical", fontsize=15.5)
    plt.subplots_adjust(left=0.06, bottom=0.08, top=0.95, wspace=0.05, right=0.909)
    plt.legend(bbox_to_anchor=[1.36, 4.65], fontsize=9)
    plt.savefig("../plots/supplementary_information/simulation_2d_grid.pdf")


def grid_plot_2d_residuals(dimensionality_range, exponent_range, n_points, half_n_harmonic_number, plotting_sumrange,
                           means, start_point, end_point) -> None:
    """
    Create a grid of residual plots for different values of the exponent and dimensionality constants
    @param dimensionality_range:
    @param exponent_range:
    @param n_points:
    @param half_n_harmonic_number:
    @param plotting_sumrange:
    @param means:
    @param start_point:
    @param end_point:
    @return:
    """
    fig = plt.figure(figsize=(16, 10))
    ax = fig.subplots(len(exponent_range), len(dimensionality_range), sharex=True, sharey=True)
    for row in range(len(exponent_range)):
        for col in range(len(dimensionality_range)):
            residuals = theory_functions.plotting_statistics(dimensionality_range[col], exponent_range[row], n_points,
                                                             half_n_harmonic_number, plotting_sumrange, means)[0]
            ax[row][col].scatter(plotting_sumrange, residuals, marker=".", c=_COLOUR_PALETTE["RESIDUALS"],
                                 label="Residuals", zorder=10)
            ax[row][col].hlines(0, plotting_sumrange[0], plotting_sumrange[-1], color=_COLOUR_PALETTE["THEORY"],
                                label="Theory")
            ax[row][-1].set_ylabel(f"a = {exponent_range[row]}", fontsize=13, rotation=0, labelpad=21)
            ax[row][-1].yaxis.set_label_position("right")
            ax[0][col].set_title(f"A = {dimensionality_range[col]:2f}")
            ax[row][col].legend(fontsize=9)
            ax[row][col].set_xlim(start_point, end_point)

    fig.text(0.5, 0.025, "s", ha="center", fontsize=15.5)
    fig.text(0.005, 0.5, "Residuals", va="center", rotation="vertical", fontsize=15.5)
    plt.subplots_adjust(left=0.06, bottom=0.08, top=0.95, wspace=0.05, right=0.909)
    plt.savefig("../plots/supplementary_information/simulation_2d_residuals_grid.pdf")


def create_2d_grid_plots(arguments: argparse.Namespace) -> None:
    """
    Bring together everything needed to plot grid plots for 2D simulation amino acid distance distributions
    @param arguments: command line arguments
    @return: None
    """
    get_2d_stats()
    dataframe = pd.read_csv("../data/simulations/2d/simulation_stats.csv")
    data_tuple = get_2d_plotting_data(dataframe)
    distances = data_tuple[1]
    chain_length = int(distances[-1] + 1)
    plotting_sum_range = range(int(distances[0]), chain_length)
    half_n_harmonic = theory_functions.harmonic_number(chain_length // 2)
    dimensionality_range = np.arange(arguments.start_dimensionality,
                                     arguments.end_dimensionality,
                                     arguments.step_dimensionality)
    exponent_range = np.arange(arguments.start_exponent,
                               arguments.end_exponent,
                               arguments.step_exponent)
    grid_plot_2d_distances(dimensionality_range=dimensionality_range,
                           exponent_range=exponent_range,
                           chain_length=chain_length,
                           half_n_harmonic=half_n_harmonic,
                           plotting_sum_range=plotting_sum_range,
                           distances=distances,
                           means=data_tuple[0],
                           lower_cl=data_tuple[2],
                           upper_cl=data_tuple[3],
                           start_point=arguments.start_point,
                           end_point=arguments.end_point)
    grid_plot_2d_residuals(dimensionality_range=dimensionality_range,
                           exponent_range=exponent_range,
                           n_points=chain_length,
                           half_n_harmonic_number=half_n_harmonic,
                           plotting_sumrange=plotting_sum_range,
                           means=data_tuple[0],
                           start_point=arguments.start_point,
                           end_point=arguments.end_point)
    plt.show()


def get_3d_histogram(arguments: argparse.Namespace) -> np.ndarray:
    """
    Return histogram of 3D data for given length range
    @param arguments: command line arguments
    @return: numpy array of histogram
    """
    return np.load(f"../data/simulations/3d/histogram_{arguments.length_range}_not_normed.npy", allow_pickle=True)


def get_3d_data_for_plotting(histogram: np.ndarray, arguments: argparse.Namespace) -> tuple:
    """
    Get data for plotting
    @param histogram: numpy array of data
    @param arguments: command line arguments
    @return: tuple of number of datapoints, measure of central tendency and confidence level bounds
    """
    distances = np.linspace(start=1, stop=300, num=300)[:-1]
    lower_bound, upper_bound = theory_functions.get_confidence_interval(histogram, arguments.quantile)
    measure = theory_functions.get_measure_of_central_tendency(histogram, arguments.measure)
    normalised_measure = measure / np.sum(measure)
    normalised_lower_bound = lower_bound / np.sum(measure)
    normalised_upper_bound = upper_bound / np.sum(measure)
    chain_length = int(arguments.length_range)
    return chain_length, distances, normalised_measure, normalised_lower_bound, normalised_upper_bound


def plot_3d_distances(distances: np.ndarray, measure: np.ndarray,
                      lower_confidence: np.ndarray, upper_confidence: np.ndarray,
                      chain_length: int, half_n_harmonic: float,
                      exponent_range: np.ndarray, dimensionality_range: np.ndarray,
                      start_point: int, end_point: int, arguments: argparse.Namespace) -> None:
    """
    Plot amino acid distance distribution from 3D simulations
    @param distances: amino acid distances
    @param measure: mean or median
    @param lower_confidence: lower confidence bound
    @param upper_confidence: upper confidence bound
    @param chain_length: length of chain
    @param half_n_harmonic: the "N/2"-th harmonic number
    @param exponent_range: range for exponent constant
    @param dimensionality_range: range for dimensionality constant
    @param start_point: point from which to start plotting
    @param end_point: point at which to stop plotting
    @param arguments: command line arguments
    @return: None
    """
    random.seed(1234)
    exponent_index = random.randint(0, len(exponent_range) - 1)
    dimensionality_index = random.randint(0, len(dimensionality_range) - 1)
    print(exponent_index, dimensionality_index)
    # parameters, covariance = scipy.optimize.curve_fit(f=theory_functions.amino_acid_distance_distribution,
    #                                                   xdata=distances, ydata=measure,
    #                                                   p0=[chain_length, half_n_harmonic,
    #                                                       exponent_range[exponent_index],
    #                                                       dimensionality_range[dimensionality_index]],
    #                                                   sigma=upper_confidence-lower_confidence, absolute_sigma=True)

    parameters = [chain_length, half_n_harmonic,
                  exponent_range[exponent_index], dimensionality_range[dimensionality_index]]
    theory_x = np.linspace(start_point, int(chain_length//2), int(chain_length//2))[:-1]
    least_squares_result = scipy.optimize.least_squares(theory_functions.vector_of_residuals, parameters,
                                                        args=(theory_x, measure[:len(theory_x)]),
                                                        bounds=([0, 0, exponent_range[0], dimensionality_range[0]],
                                                                [chain_length, half_n_harmonic,
                                                                 exponent_range[-1], dimensionality_range[-1]]))
    print(least_squares_result)
    optimised_parameters = least_squares_result.x
    jacobian = least_squares_result.jac
    covariance_matrix = theory_functions.covariance_matrix(jacobian)
    print(covariance_matrix)
    sigma = theory_functions.error_on_least_squares(covariance_matrix)
    print(sigma)
    fig = plt.figure()
    fig.set_size_inches((8, 8))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=1.8, font="Helvetica")

    plt.scatter(distances, measure, label=f"3D Simulation {arguments.length_range}",
                color=_COLOUR_PALETTE["3D_SIM_SCATTER"], marker="o")
    plt.fill_between(distances, upper_confidence, lower_confidence,
                     color=_COLOUR_PALETTE["CL"], alpha=0.4,
                     label=r"3D Simulation %i$\sigma$ C.L." % arguments.quantile, zorder=-99)
    plt.plot(theory_x, theory_functions.amino_acid_distance_distribution(theory_x, *optimised_parameters), label="Theory", color=_COLOUR_PALETTE["3D_SIM_SCATTER"], lw=1.5, ls="--")
    print("----------------Parameters----------------")
    print(f"Chain length: {optimised_parameters[0]}")
    print(f"N/2 Harmonic: {optimised_parameters[1]}")
    print(f"Exponent: {optimised_parameters[2]}")
    print(f"Dimensionality: {optimised_parameters[3]}")
    print("----------------Covariance----------------")
    print(covariance_matrix)
    print("-------------------Sigma------------------")
    print(f"Chain length: {sigma[0]}")
    print(f"N/2 Harmonic: {sigma[1]}")
    print(f"Exponent: {sigma[2]}")
    print(f"Dimensionality: {sigma[3]}")
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("s")
    plt.ylabel("P(s)")
    plt.ylim(0.00001, 10)
    # plt.xlim(start_point, end_point)
    plt.legend(loc="upper right")
    sns.despine()
    plt.tight_layout()


def create_3d_simulation_plot(arguments):
    histogram = get_3d_histogram(arguments)
    data_tuple = get_3d_data_for_plotting(histogram, arguments)
    distances = data_tuple[1]
    chain_length = data_tuple[0]
    half_n_harmonic = np.float64(theory_functions.harmonic_number(chain_length // 2))
    dimensionality_range = np.arange(arguments.start_dimensionality,
                                     arguments.end_dimensionality,
                                     arguments.step_dimensionality)
    exponent_range = np.arange(arguments.start_exponent,
                               arguments.end_exponent,
                               arguments.step_exponent).astype(np.float)

    plot_3d_distances(distances=distances,
                      measure=data_tuple[2],
                      lower_confidence=data_tuple[3],
                      upper_confidence=data_tuple[4],
                      chain_length=chain_length,
                      half_n_harmonic=half_n_harmonic,
                      exponent_range=exponent_range,
                      dimensionality_range=dimensionality_range,
                      start_point=np.float64(arguments.start_point),
                      end_point=np.float64(arguments.end_point),
                      arguments=arguments)
    plt.show()
