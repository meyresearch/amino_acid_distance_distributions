"""Functions for plotting 2D and 3D simulation amino acid distances."""
import argparse
import glob
import numpy as np
import pandas as pd
import scipy.optimize
import matplotlib.pyplot as plt
import seaborn as sns
import theory_functions
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
    mean_frequencies = simulation_dataframe["mean"].to_numpy()
    normalised_means = mean_frequencies / np.sum(mean_frequencies)
    distances = simulation_dataframe["variable"].to_numpy()
    lower_bound = simulation_dataframe["lower_bound"].to_numpy()
    normalised_lower_bound = lower_bound / np.sum(mean_frequencies)
    upper_bound = simulation_dataframe["upper_bound"].to_numpy()
    normalised_upper_bound = upper_bound / np.sum(mean_frequencies)
    return normalised_means, distances, normalised_lower_bound, normalised_upper_bound


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
    weights = upper_cl - lower_cl
    parameters, covariance = scipy.optimize.curve_fit(f=theory_functions.amino_acid_distance_distribution,
                                                      xdata=distances[start_point:end_point],
                                                      ydata=simulation_means[start_point:end_point],
                                                      sigma=weights[start_point:end_point],
                                                      bounds=([chain_length - 1, half_n_harmonic - 0.0001,
                                                               exponent_range[0],
                                                               dimensionality_range[0]],
                                                              [chain_length, half_n_harmonic,
                                                               exponent_range[-1],
                                                               dimensionality_range[-1]]))

    theory = theory_functions.amino_acid_distance_distribution(distances[start_point:end_point], *parameters)
    sigma = np.sqrt(np.diag(covariance))
    plt.figure(figsize=(8, 8))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2.88, font="Helvetica")
    plt.scatter(distances, simulation_means, label="2D Simulation", color=_COLOUR_PALETTE["2D_SIM_SCATTER"])
    plt.plot(distances[start_point:end_point], theory, linestyle="--", label="Theory", color=_COLOUR_PALETTE["THEORY"],
             lw=1.5)
    plt.fill_between(distances, upper_cl, lower_cl, color=_COLOUR_PALETTE["CL"], zorder=-10, label="95% C.L.",
                     alpha=0.4)
    print("\n")
    print("----------------Parameters----------------")
    print(f"Chain length: {parameters[0]} +/- {sigma[0]}")
    print(f"N/2 Harmonic: {parameters[1]} +/- {sigma[1]}")
    print(f"Exponent: {parameters[2]} +/- {sigma[2]}")
    print(f"Dimensionality: {parameters[3]} +/- {sigma[3]}")
    print("\n")
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("s")
    plt.ylabel("P(s)")
    plt.xlim(start_point, end_point)
    plt.legend(loc="upper right", fontsize=24, frameon=False)
    sns.despine()
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
    chain_length = distances[-1] + 1
    half_n_harmonic = theory_functions.harmonic_number(chain_length // 2)
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
    distances = np.linspace(start=1, stop=350, num=350)[:-1]
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
    weights = upper_confidence - lower_confidence
    parameters, covariance = scipy.optimize.curve_fit(f=theory_functions.amino_acid_distance_distribution,
                                                      xdata=distances[start_point:end_point],
                                                      ydata=measure[start_point:end_point],
                                                      sigma=weights[start_point:end_point],
                                                      bounds=([chain_length - 1, half_n_harmonic - 0.00001,
                                                               exponent_range[0], dimensionality_range[1]],
                                                              [chain_length, half_n_harmonic,
                                                               exponent_range[-1], dimensionality_range[-1]]))

    theory = theory_functions.amino_acid_distance_distribution(distances[start_point:end_point], *parameters)
    sigma = np.sqrt(np.diag(covariance))
    fig = plt.figure()
    fig.set_size_inches((8, 8))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2.88, font="Helvetica")
    plt.scatter(distances, measure, label=f"3D Simulation {arguments.length_range}",
                color=_COLOUR_PALETTE["3D_SIM_SCATTER"], marker="o")
    plt.plot(distances[start_point:end_point], theory, label="Theory", color=_COLOUR_PALETTE["3D_SIM_SCATTER"],
             lw=1.5, ls="--")
    plt.fill_between(distances, upper_confidence, lower_confidence,
                     color=_COLOUR_PALETTE["CL"], alpha=0.4,
                     label="95% C.L.", zorder=-99)
    print("\n")
    print("----------------Parameters----------------")
    print(f"Chain length: {parameters[0]} +/- {sigma[0]}")
    print(f"N/2 Harmonic: {parameters[1]} +/- {sigma[1]}")
    print(f"Exponent: {parameters[2]} +/- {sigma[2]}")
    print(f"Dimensionality: {parameters[3]} +/- {sigma[3]}")
    print("\n")
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("s")
    plt.ylabel("P(s)")
    plt.ylim(0.001, 1)
    plt.xlim(start_point, end_point)
    plt.legend(loc="upper right", fontsize=24, frameon=False)
    sns.despine()
    plt.tight_layout()


def create_3d_simulation_plot(arguments):
    histogram = get_3d_histogram(arguments)
    data_tuple = get_3d_data_for_plotting(histogram, arguments)
    distances = data_tuple[1]
    chain_length = data_tuple[0]
    half_n_harmonic = theory_functions.harmonic_number(chain_length // 2)
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
                      start_point=arguments.start_point,
                      end_point=arguments.end_point,
                      arguments=arguments)
    plt.show()
