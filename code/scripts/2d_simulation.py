import glob
import numpy as np
import pandas as pd
import sequence_distance_distribution.code.scripts.theory_functions as theory_functions
import matplotlib.pyplot as plt
import seaborn as sns
from sequence_distance_distribution.code.scripts.plot_functions import _COLOUR_PALETTE


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


def plot_2d_simulation_distances(distances: np.ndarray, simulation_means: np.ndarray, lower_cl: np.ndarray,
                                 upper_cl: np.ndarray, plot_sum_range: np.ndarray, n_points: int,
                                 half_n_harmonic: float, exponent: int, dimensionality: float) -> None:
    """
    Plot 2D simulation distance distribution and compare to theory
    @param distances: amino acid distances in chain length range
    @param simulation_means: mean frequencies of amino acid distances
    @param lower_cl: lower bound for confidence interval
    @param upper_cl: upper bound for confidence interval
    @param plot_sum_range: range to sum over
    @param n_points: number of points
    @param half_n_harmonic: the "N/2"-th harmonic number
    @param exponent: constant a
    @param dimensionality: constant A
    @return: None
    """
    theory = [theory_functions.amino_acid_distance_distribution(s, n_points, half_n_harmonic, exponent, dimensionality)
              for s in plot_sum_range]
    plt.figure(figsize=(8, 8))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2.4, font="Helvetica")
    plt.scatter(distances, simulation_means, label="2D Simulation", color=_COLOUR_PALETTE["2D_SIM_SCATTER"])
    plt.plot(plot_sum_range, theory, linestyle="--", label="Theory", color=_COLOUR_PALETTE["THEORY"], lw=1.5)
    plt.fill_between(distances, upper_cl, lower_cl, color=_COLOUR_PALETTE["CL"], zorder=-10, label="95% C.L.",
                     alpha=0.4)
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("s")
    plt.ylabel("P(s)")
    plt.xlim(2, 250)
    plt.legend(loc="upper right")
    sns.despine()
    plt.tight_layout()


def plot_2d_simulation_residuals(simulation_means: np.ndarray, n_points: int, half_n_harmonic: float,
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


def create_2d_plots(exponent: int, dimensionality: float):
    """
    @param exponent: constant a
    @param dimensionality: constant A
    @return:
    """
    dataframe = pd.DataFrame(f"../data/simulations/2d/simulation_stats.csv")
    data_tuple = get_2d_plotting_data(dataframe)
    distances = data_tuple[1]
    n_points = int(distances[-1] + 1)
    plotting_sum_range = range(int(distances[0]), n_points)
    half_n_harmonic = theory_functions.harmonic_number(n_points // 2)
    plot_2d_simulation_distances(distances=distances,
                                 simulation_means=data_tuple[0],
                                 lower_cl=data_tuple[2],
                                 upper_cl=data_tuple[3],
                                 plot_sum_range=plotting_sum_range,
                                 n_points=n_points,
                                 half_n_harmonic=half_n_harmonic,
                                 exponent=exponent,
                                 dimensionality=dimensionality)
    plot_2d_simulation_residuals(simulation_means=data_tuple[0],
                                 n_points=n_points,
                                 half_n_harmonic=half_n_harmonic,
                                 plot_sum_range=plotting_sum_range,
                                 exponent=exponent,
                                 dimensionality=dimensionality)
    plt.show()
