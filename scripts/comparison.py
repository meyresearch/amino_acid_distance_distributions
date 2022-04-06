"""Functions for plotting comparison plot (Fig.4.)"""
import argparse
import numpy as np
import scipy.optimize
import distances
import theory_functions
from colour_palette import _COLOUR_PALETTE
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


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
    filtered_rcsb_histogram = distances.filter_histogram(rcsb_histogram, arguments.start_point, arguments.end_point)
    rcsb_plotting_tuple = distances.get_data_for_plotting(filtered_rcsb_histogram, arguments)
    rcsb_distances = rcsb_plotting_tuple[1]
    rcsb_n_points = rcsb_plotting_tuple[0]
    rcsb_half_n_harmonic = theory_functions.harmonic_number(n_numbers=(rcsb_n_points // 2))
    rcsb_plotting_sum_range = rcsb_distances
    filtered_alpha_histogram = distances.filter_histogram(alpha_histogram, arguments.start_point, arguments.end_point)
    alpha_plotting_tuple = distances.get_data_for_plotting(filtered_alpha_histogram, arguments)
    alpha_distances = alpha_plotting_tuple[1]
    alpha_n_points = alpha_plotting_tuple[0]
    alpha_half_n_harmonic = theory_functions.harmonic_number(n_numbers=(alpha_n_points // 2))
    alpha_plotting_sum_range = alpha_distances
    filtered_sim_histogram = distances.filter_histogram(sim_histogram, arguments.start_point, arguments.end_point)
    sim_plotting_tuple = distances.get_data_for_plotting(filtered_sim_histogram, arguments)
    sim_distances = sim_plotting_tuple[1]
    # ----------------------------------------------------------------------------
    # need to change command line arguments to have rcsb and alpha startd, e, etc.
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
    # ----------------------------------------------------------------------------
    fig = plt.figure()
    fig.set_size_inches((8, 8))
    sns.set(context="notebook", palette="colorblind", style='ticks', font_scale=1.8, font="Helvetica")
    sim_err = sim_plotting_tuple[4] - sim_plotting_tuple[3]
    plt.scatter(sim_distances, sim_plotting_tuple[2], label="Simulation",
                 color=_COLOUR_PALETTE["3D_SIM_SCATTER"], marker="^")
    plt.fill_between(sim_distances, sim_plotting_tuple[4], sim_plotting_tuple[3], color=_COLOUR_PALETTE["3D_SIM_SCATTER"],
                      alpha=0.4, label="Simulation 95% C.L.", zorder=-100)
    rcsb_err = rcsb_plotting_tuple[4] - rcsb_plotting_tuple[3]

    plt.scatter(rcsb_distances, rcsb_plotting_tuple[2], label=f"RCSB {arguments.length_range}",
                color=_COLOUR_PALETTE["PDB_SCATTER"], marker="o")

    # (_, caps, _) = plt.errorbar(rcsb_distances, rcsb_plotting_tuple[2], yerr=rcsb_err, label=f"RCSB {arguments.length_range}",
    #             color=_COLOUR_PALETTE["PDB_SCATTER"], fmt="o", capsize=5)
    # for cap in caps:
    #     cap.set_color(_COLOUR_PALETTE["PDB_SCATTER"])
    #     cap.set_markeredgewidth(2)
    # plt.fill_between(rcsb_distances, rcsb_plotting_tuple[4], rcsb_plotting_tuple[3],
    #                  color=_COLOUR_PALETTE["PDB_SCATTER"], alpha=0.2, label="RCSB 95% C.L.", zorder=-99)
    alpha_err = alpha_plotting_tuple[4] - alpha_plotting_tuple[3]
    plt.scatter(alpha_distances, alpha_plotting_tuple[2], label=f"AlphaFold 2 {arguments.length_range}",
                color=_COLOUR_PALETTE["ALPHA_SCATTER"], marker="s")
    # plt.fill_between(alpha_distances, alpha_plotting_tuple[4], alpha_plotting_tuple[3],
    #                  color=_COLOUR_PALETTE["ALPHA_SCATTER"], alpha=0.2, label="AlphaFold 2 95% C.L.", zorder=-98)

    rcsb_parameters, rcsb_cov = scipy.optimize.curve_fit(
        f=theory_functions.amino_acid_distance_distribution,
        xdata=rcsb_plotting_sum_range, ydata=rcsb_plotting_tuple[2],
        bounds=([rcsb_n_points - 1, rcsb_half_n_harmonic - 0.00001,
                 rcsb_exponent_range[0], rcsb_dimensionality_range[0]],
                [rcsb_n_points, rcsb_half_n_harmonic,
                 rcsb_exponent_range[-1], rcsb_dimensionality_range[-1]]))
    rcsb_power, rcsb_power_cov = scipy.optimize.curve_fit(f=theory_functions.power_law,
                                                          xdata=rcsb_plotting_sum_range, ydata=rcsb_plotting_tuple[2])
    rcsb_theory = theory_functions.amino_acid_distance_distribution(rcsb_plotting_sum_range, *rcsb_parameters)
    rcsb_powerlaw = theory_functions.power_law(rcsb_plotting_sum_range, *rcsb_power)
    alpha_parameters, alpha_cov = scipy.optimize.curve_fit(
        f=theory_functions.amino_acid_distance_distribution,
        xdata=alpha_plotting_sum_range, ydata=alpha_plotting_tuple[2],
        bounds=([alpha_n_points - 1, alpha_half_n_harmonic - 0.00001,
                 alpha_exponent_range[0], alpha_dimensionality_range[0]],
                [alpha_n_points, alpha_half_n_harmonic,
                 alpha_exponent_range[-1], alpha_dimensionality_range[-1]]))
    # alpha_theory = theory_functions.amino_acid_distance_distribution(alpha_plotting_sum_range, *alpha_parameters)
    # alpha_power, alpha_power_cov = scipy.optimize.curve_fit(f=theory_functions.power_law,
    #                                                         xdata=alpha_plotting_sum_range,
    #                                                         ydata=alpha_plotting_tuple[2])

    # alpha_powerlaw = theory_functions.power_law(alpha_plotting_sum_range, *alpha_power)
    plt.plot(rcsb_plotting_sum_range, rcsb_theory, label="Theory RCSB", color="k", lw=1.5)
    plt.plot(rcsb_plotting_sum_range, rcsb_powerlaw, label="Power law RCSB", color="k", lw=1.5,
             ls="--")

    # plt.plot(alpha_plotting_sum_range, alpha_theory, label="Theory AlphaFold 2", color=_COLOUR_PALETTE["ALPHA_SCATTER"],
    #          lw=1.5)
    # plt.plot(alpha_plotting_sum_range, alpha_powerlaw, label="Power law AlphaFold 2",
    #          color=_COLOUR_PALETTE["POWER"], lw=1.5, ls="-.")

    print("-----------------------RCSB THEORY-----------------------")
    print(f"N: {rcsb_parameters[0]}")
    print(f"H-N/2: {rcsb_parameters[1]}")
    print(f"a: {rcsb_parameters[2]}")
    print(f"A: {rcsb_parameters[3]}")
    print("----------------------RCSB POWER LAW----------------------")
    print(f"gamma: {rcsb_power[0]}")
    print(f"constant: {rcsb_power[1]}")
    # print("--------------------AlphaFold 2 THEORY--------------------")
    # print(f"N: {alpha_parameters[0]}")
    # print(f"H-N/2: {alpha_parameters[1]}")
    # print(f"a: {alpha_parameters[2]}")
    # print(f"A: {alpha_parameters[3]}")
    # print("-------------------AlphaFold 2 POWER LAW------------------")
    # print(f"gamma: {alpha_power[0]}")
    # print(f"constant: {alpha_power[1]}")

    plt.yscale("log")
    plt.xscale("log")
    plt.xlim(10, int(arguments.length_range)//2)
    plt.ylim(0.001, 0.2)
    plt.xlabel("s")
    plt.ylabel("P(s)")
    plt.legend(loc="upper right")
    sns.despine()
    plt.tight_layout()
    plt.show()
