"""Functions for plotting comparison plot (Fig.4.)"""
import argparse
import numpy as np
import sequence_distance_distribution.code.scripts.theory_functions as theory_functions
from sequence_distance_distribution.code.scripts.plot_functions import _COLOUR_PALETTE
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
                                                                     arguments.e_rcsb,
                                                                     arguments.d_rcsb) for s in rcsb_sum_range]
    alphafold_theory = [theory_functions.amino_acid_distance_distribution(s,
                                                                          alphafold_n_points,
                                                                          alphafold_half_n_harmonic,
                                                                          arguments.e_c,
                                                                          arguments.d_c) for s in alphafold_sum_range]
    plt.figure(figsize=(8, 8))
    sns.set(context="notebook", palette="colorblind", style='ticks', font_scale=1.8, font='Helvetica')
    plt.scatter(sim_tuple[0], sim_tuple[1], label="3D Simulation", color=_COLOUR_PALETTE["SIM_SCATTER"], marker="^")
    plt.fill_between(sim_tuple[0], sim_tuple[3], sim_tuple[2], color=_COLOUR_PALETTE["SIM_SCATTER"], alpha=0.4,
                     label="3D Simulation 95% C.L.", zorder=-100)
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
