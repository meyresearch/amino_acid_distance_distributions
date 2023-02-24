"""Functions for plotting comparison plot (Fig.4.b.)"""
import argparse
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import seaborn as sns
import plot_functions
import theory_functions
from colour_palette import _COLOUR_PALETTE


def create_comparison_plot(arguments: argparse.Namespace, rcsb_histogram: np.ndarray,
                           alpha_histogram: np.ndarray) -> None:
    """
    Create amino acid distance distribution plot (Fig.4.b) and compare to theory and power law
    @param arguments: command line arguments from user
    @param rcsb_histogram: numpy array containing RCSB PDB data
    @param alpha_histogram: numpy array containing AlphaFold 2 PDB data
    @return: None
    """
    rcsb_plotting_tuple = plot_functions.get_data_for_plotting(rcsb_histogram, arguments, arguments.length_range)
    rcsb_chain_length = rcsb_plotting_tuple[0]
    rcsb_distances = rcsb_plotting_tuple[1]
    rcsb_measure = rcsb_plotting_tuple[2]
    rcsb_lower_bound, rcsb_upper_bound = rcsb_plotting_tuple[3], rcsb_plotting_tuple[4]
    half_n_harmonic = theory_functions.harmonic_number(n_numbers=(rcsb_chain_length // 2))
    alpha_plotting_tuple = plot_functions.get_data_for_plotting(alpha_histogram, arguments, arguments.length_range)
    alpha_distances = alpha_plotting_tuple[1]
    alpha_measure = alpha_plotting_tuple[2]

    dimensionality_range = np.arange(arguments.rcsb_startd,
                                     arguments.rcsb_endd,
                                     arguments.step_dimensionality)
    exponent_range = np.arange(arguments.rcsb_starte,
                               arguments.rcsb_ende,
                               arguments.step_exponent)
    start = arguments.start_point
    end = arguments.end_point
    weights = rcsb_plotting_tuple[4] - rcsb_plotting_tuple[3]
    rcsb_parameters, rcsb_cov = scipy.optimize.curve_fit(f=theory_functions.amino_acid_distance_distribution,
                                                         xdata=rcsb_distances[start:end], ydata=rcsb_measure[start:end],
                                                         sigma=weights[start:end],
                                                         bounds=([rcsb_chain_length - 1, half_n_harmonic - 0.00001,
                                                                  exponent_range[0], dimensionality_range[0]],
                                                                 [rcsb_chain_length, half_n_harmonic,
                                                                  exponent_range[-1], dimensionality_range[-1]]))
                                                                  
    rcsb_theory = theory_functions.amino_acid_distance_distribution(rcsb_distances[start:end], *rcsb_parameters)
    rcsb_sigma = np.sqrt(np.diag(rcsb_cov))

    ks_index = int(rcsb_chain_length) + 15
    ks_statistics = scipy.stats.ks_2samp(rcsb_measure[:ks_index], alpha_measure[:ks_index])
    alpha_confidence = 0.01
    is_accepted = theory_functions.accept_null_hypothesis(ks_statistics.pvalue, alpha_confidence)
    fig = plt.figure()
    fig.set_size_inches((6, 6))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=1.88)

    plt.scatter(rcsb_distances, rcsb_measure, label=f"RCSB {arguments.length_range}",
                color=_COLOUR_PALETTE["PDB_SCATTER"], marker="o")
    plt.fill_between(rcsb_distances, rcsb_upper_bound, rcsb_lower_bound,
                     color=_COLOUR_PALETTE["PDB_SCATTER"], alpha=0.25, label="RCSB 95% C.L.", zorder=-99)
    # plt.scatter(alpha_distances, alpha_measure, label=f"AlphaFold 2 {arguments.length_range}",
    #             color=_COLOUR_PALETTE["ALPHA_SCATTER"], marker="s")

    plt.plot(rcsb_distances[start:end], rcsb_theory, label="Theory RCSB", color="k", lw=1.5)

    print("\n")
    print("-----------------------RCSB THEORY-----------------------")
    print(f"N: {rcsb_parameters[0]:f} +/- {rcsb_sigma[0]:f} ")
    print(f"H-N/2: {rcsb_parameters[1]:f} +/- {rcsb_sigma[1]:f}")
    print(f"a: {rcsb_parameters[2]:f} +/- {rcsb_sigma[2]:f}")
    print(f"A: {rcsb_parameters[3]:f} +/- {rcsb_sigma[3]:f}")
    print("\n")
    print("-----------------------KS STATISTICS----------------------")
    print(f"KS statistic: {ks_statistics.statistic}")
    print(f"p value: {ks_statistics.pvalue}")
    print(f"Null hypothesis accepted at alpha={alpha_confidence}: {is_accepted}")
    print("\n")
    plt.yscale("log")
    plt.xscale("log")
    plt.xlim(4, arguments.end_point)
    plt.ylim(0.0005, 0.1)
    plt.xlabel("s")
    plt.ylabel("P(s)")
    plt.legend(loc="upper right", fontsize=16)
    sns.despine()
    plt.tight_layout()
    plt.show()
    plt.savefig(f"../plots/supplementary_material/comparison_{rcsb_chain_length}.pdf")


def create_low_condidence_comparison_plot(arguments: argparse.Namespace, rcsb_histogram: np.ndarray,
                                          alpha_histogram: np.ndarray):

    rcsb_plotting_tuple = plot_functions.get_data_for_plotting(rcsb_histogram, arguments, arguments.length_range)
    rcsb_chain_length = rcsb_plotting_tuple[0]
    rcsb_distances = rcsb_plotting_tuple[1]
    rcsb_measure = rcsb_plotting_tuple[2]
    rcsb_lower_bound, rcsb_upper_bound = rcsb_plotting_tuple[3], rcsb_plotting_tuple[4]

    alpha_plotting_tuple = plot_functions.get_data_for_plotting(alpha_histogram, arguments, arguments.length_range)
    alpha_distances = alpha_plotting_tuple[1]
    alpha_measure = alpha_plotting_tuple[2]

    ks_index = int(rcsb_chain_length) + 15
    ks_statistics = scipy.stats.ks_2samp(rcsb_measure[:ks_index], alpha_measure[:ks_index])
    alpha_confidence = 0.01
    is_accepted = theory_functions.accept_null_hypothesis(ks_statistics.pvalue, alpha_confidence)

    fig = plt.figure()
    fig.set_size_inches((6, 6))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=1.88)

    plt.scatter(rcsb_distances, rcsb_measure, label=f"RCSB {arguments.length_range}",
                color=_COLOUR_PALETTE["PDB_SCATTER"], marker="o")
    plt.fill_between(rcsb_distances, rcsb_upper_bound, rcsb_lower_bound,
                     color=_COLOUR_PALETTE["PDB_SCATTER"], alpha=0.25, label="RCSB 95% C.L.", zorder=-99)
    plt.scatter(alpha_distances, alpha_measure, label=f"AlphaFold 2 {arguments.length_range} low-conf.",
                color=_COLOUR_PALETTE["ALPHA_SCATTER"], marker="s")
    print("\n")
    print("-----------------------KS STATISTICS----------------------")
    print(f"KS statistic: {ks_statistics.statistic}")
    print(f"p value: {ks_statistics.pvalue}")
    print(f"Null hypothesis accepted at alpha={alpha_confidence}: {is_accepted}")
    print("\n")
    plt.yscale("log")
    plt.xscale("log")
    plt.xlim(4, arguments.end_point)
    plt.ylim(0.000005, 0.1)
    plt.xlabel("s")
    plt.ylabel("P(s)")
    plt.legend(loc="lower left", fontsize=16)
    sns.despine()
    plt.tight_layout()
    plt.show()
    plt.savefig("../plots/supplementary_material/low_confidence_comparison.jpeg", dpi=2600)

