"""
Functions used for plotting amino acid distance distributions.
"""
import matplotlib.pyplot as plt
import numpy as np
import theory_functions
import seaborn as sns

_COLOUR_PALETTE = {"PDB_SCATTER": "#006374",
                   "SIM_SCATTER": "#fbafe4",
                   "THEORY": "#006374",
                   "RESIDUALS": "#fbafe4"}


def normalise_data(data_array: np.ndarray) -> np.ndarray:
    """
    Take a data array and width of bins and normalise the data array.
    @return: normalised_array
    """
    data_sum = np.sum(data_array)
    return data_array / data_sum


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


def grid_plot_distances(dimensionality_range: np.ndarray, exponent_range: np.ndarray,
                        n_points: int, half_n_harmonic: float, plotting_sumrange: np.ndarray,
                        normalised_means: np.ndarray, distances: np.ndarray,
                        normalised_sim_means: np.ndarray, sim_distances: np.ndarray,
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
    sns.set(context="notebook", palette="colorblind", style='ticks', font_scale=2.4, font='Helvetica')
    ax = fig.subplots(len(exponent_range), len(dimensionality_range), sharex=True, sharey=True)
    for row in range(len(exponent_range)):
        for col in range(len(dimensionality_range)):
            theory = [
                theory_functions.amino_acid_distance_distribution(s, n_points, half_n_harmonic, exponent_range[row],
                                                                  dimensionality_range[col]) for s in plotting_sumrange]
            ax[row][col].scatter(distances, normalised_means, s=10, c=_COLOUR_PALETTE["PDB_SCATTER"], label=plot_label)
            ax[row][col].fill_between(distances, upper_cl, lower_cl, c=_COLOUR_PALETTE["PDB_SCATTER"], alpha=0.4,
                                      label="95% CL", zorder=-100)
            ax[row][col].scatter(sim_distances, normalised_sim_means, s=10, m="^", c=_COLOUR_PALETTE["SIM_SCATTER"],
                                 label=f"SIM {length_range}s", zorder=-50)
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
    sns.set(context="notebook", palette="colorblind", style='ticks', font_scale=2.4, font='Helvetica')
    ax = fig.subplots(len(exponent_range), len(dimensionality_range), sharex=True, sharey=True)
    for row in range(len(exponent_range)):
        for col in range(len(dimensionality_range)):
            residuals = theory_functions.plotting_statistics(dimensionality_range[col], exponent_range[row],
                                                             n_points, half_n_harmonic_number, plotting_sumrange,
                                                             normalised_means)[0]

            ax[row][col].scatter(plotting_sumrange, residuals, m=".", c=_COLOUR_PALETTE["RESIDUALS"], label="Residuals",
                                 zorder=10)
            ax[row][col].hlines(0, plotting_sumrange[0], plotting_sumrange[-1], c=_COLOUR_PALETTE["THEORY"],
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
