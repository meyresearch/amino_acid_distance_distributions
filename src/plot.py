"""
plot_link_lengths.py
"""
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from numpy import ndarray
from sklearn.metrics import r2_score


def normalise_data(data_array: np.ndarray) -> np.ndarray:
    """
    Take a data array and width of bins and normalise the data array.
    @return: normalised_array
    """
    data_sum: ndarray = np.sum(data_array)
    normalised_array = data_array / data_sum
    return normalised_array


def nth_harmonic(n: int) -> float:
    """
    @param n: number of natural numbers
    @return: harmonic: nth harmonic number
    """
    n_int = int(n)
    harmonic = 1.00
    for i in range(2, n_int + 1):
        harmonic += 1 / i
    return harmonic


def f_high(link_length_s: int, chain_length: int, half_n_harmonic_number: float) -> float:
    """
    Calculate f for k >> 1
    @param link_length_s: separation between each link, 2 <= s < N/2
    @param chain_length: number of C-alphas in a chain
    @param half_n_harmonic_number: the 'N/2'-th harmonic number
    @return: sequence distribution evaluated at a link length with k >> 1
    """
    harmonic_number = nth_harmonic(link_length_s)
    constant_multiplier = 2 / chain_length
    s_term = ((half_n_harmonic_number - harmonic_number + 1) / (half_n_harmonic_number - 1)) * link_length_s
    constant_term = half_n_harmonic_number / (half_n_harmonic_number - 1)
    return constant_multiplier * (s_term - constant_term)


def f_low(link_length_s: int, chain_length: int) -> float:
    """
    @param link_length_s: separation between each link, 2 <= s < N/2
    @param chain_length: number of C-alphas in a chain
    @return: sequence distribution evaluated at a given link_length with k << N/2
    """
    s_square_term = (-2 / chain_length ** 2) * link_length_s ** 2
    s_term = (2 / chain_length + 6 / chain_length ** 2) * link_length_s
    constant_term = (2 / chain_length) + (4 / chain_length ** 2)
    return s_square_term + s_term - constant_term


# noinspection PyPep8Naming,PyShadowingNames
def link_length_distribution(link_length: int, chain_length: int, half_n_harmonic: float,
                             constant_a: int, constant_A: float) -> float:
    """

    @param link_length: separation between each link 2 <= s < N/2
    @param chain_length: number of C-alphas
    @param half_n_harmonic: the 'N/2'-th harmonic number
    @param constant_a: number of steps used from Eq. 7
    @param constant_A: geometric series factor from Eq. 12, i.e. C_k = A
    @return: probability distribution of the realised residue distances
    """
    f_low_k = f_low(link_length, chain_length)
    f_high_k = f_high(link_length, chain_length, half_n_harmonic)
    return (((1 - f_low_k) / (1 - f_high_k)) ** constant_a) * (constant_A / f_high_k)


def distribution_of_links(f_low_k: float, f_high_k: float, chain_length: int,
                          constant_a: int, k: int) -> float:
    """
    @param f_low_k: sequence distribution evaluated at a given link_length with k << N/2
    @param f_high_k: sequence distribution evaluated at a link length with k >> 1
    @param chain_length: number of C-alphas
    @param constant_a: number of steps used from Eq. 7
    @param k: step number in Eq. 11
    @return: expected distribution of links from Eq. 11
    """
    return chain_length * (((1 - f_low_k) / (1 - f_high_k)) ** constant_a) * (1 - f_high_k) ** k


# noinspection PyShadowingNames
def plotting_command_line_arguments() -> argparse.Namespace:
    """
    Parser for argparse used for plotting # move to functions.py
    @return: Namespace containing command line arguments
    """
    parser = argparse.ArgumentParser(description="Plot link lengths and residuals.")
    parser.add_argument("link_length", type=str, choices=["100", "200", "300"], help="link length to be plotted")
    parser.add_argument("algorithm", type=str, choices=["BS", "C"],
                        help="get link lengths from bootstrapping (BS) or chunking (C)")
    parser.add_argument("--A-begin", dest="start_A", type=float, help="starting value for constant A", required=True)
    parser.add_argument("--A-end", dest="end_A", type=float, help="last value for constant A", required=True)
    parser.add_argument("--a-begin", dest="start_a", type=int, help="starting value for constant a", required=True)
    parser.add_argument("--a-end", dest="end_a", type=int, help="starting value for constant a", required=True)
    parser.add_argument("--a-step", dest="step_a", type=int, help="step size for constant a range, default=1")
    parser.add_argument("--A-step", dest="step_A", type=float, help="step size for constant A range, default=0.0001")
    parser.set_defaults(a_step=1, A_step=0.0001)

    arguments = parser.parse_args()

    return arguments


arguments = plotting_command_line_arguments()
link_length = arguments.link_length
algorithm = arguments.algorithm
start_A = arguments.start_A
end_A = arguments.end_A
start_a = arguments.start_a
end_a = arguments.end_a
step_A = arguments.step_A
step_a = arguments.step_a

dataframe = None
if algorithm == "BS":
    dataframe = pd.read_csv(f"../data/pdb_data/bootstrapped_{link_length}s_with_stats.csv")
elif algorithm == "C":
    dataframe = pd.read_csv(f"../data/alphafold_data/chunked_{link_length}s_with_stats.csv")

means = dataframe["mean"].to_numpy()
lower_b = dataframe["lower_bound"].to_numpy()
upper_b = dataframe["upper_bound"].to_numpy()
var = dataframe["variable"].to_numpy()

normed_means = normalise_data(means)
normed_lower_b = normalise_data(lower_b)
normed_upper_b = normalise_data(upper_b)

n_datapoints = int(var[-1] + 1)
half_n_harmonic = nth_harmonic(n_datapoints)
sumrange = np.array(range(int(var[0]), n_datapoints))

if start_A == end_A and start_a == end_a:
    A = start_A
    a = start_a
    link_length_fig = plt.figure()
    ax = link_length_fig.subplots(1, 1)
    ax.scatter(var, normed_means, s=10, c="k", label="PDB 200s")
    ax.fill_between(var, normed_lower_b, normed_upper_b, color="gray", alpha=1, label="95% C.L.")
    ax.plot(sumrange, [link_length_distribution(s, n_datapoints, half_n_harmonic, a, A) for s in sumrange],
            c="#984ea3", label="Theory")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylabel(f"a = {a}", fontsize=13, rotation=0, labelpad=21)
    ax.yaxis.set_label_position("right")
    ax.set_title(f"A = {A}")
    plt.legend()
    link_length_fig.text(0.5, 0.025, "s / a.u. ", ha="center", fontsize=15.5)
    link_length_fig.text(0.005, 0.5, "P(s)", va="center", rotation="vertical", fontsize=15.5)

    residual_fig = plt.figure()
    ax = residual_fig.subplots(1, 1)
    theory_function = [link_length_distribution(s, n_datapoints, half_n_harmonic, a, A) for s in sumrange]
    residuals = normed_means - theory_function
    mean_of_residuals = np.mean(residuals)
    std_of_residuals = np.std(residuals)
    sum_of_residuals = np.sum(residuals)
    r_squared = r2_score(normed_means, theory_function)
    RSS = np.sum(residuals ** 2)
    ax.scatter(sumrange, residuals, s=10, marker=".", c="red", label="Residuals")
    ax.hlines(0, sumrange[0], sumrange[-1], ls="--", color="k", label="Zero")
    ax.plot([], [], ls=' ', label=f"RSS: {RSS:.4f}")
    ax.plot([], [], ls=' ', label=f"$R^2$: {r_squared:.4f}")
    ax.plot([], [], ls=' ', label=f"<R>: {mean_of_residuals:.4f}")
    ax.set_ylabel(f"a = {a}", fontsize=13, rotation=0, labelpad=21)
    ax.yaxis.set_label_position("right")
    ax.set_title(f"A = {A}")
    ax.legend()
    residual_fig.text(0.5, 0.025, "s / a.u. ", ha="center", fontsize=15.5)
    residual_fig.text(0.005, 0.5, "Residuals", va="center", rotation="vertical", fontsize=15.5)
    plt.show()
else:
    A_range = np.arange(start_A, end_A, step_A)
    a_range = np.arange(start_a, end_a, step_a)
    link_length_fig = plt.figure()
    ax = link_length_fig.subplots(len(a_range), len(A_range), sharex=True, sharey=True)
    for row in range(len(a_range)):
        for col in range(len(A_range)):
            ax[row][col].scatter(var, normed_means, s=10, c="k", label="AF PDB 300s")
            ax[row][col].fill_between(var, normed_lower_b, normed_upper_b, color="gray", alpha=0.6, label="95% C.L.")
            ax[row][col].plot(sumrange, [link_length_distribution(s, n_datapoints, half_n_harmonic, a_range[row],
                                                                  A_range[col]) for s in sumrange],
                              c="#984ea3", label="Theory")
            ax[row][col].set_yscale("log")
            ax[row][col].set_xscale("log")
            ax[row][-1].set_ylabel(f"a = {a_range[row]}", fontsize=13, rotation=0, labelpad=21)
            ax[row][-1].yaxis.set_label_position("right")
            ax[0][col].set_title(f"A = {A_range[col]:2f}")
    plt.legend()
    link_length_fig.text(0.5, 0.025, "s / a.u. ", ha="center", fontsize=15.5)
    link_length_fig.text(0.005, 0.5, "P(s)", va="center", rotation="vertical", fontsize=15.5)
    plt.subplots_adjust(left=0.06, bottom=0.08, top=0.95, wspace=0.1, right=0.95)

    fig = plt.figure()
    ax = fig.subplots(len(a_range), len(A_range), sharex=True, sharey=True)
    for row in range(len(a_range)):
        for col in range(len(A_range)):
            theory_function = [link_length_fig(s, n_datapoints, half_n_harmonic, a_range[row], A_range[col])
                               for s in sumrange]
            residuals = normed_means - theory_function
            mean_of_residuals = np.mean(residuals)
            std_of_residuals = np.std(residuals)
            sum_of_residuals = np.sum(residuals)
            r_squared = r2_score(normed_means, theory_function)
            RSS = np.sum(residuals ** 2)
            ax[row][col].scatter(sumrange, residuals, s=10, marker=".", c="red", label="Residuals")
            ax[row][col].hlines(0, sumrange[0], sumrange[-1], ls="--", color="k", label="Zero")
            ax[row][col].plot([], [], ls=" ", label=f"RSS: {RSS:.4f}")
            ax[row][col].plot([], [], ls=" ", label=f'$R^2$: {r_squared:.4f}')
            ax[row][col].plot([], [], ls=" ", label=f"<R>: {mean_of_residuals:.4f}")
            ax[row][-1].set_ylabel(f"a = {a_range[row]}", fontsize=13, rotation=0, labelpad=21)
            ax[row][-1].yaxis.set_label_position("right")
            ax[0][col].set_title(f"A = {A_range[col]:2f}")
            ax[row][col].legend()
    fig.text(0.5, 0.025, "s / a.u. ", ha="center", fontsize=15.5)
    fig.text(0.005, 0.5, "Residuals", va="center", rotation="vertical", fontsize=15.5)
    plt.subplots_adjust(left=0.06, bottom=0.08, top=0.95, wspace=0.1, right=0.95)
    plt.show()
