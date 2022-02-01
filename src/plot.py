"""
plot.py
"""
import argparse
import glob
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
from numpy import ndarray
from sklearn.metrics import r2_score
_COLOUR_PALETTE = {"PDB_SCATTER_COLOUR": "#009E73", "SIM_SCATTER_COLOUR": "#E69F00", "CL_COLOUR": "#56B4E9",
                   "THEORY_COLOUR": "k", "RESIDUALS_COLOUR": "#009E73"}


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


def link_length_distribution(link_length_s: int, chain_length: int, half_n_harmonic_number: float, constant_a: int,
                             constant_A: float) -> float:
    """

    @param link_length_s: separation between each link 2 <= s < N/2
    @param chain_length: number of C-alphas
    @param half_n_harmonic_number: the 'N/2'-th harmonic number
    @param constant_a: number of steps used from Eq. 7
    @param constant_A: geometric series factor from Eq. 12, i.e. C_k = A
    @return: probability distribution of the realised residue distances
    """
    f_low_k = f_low(link_length_s, chain_length)
    f_high_k = f_high(link_length_s, chain_length, half_n_harmonic_number)
    return (((1 - f_low_k) / (1 - f_high_k)) ** constant_a) * (constant_A / f_high_k)


def distribution_of_links(f_low_k: float, f_high_k: float, chain_length: int, constant_a: int, k: int) -> float:
    """
    @param f_low_k: sequence distribution evaluated at a given link_length with k << N/2
    @param f_high_k: sequence distribution evaluated at a link length with k >> 1
    @param chain_length: number of C-alphas
    @param constant_a: number of steps used from Eq. 7
    @param k: step number in Eq. 11
    @return: expected distribution of links from Eq. 11
    """
    return chain_length * (((1 - f_low_k) / (1 - f_high_k)) ** constant_a) * (1 - f_high_k) ** k


def get_dataframe(command_line_arguments: argparse.Namespace) -> pd.DataFrame:
    """
    Get link length dataframe for plotting
    @param command_line_arguments:
    @return: dataframe containing bootstrapped/chunked data and confidence level upper and lower limits
    """
    link_length_dataframe = None
    linklength = command_line_arguments.link_length
    if command_line_arguments.algorithm == "BS":
        link_length_dataframe = pd.read_csv(f"../data/rcsb_data/bs_{linklength}_stats.csv")
    elif command_line_arguments.algorithm == "C":
        link_length_dataframe = pd.read_csv(f"../data/alphafold_data/chunked_{linklength}s_with_stats.csv")
    return link_length_dataframe


def get_pdb_data_for_plotting(dataframe_for_plotting: pd.DataFrame) -> tuple:
    """
    Get data and confidence level limits from link length dataframe
    @param dataframe_for_plotting:
    @return: tuple containing number of datapoints, normalised means, and lower and upper confidence level limits
    """
    means = dataframe_for_plotting["mean"].to_numpy()
    lower_b = dataframe_for_plotting["lower_bound"].to_numpy()
    upper_b = dataframe_for_plotting["upper_bound"].to_numpy()
    variables = dataframe_for_plotting["variable"].to_numpy()
    number_of_datapoints = int(variables[-1] + 1)
    return number_of_datapoints, normalise_data(means), normalise_data(lower_b), normalise_data(upper_b)


def get_simulation_filename(command_line_length: str) -> str:
    """

    @param command_line_length:
    @return:
    """
    length = " "
    if command_line_length == "100":
        length = command_line_length
    elif command_line_length == "200":
        length = "209"
    elif command_line_length == "300":
        length = "302"
    return length


def get_simulation_files(command_line_length: str) -> list:
    """

    @param command_line_length:
    @return:
    """
    file_link_length = get_simulation_filename(command_line_length)
    simulation_files = glob.glob(f"../data/simulation_data/Nora50_matrix_{file_link_length}_*")
    return simulation_files


def get_simulation_link_lengths(command_line_link_length: str) -> tuple:
    """

    @param command_line_link_length:
    @return:
    """
    files = get_simulation_files(command_line_link_length)
    link_length_list = []
    histogram, edges = None, None
    for file in files:
        matrix = np.loadtxt(file)
        matrix[matrix > 1] = 0
        protein_graph = nx.from_numpy_matrix(matrix)
        links = list(protein_graph.edges())
        for link in links:
            link_length_list.append(abs(link[0] - link[1]))
        max_link_length = max(link_length_list)
        histogram, edges = np.histogram(link_length_list, bins=range(2, max_link_length), density=True)
    # normalised_histogram = histogram / len(files)
    return histogram, edges


# noinspection PyPep8Naming
def parse_command_line_arguments() -> argparse.Namespace:
    """
    Parse command line arguments for plotting # move to plotting_functions.py
    @return: Namespace containing command line arguments
    """
    parser = argparse.ArgumentParser(description="Plot link lengths and residuals.")
    parser.add_argument("link_length", type=str, choices=["100", "200", "300"], help="link length to be plotted")
    parser.add_argument("algorithm", type=str, choices=["BS", "C", "BOTH"],
                        help="get link lengths from bootstrapping (BS), chunking (C) or compare both")
    parser.add_argument("--A-begin", dest="start_A", type=float, help="starting value for constant A")
    parser.add_argument("--A-end", dest="end_A", type=float, help="last value for constant A")
    parser.add_argument("--a-begin", dest="start_a", type=int, help="starting value for constant a")
    parser.add_argument("--a-end", dest="end_a", type=int, help="starting value for constant a")
    parser.add_argument("--a-step", dest="step_a", type=int, nargs="?", const=1, default=1,
                        help="step size for constant a range, default=1")
    parser.add_argument("--A-step", dest="step_A", type=float, nargs="?", const=0.0001, default=0.0001,
                        help="step size for constant A range, default=0.0001")
    parser.add_argument("--A-BS", dest="A_bs", type=float, help="constant A for BS algorithm in comparison plot")
    parser.add_argument("--A-C", dest="A_c", type=float, help="constant A for C algorithm in comparison plot")
    parser.add_argument("--a-BS", dest="a_bs", type=float, help="constant a for BS algorithm in comparison plot")
    parser.add_argument("--a-C", dest="a_c", type=float, help="constant a for C algorithm in comparison plot")

    command_line_arguments = parser.parse_args()
    check_algorithm = command_line_arguments.algorithm
    check_start_A = command_line_arguments.start_A
    check_end_A = command_line_arguments.end_A
    check_start_a = command_line_arguments.start_a
    check_end_a = command_line_arguments.end_a
    check_A_BS = command_line_arguments.A_bs
    check_A_C = command_line_arguments.A_c
    check_a_BS = command_line_arguments.a_bs
    check_a_C = command_line_arguments.a_c
    if (check_algorithm == "BS" and check_start_A is None) or (check_algorithm == "BS" and check_end_A is None) or \
            (check_algorithm == "BS" and check_start_a is None) or (check_algorithm == "BS" and check_end_a is None):
        parser.error("BS requires --A-begin, --A-end, --a-begin, --a-end")
    elif (check_algorithm == "C" and check_start_A is None) or (check_algorithm == "C" and check_end_A is None) or \
            (check_algorithm == "C" and check_start_a is None) or (check_algorithm == "C" and check_end_a is None):
        parser.error("C requires --A-begin, --A-end, --a-begin, --a-end")
    elif (check_algorithm == "BOTH" and check_A_BS is None) or (check_algorithm == "BOTH" and check_A_C is None) \
            or (check_algorithm == "BOTH" and check_a_BS is None) or (check_algorithm == "BOTH" and
                                                                      check_a_C is None):
        parser.error("BOTH requires --A-BS, --A-C, --a-BS, --a-C")

    return command_line_arguments


def get_command_line_arguments(command_line_arguments: argparse.Namespace) -> tuple:
    """

    @param command_line_arguments:
    @return:
    """
    return command_line_arguments.link_length, command_line_arguments.algorithm, command_line_arguments.start_A, \
           command_line_arguments.end_A, command_line_arguments.start_a, command_line_arguments.end_a, \
           command_line_arguments.step_A, command_line_arguments.step_a, command_line_arguments.A_bs, \
           command_line_arguments.a_bs, command_line_arguments.A_c, command_line_arguments.a_c


def create_link_length_plot_label(command_line_link_length: str, command_line_algorithm: str) -> str:
    """
    Create custom labels for Alphafold PDB data and RCSB PDB data with the given link length
    @param command_line_link_length:
    @param command_line_algorithm:
    @return:
    """
    plot_label = " "
    if command_line_algorithm == "BS" or command_line_algorithm == "BOTH":
        plot_label = f"PDB {command_line_link_length}s"
    elif command_line_algorithm == "C":
        plot_label = f"AF PDB {command_line_link_length}"
    return plot_label


def get_stats_for_plotting(constant_A: float, constant_a: int, n_points: int, half_n_harmonic_number: float,
                           plotting_sumrange: np.ndarray, normed_means: np.ndarray) -> tuple:
    """

    @param constant_A:
    @param constant_a:
    @param n_points:
    @param half_n_harmonic_number:
    @param plotting_sumrange:
    @param normed_means:
    @return: tuple containing array of residuals, their mean, R^2 statistic, and the residual sum of squares
    """
    theory_function = [link_length_distribution(s, n_points, half_n_harmonic_number, constant_a, constant_A)
                       for s in plotting_sumrange]
    link_length_residuals = normed_means - theory_function
    residuals_mean = np.mean(link_length_residuals)
    residuals_sum = np.sum(link_length_residuals)
    r_square_value = r2_score(normed_means, theory_function)
    residual_sum_of_squares = residuals_sum ** 2

    return link_length_residuals, residuals_mean, r_square_value, residual_sum_of_squares


def plot_single_link_length_plot(constant_A: float, constant_a: int, n_points: int, half_n_harmonic_number: float,
                                 plotting_sumrange: np.ndarray, normed_means: np.ndarray, variables: np.ndarray,
                                 bin_centres: np.ndarray, normed_sim_data: np.ndarray, lower_conf_lev: np.ndarray,
                                 upper_conf_lev: np.ndarray, command_line_link_length: str,
                                 command_line_algorithm: str) -> None:
    """
    @param normed_sim_data:
    @param bin_centres:
    @param variables:
    @param constant_A:
    @param constant_a:
    @param half_n_harmonic_number:
    @param plotting_sumrange:
    @param n_points:
    @param upper_conf_lev:
    @param lower_conf_lev:
    @param normed_means:
    @param command_line_link_length:
    @param command_line_algorithm:
    @return:
    """
    plot_label = create_link_length_plot_label(command_line_link_length, command_line_algorithm)
    fig = plt.figure()
    ax = fig.subplots(1, 1)
    ax.scatter(variables, normed_means, s=15, c=_COLOUR_PALETTE["PDB_SCATTER_COLOUR"], label=plot_label)
    ax.scatter(bin_centres, normed_sim_data, s=15, marker="^", c=_COLOUR_PALETTE["SIM_SCATTER_COLOUR"],
               label=f"SIM {command_line_link_length}s",
               zorder=-50)
    ax.fill_between(variables, lower_conf_lev, upper_conf_lev, color=_COLOUR_PALETTE["CL_COLOUR"], alpha=0.6,
                    label=r"$95\%$ C.L.",
                    zorder=-100)
    ax.plot(plotting_sumrange, [link_length_distribution(s, n_points, half_n_harmonic_number, constant_a, constant_A)
                                for s in plotting_sumrange], c=_COLOUR_PALETTE["THEORY_COLOUR"], label="Theory")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("s", fontsize=14)
    ax.set_ylabel("P(s)", fontsize=14)
    ax.set_ylim(0.0001, 0.5)
    plt.legend(fontsize=14)
    plt.subplots_adjust(left=0.075, bottom=0.08, top=0.95, wspace=0.05, right=0.98)


def compare_link_length_plots(command_line_link_length: str, command_line_algorithm: str, bin_centres: np.ndarray,
                              normed_sim_data: np.ndarray, bs_variables: np.ndarray, bs_means: np.ndarray,
                              bs_lower_conf_lev: np.ndarray, bs_upper_conf_lev: np.ndarray,
                              c_variables: np.ndarray, c_means: np.ndarray, c_lower_conf_lev: np.ndarray,
                              c_upper_conf_lev: np.ndarray, bs_plotting_sumrange: np.ndarray, bs_n_points: int,
                              bs_half_n_harmonic_number: float, bs_constant_a: int, bs_constant_A: float,
                              c_plotting_sumrange: np.ndarray, c_n_points: int, c_half_n_harmonic_number: float,
                              c_constant_a: int, c_constant_A: float) -> None:
    plot_label = create_link_length_plot_label(command_line_link_length, command_line_algorithm)
    fig = plt.figure()
    ax = fig.subplots(1, 1)
    ax.scatter(bs_variables, bs_means, s=15, c=_COLOUR_PALETTE["PDB_SCATTER_COLOUR"], label=plot_label)
    ax.scatter(c_variables, c_means, s=15, c="#0071b2", marker="s", label="AF " + plot_label)
    ax.scatter(bin_centres, normed_sim_data, s=15, marker="^", c=_COLOUR_PALETTE["SIM_SCATTER_COLOUR"],
               label=f"SIM {command_line_link_length}s",
               zorder=-50)
    ax.fill_between(bs_variables, bs_lower_conf_lev, bs_upper_conf_lev, color=_COLOUR_PALETTE["PDB_SCATTER_COLOUR"],
                    alpha=0.6, label=r"$95\%$ C.L.", zorder=-100)
    ax.fill_between(c_variables, c_lower_conf_lev, c_upper_conf_lev, color="#0071b2", alpha=0.6,
                    label=r"AF $95\%$ C.L.", zorder=-100)
    ax.plot(bs_plotting_sumrange, [link_length_distribution(s, bs_n_points, bs_half_n_harmonic_number, bs_constant_a,
                                                            bs_constant_A) for s in bs_plotting_sumrange],
            c=_COLOUR_PALETTE["THEORY_COLOUR"], label="Theory PDB")
    ax.plot(c_plotting_sumrange, [link_length_distribution(s, c_n_points, c_half_n_harmonic_number, c_constant_a,
                                                           c_constant_A) for s in c_plotting_sumrange], c="gray",
            label="Theory AF")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("s", fontsize=14)
    ax.set_ylabel("P(s)", fontsize=14)
    ax.set_ylim(0.0001, 0.5)
    plt.legend(fontsize=14)
    plt.subplots_adjust(left=0.075, bottom=0.08, top=0.95, wspace=0.05, right=0.98)


def plot_link_lengths_grid(constant_A_range: np.ndarray, constant_a_range: np.ndarray, n_points: int,
                           half_n_harmonic_number: float, plotting_sumrange: np.ndarray, normed_means: np.ndarray,
                           variables: np.ndarray, bin_centres: np.ndarray, normed_sim_data: np.ndarray,
                           lower_conf_lev: np.ndarray, upper_conf_lev: np.ndarray, command_line_link_length: str,
                           command_line_algorithm: str) -> None:
    """
    @param normed_sim_data:
    @param bin_centres:
    @param constant_A_range:
    @param constant_a_range:
    @param n_points:
    @param half_n_harmonic_number:
    @param plotting_sumrange:
    @param normed_means:
    @param variables:
    @param lower_conf_lev:
    @param upper_conf_lev:
    @param command_line_link_length:
    @param command_line_algorithm:
    @return:
    """
    plot_label = create_link_length_plot_label(command_line_link_length, command_line_algorithm)
    fig = plt.figure()
    ax = fig.subplots(len(constant_a_range), len(constant_A_range), sharex=True, sharey=True)
    for row in range(len(constant_a_range)):
        for col in range(len(constant_A_range)):
            ax[row][col].scatter(variables, normed_means, s=10, c=_COLOUR_PALETTE["PDB_SCATTER_COLOUR"],
                                 label=plot_label)
            ax[row][col].fill_between(variables, lower_conf_lev, upper_conf_lev, color=_COLOUR_PALETTE["CL_COLOUR"],
                                      alpha=0.6, label="95% C.L.", zorder=-100)
            ax[row][col].scatter(bin_centres, normed_sim_data, s=10, marker="^",
                                 c=_COLOUR_PALETTE["SIM_SCATTER_COLOUR"], label=f"SIM {command_line_link_length}s",
                                 zorder=-50)

            ax[row][col].plot(plotting_sumrange, [link_length_distribution(s, n_points, half_n_harmonic_number,
                                                                           constant_a_range[row],
                                                                           constant_A_range[col])
                                                  for s in plotting_sumrange], c=_COLOUR_PALETTE["THEORY_COLOUR"],
                              label="Theory")
            ax[row][col].set_yscale("log")
            ax[row][col].set_xscale("log")
            ax[row][col].set_ylim(0.0001, 1)
            ax[row][-1].set_ylabel(f"a = {constant_a_range[row]}", fontsize=13, rotation=0, labelpad=21)
            ax[row][-1].yaxis.set_label_position("right")
            ax[0][col].set_title(f"A = {constant_A_range[col]:2f}")
    # bbox_to_anchor = [1.43, 4.65],
    plt.legend(bbox_to_anchor=[1.57, 4.65], fontsize=9)
    fig.text(0.5, 0.025, "s / a.u. ", ha="center", fontsize=15.5)
    fig.text(0.005, 0.5, "P(s)", va="center", rotation="vertical", fontsize=15.5)
    plt.subplots_adjust(left=0.06, bottom=0.08, top=0.95, wspace=0.05, right=0.909)


def plot_single_residual_plot(constant_A: float, constant_a: int, n_points: int, half_n_harmonic_number: float,
                              plotting_sumrange: np.ndarray, normed_means: np.ndarray):
    """

    @param constant_A:
    @param constant_a:
    @param n_points:
    @param half_n_harmonic_number:
    @param plotting_sumrange:
    @param normed_means:
    @return:
    """
    residuals = get_stats_for_plotting(constant_A, constant_a, n_points, half_n_harmonic_number, plotting_sumrange,
                                       normed_means)[0]
    # mean_of_residuals = get_stats_for_plotting(constant_A, constant_a, n_points, half_n_harmonic_number,
    #                                            plotting_sumrange, normed_means)[1]
    # r_squared = get_stats_for_plotting(constant_A, constant_a, n_points, half_n_harmonic_number, plotting_sumrange,
    #                                    normed_means)[2]
    # rss = get_stats_for_plotting(constant_A, constant_a, n_points, half_n_harmonic_number, plotting_sumrange,
    #                              normed_means)[3]
    fig = plt.figure()
    ax = fig.subplots(1, 1)
    ax.scatter(plotting_sumrange, residuals, marker=".", c=_COLOUR_PALETTE["RESIDUALS_COLOUR"], label="Residuals",
               zorder=10)
    ax.hlines(0, plotting_sumrange[0], plotting_sumrange[-1], color="gray", label="Zero")
    ax.set_xlabel("s", fontsize=14)
    ax.set_ylabel("P(s)", fontsize=14)
    ax.set_ylim(-0.015, 0.015)
    ax.legend(fontsize=14)
    plt.subplots_adjust(left=0.075, bottom=0.08, top=0.95, wspace=0.05, right=0.98)


def plot_residuals_grid(A_range: np.ndarray, a_range: np.ndarray, n_points: int, half_n_harmonic_number: float,
                        plotting_sumrange: np.ndarray, normed_means: np.ndarray):
    fig = plt.figure()
    ax = fig.subplots(len(a_range), len(A_range), sharex=True, sharey=True)
    for row in range(len(a_range)):
        for col in range(len(A_range)):
            residuals = get_stats_for_plotting(A_range[col], a_range[row], n_points, half_n_harmonic_number,
                                               plotting_sumrange, normed_means)[0]
            mean_of_residuals = get_stats_for_plotting(A_range[col], a_range[row], n_points, half_n_harmonic_number,
                                                       plotting_sumrange, normed_means)[1]
            r_squared = get_stats_for_plotting(A_range[col], a_range[row], n_points, half_n_harmonic_number,
                                               plotting_sumrange, normed_means)[2]
            rss = get_stats_for_plotting(A_range[col], a_range[row], n_points, half_n_harmonic_number,
                                         plotting_sumrange, normed_means)[3]
            ax[row][col].scatter(plotting_sumrange, residuals, marker=".", c=_COLOUR_PALETTE["RESIDUALS_COLOUR"],
                                 label="Residuals", zorder=10)
            ax[row][col].hlines(0, plotting_sumrange[0], plotting_sumrange[-1], color="grey", label="Zero")
            ax[row][-1].set_ylabel(f"a = {a_range[row]}", fontsize=13, rotation=0, labelpad=21)
            ax[row][-1].yaxis.set_label_position("right")
            ax[row][col].set_ylim(-0.025, 0.025)
            ax[0][col].set_title(f"A = {A_range[col]:2f}")
            textstring = "\n".join((r"$\Sigma \sigma = $" + f"{rss:.4f}",
                                    r"$R^{2} = $" + f"{r_squared:.4f}",
                                    r"$\bar{\sigma} = $" + f"{mean_of_residuals:.4f}"))

            text_box_properties = dict(boxstyle="square", facecolor="white", edgecolor="lightgrey", alpha=0.5)

            ax[row][col].text(0.05, 0.95, textstring, transform=ax[row][col].transAxes, fontsize=9,
                              verticalalignment='top', bbox=text_box_properties)

            ax[row][col].legend(fontsize=9)

    fig.text(0.5, 0.025, "s / a.u. ", ha="center", fontsize=15.5)
    fig.text(0.005, 0.5, "Residuals", va="center", rotation="vertical", fontsize=15.5)
    plt.subplots_adjust(left=0.06, bottom=0.08, top=0.95, wspace=0.1, right=0.95)


# noinspection PyPep8Naming
def create_plots() -> None:
    """

    """
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    arguments = parse_command_line_arguments()
    link_length = get_command_line_arguments(arguments)[0]
    algorithm = get_command_line_arguments(arguments)[1]

    if algorithm == "BOTH":
        bs_A = get_command_line_arguments(arguments)[8]
        bs_a = get_command_line_arguments(arguments)[9]
        c_A = get_command_line_arguments(arguments)[10]
        c_a = get_command_line_arguments(arguments)[11]

        # get both dataframes
        bs_dataframe = pd.read_csv(f"../data/rcsb_data/bs_{link_length}_stats.csv")
        c_dataframe = pd.read_csv(f"../data/alphafold_data/chunked_{link_length}s_with_stats.csv")

        # get possible link lengths ("variable") for both
        bs_possible_link_lengths = bs_dataframe["variable"].to_numpy()
        c_possible_link_lengths = c_dataframe["variable"].to_numpy()

        # get starting link length for both
        bs_starting_link_length = int(bs_possible_link_lengths[0])
        c_starting_link_length = int(c_possible_link_lengths[0])

        # get number of points for both
        bs_n_datapoints = int(bs_possible_link_lengths[-1] + 1)
        c_n_datapoints = int(c_possible_link_lengths[-1] + 1)

        # normalise means for both
        bs_means = bs_dataframe["mean"].to_numpy()
        c_means = c_dataframe["mean"].to_numpy()
        bs_normalised_means = normalise_data(bs_means)
        c_normalised_means = normalise_data(c_means)

        # get confidence level limits for both
        bs_lower_confidence_level = normalise_data(bs_dataframe["lower_bound"].to_numpy())
        c_lower_confidence_level = normalise_data(c_dataframe["lower_bound"].to_numpy())
        bs_upper_confidence_level = normalise_data(bs_dataframe["upper_bound"].to_numpy())
        c_upper_confidence_level = normalise_data(c_dataframe["upper_bound"].to_numpy())

        # get simulation data
        simulation_data, simulation_bin_edges = get_simulation_link_lengths(link_length)
        sim_bin_centres = simulation_bin_edges[:-1] + np.diff(simulation_bin_edges) / 2
        normalised_simulation_data = normalise_data(simulation_data)

        # get N/2 harmonic for both
        bs_half_n_harmonic_number = nth_harmonic(bs_n_datapoints // 2)
        c_half_n_harmonic_number = nth_harmonic(c_n_datapoints // 2)

        # get sumrange for both
        bs_sumrange = np.array(range(bs_starting_link_length, bs_n_datapoints))
        c_sumrange = np.array(range(c_starting_link_length, c_n_datapoints))

        # plot both on the same graph
        compare_link_length_plots(link_length, algorithm, sim_bin_centres, normalised_simulation_data,
                                  bs_possible_link_lengths, bs_normalised_means, bs_lower_confidence_level,
                                  bs_upper_confidence_level, c_possible_link_lengths, c_normalised_means,
                                  c_lower_confidence_level, c_upper_confidence_level, bs_sumrange, bs_n_datapoints,
                                  bs_half_n_harmonic_number, bs_a, bs_A, c_sumrange, c_n_datapoints,
                                  c_half_n_harmonic_number, c_a, c_A)
        plt.show()

    else:
        start_A = get_command_line_arguments(arguments)[2]
        end_A = get_command_line_arguments(arguments)[3]
        start_a = get_command_line_arguments(arguments)[4]
        end_a = get_command_line_arguments(arguments)[5]
        step_A = get_command_line_arguments(arguments)[6]
        step_a = get_command_line_arguments(arguments)[7]
        dataframe = get_dataframe(arguments)
        possible_link_lengths = dataframe["variable"].to_numpy()
        starting_link_length = int(possible_link_lengths[0])
        n_datapoints = get_pdb_data_for_plotting(dataframe)[0]
        normalised_means = get_pdb_data_for_plotting(dataframe)[1]
        lower_confidence_level = get_pdb_data_for_plotting(dataframe)[2]
        upper_confidence_level = get_pdb_data_for_plotting(dataframe)[3]

        simulation_data, simulation_bin_edges = get_simulation_link_lengths(link_length)
        sim_bin_centres = simulation_bin_edges[:-1] + np.diff(simulation_bin_edges) / 2
        normalised_simulation_data = normalise_data(simulation_data)
        half_n_harmonic = nth_harmonic(n_datapoints // 2)
        sumrange = np.array(range(starting_link_length, n_datapoints))

        if start_A == end_A and start_a == end_a:
            A, a = start_A, start_a
            plot_single_link_length_plot(A, a, n_datapoints, half_n_harmonic, sumrange, normalised_means,
                                         possible_link_lengths, sim_bin_centres, normalised_simulation_data,
                                         lower_confidence_level, upper_confidence_level, link_length, algorithm)

            plot_single_residual_plot(A, a, n_datapoints, half_n_harmonic, sumrange, normalised_means)
            plt.show()
        elif start_A != end_A and start_a != end_a:
            A_range, a_range = np.arange(start_A, end_A, step_A), np.arange(start_a, end_a, step_a)
            plot_link_lengths_grid(A_range, a_range, n_datapoints, half_n_harmonic, sumrange, normalised_means,
                                   possible_link_lengths, sim_bin_centres, normalised_simulation_data,
                                   lower_confidence_level, upper_confidence_level, link_length, algorithm)
            plot_residuals_grid(A_range, a_range, n_datapoints, half_n_harmonic, sumrange, normalised_means)
            plt.show()


def main():
    create_plots()


if __name__ == '__main__':
    main()
