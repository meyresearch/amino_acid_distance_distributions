"""
Get amino acid distance distributions for a user-specified value of C
"""
import argparse
import functions
import plot_functions
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from colour_palette import _COLOUR_PALETTE


def commandline_arguments():
    parser = argparse.ArgumentParser(description="Get amino acid distance distributions for a user-specified value of C")
    parser.add_argument("-a", "--algorithm", type=str, choices=["rcsb", "alpha"], default="alpha")
    parser.add_argument("-c", "--cutoff", type=str, help="Cutoff value")
    parser.add_argument("-r", "--range", type=str, choices=["100", "200", "300"], default="300", help="Chain length range")
    parser.add_argument("-p", "--path", type=str, help="Path to csv file containing PDB data")
    parser.add_argument("-d", "--distances", action="store_true", help="Compute distances")
    parser.add_argument("-pd", "--plot-distances", type=str, help="npy file for plotting distances")
    parser.add_argument("-q", "--quantile", type=int, nargs="?", choices=[1, 2, 3], const=2, default=2, help="confidence interval in multipltes of std.")
    parser.add_argument("-m", "--measure", type=str, choices=["median", "mean"], nargs="?", const="mean", default="mean", help="measure of central tendency")

    return parser.parse_args()
    

def create_plot(arguments: argparse.Namespace, histogram: np.array):
    
    plotting_tuple = plot_functions.get_data_for_plotting(histogram, arguments, arguments.range)
    chain_length = plotting_tuple[0]
    distances = plotting_tuple[1]
    mean = plotting_tuple[2]
    lower_bound, upper_bound = plotting_tuple[3], plotting_tuple[4]
    x1, x2 = 4, 150
    fig = plt.figure()
    fig.set_size_inches((6, 6))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=1.88)
    plt.scatter(distances, mean, label=f"RCSB PDB {chain_length} c={arguments.cutoff} Å",
                color=_COLOUR_PALETTE["PDB_SCATTER"], marker="o")
    plt.fill_between(distances, upper_bound, lower_bound,
                color=_COLOUR_PALETTE["PDB_SCATTER"], alpha=0.25, label="RCSB 95% C.L.", zorder=-99)
    plt.xlim(x1, x2)
    plt.ylim(0.0005, 0.1)
    plt.xlabel("s")
    plt.ylabel("P(s)")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(loc="upper right", fontsize=16)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"../plots/supplementary_material/distances_c_{arguments.cutoff}.png")


def main():
    arguments = commandline_arguments()
    algorithm = arguments.algorithm
    cutoff = arguments.cutoff
    length_range = arguments.range
    path = arguments.path
    distances = arguments.distances
    plot_path = arguments.plot_distances

    if distances:
        functions.get_distances_with_different_cutoff(algorithm, length_range, path, cutoff)
    if plot_path:
        histogram = plot_functions.get_histogram(length_range, algorithm=plot_path, low_confidence=False)
        create_plot(arguments, histogram)


if __name__ == "__main__":
    main()