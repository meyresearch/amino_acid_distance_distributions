"""
Compute amino acid distances using Shadow map [https://doi.org/10.1021/jp300852d]
"""
import argparse
import functions as fn
import plot_functions as pf
import matplotlib.pyplot as plt
import seaborn as sns
from colour_palette import _COLOUR_PALETTE


def shadow_commandline_options() -> None:
    """
    Read in CL arguments and make a choice depending on given option
    @param cl_arguments from user
    @return: None
    """
    parser = argparse.ArgumentParser(description="Get amino acid distance distributions from Shadow maps.")
    parser.add_argument("-p", "--path", type=str, help="Full path to input file")
    parser.add_argument("-c", "--cutoff", type=str, help="Cutoff distance in Angstrom")
    parser.add_argument("-s", "--shadow-radius", type=str, help="Shadowing radius")
    parser.add_argument("-d", "--distances", action="store_true", help="Compute distances from Shadow map")
    parser.add_argument("-pd", "--plot-distances", action="store_true", help="Plot Shadow map distances")
    parser.add_argument("-sg", "--smog", action="store_true", help="Compute Shadow maps from SMOG2")
    parser.add_argument("-r", "--length-range", type=str, choices=["100", "200", "300"], default="300", help="Chain length range")
    parser.add_argument("-m", "--measure", type=str, choices=["mean", "median"], default="mean", help="Measure of central tendency for plotting")
    parser.add_argument("-q", "--quantile", type=int, choices=[1, 2, 3], default=2, help="How many standard deviations to use in CL")
    return parser.parse_args()


def plot_shadow_comparison(cutoff: str, length_range: str, arguments: argparse.Namespace):
    """
    """
    path_s0 = f"../data/alphafold/shadow_distance_histogram_not_normed_s_0_c_{cutoff}.npy"
    path_s1 = f"../data/alphafold/shadow_distance_histogram_not_normed_s_1_c_{cutoff}.npy"
    histogram_s0 = pf.get_histogram("", algorithm=path_s0, low_confidence=False)
    histogram_s1 = pf.get_histogram("", algorithm=path_s1, low_confidence=False)
    pdb_path  = f"../data/rcsb/histogram_c_{cutoff}_{length_range}_not_normed.npy"
    pdb_histogram = pf.get_histogram("", algorithm=pdb_path, low_confidence=False)
    chain_length, bins, n_mean_pdb, n_low_bound, n_high_bound = pf.get_data_for_plotting(pdb_histogram, arguments, length_range)
    _, _, n_mean_s0, _, _ = pf.get_data_for_plotting(histogram_s0, arguments, length_range)
    _, _, n_mean_s1, _, _ = pf.get_data_for_plotting(histogram_s1, arguments, length_range)
    fig = plt.figure()
    fig.set_size_inches((6, 6))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=1.88)
    plt.scatter(bins, n_mean_pdb, label=f"RCSB 300", marker="o", color=_COLOUR_PALETTE["PDB_SCATTER"])
    plt.fill_between(bins, n_high_bound, n_low_bound, color=_COLOUR_PALETTE["PDB_SCATTER"], alpha=0.25, label="RCSB 300 95% C.L.", zorder=-100)
    plt.scatter(bins, n_mean_s0, label=f"AF2 SCM S = 0 \u212B", marker="^", color="#47b856")
    plt.scatter(bins, n_mean_s1, label=f"AF2 SCM S = 1 \u212B", marker="s", color="#15371a")
    plt.yscale("log")
    plt.xscale("log")
    plt.xlim(4, 150)
    plt.ylim(0.0005, 0.1)
    plt.xlabel("s")
    plt.ylabel("P(s)")
    plt.legend(loc="upper right", fontsize=16)
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"../plots/supplementary_material/shadow_comparison_c_{cutoff}.pdf")
    plt.show()


def main():
    arguments = shadow_commandline_options()
    path = arguments.path
    cutoff_radius = arguments.cutoff
    shadow_radius = arguments.shadow_radius
    distances = arguments.distances
    plot = arguments.plot_distances
    smog = arguments.smog
    length_range = arguments.length_range
    if smog:
        print("Computing Shadow maps from SMOG2")
        fn.run_smog(path, cutoff_radius, shadow_radius)
    if distances:
        print("Computing distances from Shadow map")
        fn.get_shadow_distance_histograms(path, cutoff_radius, shadow_radius)
    if plot:
        print("Plotting distances from Shadow map")
        plot_shadow_comparison(cutoff_radius, length_range, arguments)

        
if __name__ == "__main__":
    main()
