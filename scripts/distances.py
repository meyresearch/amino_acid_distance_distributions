"""Functions for plotting RCSB and AlphaFold amino acid distance distributions for each length in same plot"""
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import plot_functions


def create_plot_label(arguments: argparse.Namespace, length_range: str) -> str:
    """
    Create legend label for the plots
    @param arguments: command line arguments
    @param length_range: 100, 200 or 300
    @return:
    """
    label = ""
    if arguments.algorithm == "rcsb":
        label = f"RCSB {length_range}"
    elif arguments.algorithm == "alpha":
        label = f"AlphaFold 2 {length_range}"
    return label


def create_plots(arguments: argparse.Namespace) -> None:
    """
    Plot amino acid distances in each length range on the same plot
    @param arguments: command line arguments
    @return: None
    """
    length_ranges = ["100", "200", "300"]

    colours = ["#a3dbab", "#47b856", "#286830"]
    if arguments.algorithm == "rcsb":
        colours = ["#e2b6dc", "#b847a9", "#682860"]
    elif arguments.algorithm == "alpha":
        colours = ["#a3dbab", "#47b856", "#15371a"]
    markers = ["o", "s", "^"]
    plt.figure(figsize=(6, 6))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=1.8, font="Helvetica")
    for i in range(len(length_ranges)):
        histogram = plot_functions.get_histogram(length_ranges[i], arguments.algorithm, arguments.low)
        plotting_tuple = plot_functions.get_data_for_plotting(histogram, arguments, length_ranges[i])
        distance_bins = plotting_tuple[1]
        measure = plotting_tuple[2]
        plt.scatter(distance_bins, measure, label=create_plot_label(arguments, length_ranges[i]),
                    color=colours[i], marker=markers[i])
        plt.yscale("log")
        plt.xscale("log")
        plt.xlabel("s")
        plt.ylabel("P(s)")
        plt.xlim(2, 315)
        plt.ylim(0.0001, 0.1)
        plt.legend(loc="upper right", fontsize=20, frameon=False)
        sns.despine()
        plt.tight_layout()
    plt.show()
    if arguments.low:
        plt.savefig("../plots/supplementary_material/low_confidence_af2.jpeg", dpi=2600)
    elif not arguments.low:
        plt.savefig("../plots/supplementary_material/af2.jpeg", dpi=2600)


