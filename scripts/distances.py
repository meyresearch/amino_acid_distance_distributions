"""Functions for plotting RCSB and AlphaFold amino acid distance distributions"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import seaborn as sns
import theory_functions
import plot_functions
from colour_palette import _COLOUR_PALETTE
import pandas as pd
import argparse


def create_plot_label(arguments: argparse.Namespace, length_range: str) -> str:
    """

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
    Plot amino acid distances and residuals in single plots
    @param arguments: command line arguments
    @return: None
    """
    length_ranges = ["100", "200", "300"]
    colours = ["#286830", "#682860", "#fbafe4"]
    plt.figure(figsize=(8, 8))
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=2.4, font="Helvetica")
    for i in range(len(length_ranges)):
        histogram = plot_functions.get_histogram(length_ranges[i], arguments.algorithm)
        plotting_tuple = plot_functions.get_data_for_plotting(histogram, arguments, length_ranges[i])
        distance_bins = plotting_tuple[1]
        measure = plotting_tuple[2]
        plt.scatter(distance_bins, measure, label=create_plot_label(arguments, length_ranges[i]), color=colours[i])
        plt.yscale("log")
        plt.xscale("log")
        plt.xlabel("s")
        plt.ylabel("P(s)")
        plt.xlim(2, 315)
        plt.legend()
        sns.despine()
        plt.tight_layout()
    plt.show()
