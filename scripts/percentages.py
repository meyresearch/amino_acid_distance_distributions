"""
Plot alpha helix content vs beta sheet content
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib


def plot_percentages(dataframe: pd.DataFrame, use_color_map: bool) -> None:
    """
    Take secondary structure information and plot alpha percentage vs beta percentage
    @param dataframe: dataframe containing RCSB or AlphaFold data
    @param use_color_map: True for AlphaFold, False for RCSB
    @return: None
    """
    plt.figure(figsize=(8, 8))
    sns.set(context="notebook", palette="colorblind", style='ticks', font_scale=1.8, font='Helvetica')
    if use_color_map:
        axes = plt.scatter(dataframe["H"], dataframe["beta"], c=dataframe["mean"], cmap="magma_r")
        plt.colorbar(axes)
    else:
        plt.scatter(dataframe["H"], dataframe["E"]+dataframe["B"], c="#fbafe4")
    plt.xlabel(r"$\alpha$-helix %")
    plt.ylabel(r"$\beta$-sheet %")
    sns.despine()


length_ranges = ["100", "200", "300"]
for length in length_ranges:
    alphafold = pd.read_csv(f"../data/alphafold/secondary_structures_{length}.csv")
    plot_percentages(alphafold, use_color_map=True)
    plt.savefig(f"../plots/alphafold/percentages_{length}.pdf")
    plt.show()
    rcsb = pd.read_csv(f"../data/rcsb/structures_{length}_raw.csv")
    plot_percentages(rcsb, use_color_map=False)
    plt.savefig(f"../plots/rcsb/percentages_{length}.pdf")
    plt.show()
