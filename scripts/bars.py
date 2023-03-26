"""Functions for plotting bar plots"""
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from colour_palette import _COLOUR_PALETTE


def get_data_for_bars(path_to_data: str) -> tuple:
    """
    Read csv file from path and return numpy arrays as tuple
    @param path_to_data:
    @return:
    """
    dataframe = pd.read_csv(path_to_data)
    bins = dataframe["Bins"].to_numpy()
    frequencies = dataframe["Number"].to_numpy()
    return bins, frequencies


def calculate_adjusted_frequency(frequency_to_adjust: np.ndarray, bottom_frequency: np.ndarray) -> np.ndarray:
    """
    Calculate the adjusted frequency of the "top" bar to create a stack bar plot
    @param frequency_to_adjust: "top" bars i.e. either all SwissProt or all RCSB frequencies
    @param bottom_frequency: "bottom" bars i.e. the used PDBs' frequencies
    @return: difference between top and bottom bars
    """
    return frequency_to_adjust - bottom_frequency


def create_bar_plots() -> None:
    """
    Bring all stats together and plot subplots of PDB statistics
    @return: None
    """
    path = "../data/pdb_statistics/"
    alphafold_file = "alphafold_unique.csv"
    pdb_file = "RCSB_used.csv"
    swissprot_file = "UniProt_Swiss_Prot.csv"
    rcsb_file = "RCSB_by_length.csv"

    pdb_bins, pdb_frequencies = get_data_for_bars(path + pdb_file)
    bins, rcsb_frequencies = get_data_for_bars(path + rcsb_file)

    adjusted_rcsb = calculate_adjusted_frequency(rcsb_frequencies, pdb_frequencies)

    fig, ax = plt.subplots(figsize=(6, 6))

    ax.bar(bins, adjusted_rcsb,
              color=_COLOUR_PALETTE["DATABANK"],
              label="RCSB PDB",
              bottom=pdb_frequencies)
    ax.bar(bins, pdb_frequencies,
              color=_COLOUR_PALETTE["USED"],
              label="Used RCSB PDB")

    ax.tick_params(axis="x", labelrotation=90, labelsize=20)
    ax.tick_params(axis="y", labelsize=20)
    ax.legend(fontsize=18, frameon=False)
    ax.set_ylabel("Frequency\n", fontsize=24)
    ax.set_xlabel("Chain length", fontsize=24)
    plt.tight_layout()
    sns.despine()
    plt.savefig("../plots/individual_plots_for_paper/bars.pdf")
