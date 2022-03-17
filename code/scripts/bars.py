import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sequence_distance_distribution.code.scripts.plot_functions import _COLOUR_PALETTE
import numpy as np


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
    alphafold_file = "AlphaFold_used.csv"
    pdb_file = "ll_used.csv"
    swissprot_file = "UniProt_Swiss_Prot.csv"
    rcsb_file = "RCSB_by_length.csv"

    alpha_bins, alpha_frequencies = get_data_for_bars(path + alphafold_file)
    swiss_bins, swiss_frequencies = get_data_for_bars(path + swissprot_file)
    pdb_bins, pdb_frequencies = get_data_for_bars(path + pdb_file)
    bins, rcsb_frequencies = get_data_for_bars(path + rcsb_file)

    adjusted_swiss = calculate_adjusted_frequency(swiss_frequencies, alpha_frequencies)
    adjusted_rcsb = calculate_adjusted_frequency(rcsb_frequencies, pdb_frequencies)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.subplots(2, 1, sharex=True)

    ax[0].bar(bins, adjusted_swiss,
              color=_COLOUR_PALETTE["DATABANK"],
              label="Swiss-Prot Frequencies",
              bottom=alpha_frequencies)
    ax[0].bar(bins, alpha_frequencies,
              color=_COLOUR_PALETTE["USED"],
              label="Used AlphaFold Frequencies")
    ax[1].bar(bins, adjusted_rcsb,
              color=_COLOUR_PALETTE["DATABANK"],
              label="RCSB Frequencies",
              bottom=pdb_frequencies)
    ax[1].bar(bins, pdb_frequencies,
              color=_COLOUR_PALETTE["USED"],
              label="Used RCSB Frequencies")

    ax[0].tick_params(axis="x", labelrotation=90)
    ax[1].tick_params(axis="x", labelrotation=90)

    ax[0].legend(fontsize=14)
    ax[1].legend(fontsize=14)

    plt.subplots_adjust(left=0.09, bottom=0.08, top=0.99, wspace=0.05, right=1)
    plt.tight_layout()
    sns.despine()
    plt.show()
