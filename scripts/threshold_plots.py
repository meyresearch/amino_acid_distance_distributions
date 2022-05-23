""" Supplemental plots for investigating the effect of the threshold value on amino acid distance
    distributions.
"""

import functions
import numpy as np
import pandas as pd
import protein_contact_map
import matplotlib.pyplot as plt
from colour_palette import _COLOUR_PALETTE
import seaborn as sns


lengths = [100, 200, 300]

for length in lengths:
    default_threshold = np.load(f"../data/rcsb/histogram_{length}_not_normed.npy", allow_pickle=True)
    simulation = np.load(f"../data/simulations/3d/histogram_{length}_not_normed.npy", allow_pickle=True)
    threshold_10 = np.load(f"../data/rcsb/threshold_histograms/histogram_{length}_10.npy", allow_pickle=True)
    threshold_15 = np.load(f"../data/rcsb/threshold_histograms/histogram_{length}_15.npy", allow_pickle=True)
    threshold_21 = np.load(f"../data/rcsb/threshold_histograms/histogram_{length}_21.npy", allow_pickle=True)
    distances = np.linspace(1, 350, 350)[:-1]

    mean_default = np.mean(default_threshold, axis=0)
    normed_mean_default = mean_default / np.sum(mean_default)
    mean_10 = np.mean(threshold_10, axis=0)
    normed_mean_10 = mean_10 / np.sum(mean_10)
    mean_15 = np.mean(threshold_15, axis=0)
    normed_mean_15 = mean_15 / np.sum(mean_15)
    mean_21 = np.mean(threshold_21, axis=0)
    normed_mean_21 = mean_21 / np.sum(mean_21)
    simulation_mean = np.mean(simulation, axis=0)
    normed_simulation_mean = simulation_mean / np.sum(simulation_mean)

    fig = plt.figure()
    fig.set_size_inches(8, 8)
    sns.set(context="notebook", palette="colorblind", style="ticks", font_scale=1.88, font="Helvetica")
    plt.scatter(distances, normed_simulation_mean, label=f"3D Simulation {length}", c="gray")
    plt.scatter(distances, normed_mean_default, label=f"Threshold 8 Å {length}",
                c=_COLOUR_PALETTE["PDB_SCATTER"], marker="s")
    plt.scatter(distances, normed_mean_10, label=f"Threshold 10 Å {length}", c="#e2b6dc", marker="^")
    plt.scatter(distances, normed_mean_15, label=f"Threshold 15 Å {length}", c="#b847a9", marker="v")
    plt.scatter(distances, normed_mean_21, label=f"Threshold 21 Å {length}", c="#286830", marker="x")
    sns.despine()
    plt.xlim(2, int(length))
    plt.ylim(0.0001, 0.1)
    plt.loglog()
    plt.tight_layout()
    plt.xlabel("s")
    plt.ylabel("P(s)")
    plt.legend(frameon=False)
    plt.show()
