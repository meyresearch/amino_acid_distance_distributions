import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import pandas as pd
import scipy
from colour_palette import _COLOUR_PALETTE


def contour_cloud(alphas, betas, colour, plot_alpha) -> matplotlib.contour.ContourSet:
    """
    Take alphas and betas and return a contour cloud in an alpha-beta plot
    @param alphas: alpha helix content
    @param betas: beta sheet content
    @param colour: colour to use for the contour plot
    @param plot_alpha:
    @return:
    """
    x_min, x_max = 0, 1
    y_min, y_max = 0, 1

    xx, yy = np.mgrid[x_min:x_max:100j, y_min:y_max:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([alphas, betas])
    kernel = scipy.stats.gaussian_kde(values)
    z = np.reshape(kernel(positions).T, xx.shape)
    contour = plt.contour(xx, yy, z, colors=colour, alpha=plot_alpha)
    return contour


def create_contour_plots() -> None:
    """
    Create contour distribution plots of the secondary structure content in PDB files
    @return: None
    """
    lengths = ["100", "200", "300"]
    indices = [50, 50, 70]
    for i in range(len(lengths)):
        sns.set(context="notebook", style="ticks", font_scale=1.8, font="Helvetica")
        alphafold_df = pd.read_csv(f"../data/alphafold/unique_secondary_structures_{lengths[i]}.csv")
        rcsb_df = pd.read_csv(f"../data/rcsb/secondary_structures_{lengths[i]}.csv")
        rcsb_df = rcsb_df.loc[(rcsb_df[['alpha', 'beta']] != 0).all(axis=1)]
        alphafold_df = alphafold_df.loc[(alphafold_df[["alpha", "beta"]] != 0).all(axis=1)]
        rcsb_x = rcsb_df["alpha"]
        rcsb_y = rcsb_df["beta"]
        alphafold_x = alphafold_df["alpha"]
        alphafold_y = alphafold_df["beta"]
        plt.scatter(x=rcsb_x[::indices[i]], y=rcsb_y[::indices[i]],
                    color=_COLOUR_PALETTE["RCSB_CONTOUR"], alpha=0.1)
        plt.scatter(x=alphafold_x[::indices[i]], y=alphafold_y[::indices[i]],
                    color=_COLOUR_PALETTE["ALPHAFOLD_CONTOUR"], alpha=0.3)
        fig = plt.gcf()
        ax = plt.gca()
        fig.set_size_inches((5, 5))
        rcsb_contour = contour_cloud(alphas=rcsb_df["alpha"], betas=rcsb_df["beta"],
                                     colour=_COLOUR_PALETTE["RCSB_CONTOUR"], plot_alpha=0.6)
        alphafold_contour = contour_cloud(alphas=alphafold_df["alpha"], betas=alphafold_df["beta"],
                                          colour=_COLOUR_PALETTE["ALPHAFOLD_CONTOUR"], plot_alpha=0.7)
        rcsb_contour.collections[0].set_label(f"RCSB {lengths[i]}")
        alphafold_contour.collections[0].set_label(f"AlphaFold {lengths[i]}")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel(r"$\alpha$-helix %")
        ax.set_ylabel(r"$\beta$-sheet %")
        plt.xticks(np.arange(0, 1.0, 0.2))
        ax.legend()
        sns.despine()
        plt.tight_layout()
        plt.show()
