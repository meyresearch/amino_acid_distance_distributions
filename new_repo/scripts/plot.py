"""
Handle all plotting functions and actually make the plots
"""
import pandas as pd

import plot_functions
import theory_functions
import numpy as np


def plot() -> None:
    """

    @return:
    """
    arguments = plot_functions.parse_command_line_arguments()
    length_range = arguments.length_range
    algorithm = arguments.algorithm

    if algorithm == "BOTH":
        # plot comparison plot
        # Fig.4.
        pass
    else:
        dimensionality_start = arguments.start_dimensionality
        dimensionality_end = arguments.end_dimensionality
        exponent_start = arguments.start_exponent
        exponent_end = arguments.end_exponent
        pdb_dataframe = plot_functions.get_dataframe(arguments)
        sim_dataframe = pd.read_csv(f"../data/simulations/3d/lls_{length_range}.csv")

        if dimensionality_start == dimensionality_end and exponent_start == exponent_end:
            pass
        elif dimensionality_start != dimensionality_end and exponent_start != exponent_end:
            plot_functions.create_grid_plots(arguments, pdb_dataframe, sim_dataframe)


def main():
    plot()


if __name__ == "__main__":
    main()
