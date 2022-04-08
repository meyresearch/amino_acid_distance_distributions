"""
Handle all plotting functions and actually make the plots
"""
import pandas as pd
import distances
import plot_functions
import bars
import comparison
import adjacency
import simulation
import plot_comparison


def plot() -> None:
    """

    @return:
    """
    arguments = plot_functions.parse_command_line_arguments()
    length_range = arguments.length_range
    algorithm = arguments.algorithm

    if algorithm == "comp":
        rcsb_histogram = plot_comparison.get_histogram(arguments, "rcsb")
        alpha_histogram = plot_comparison.get_histogram(arguments, "alpha")
        # sim_histogram = plot_comparison.get_simulation_histogram(arguments.length_range)
        plot_comparison.create_comparison_plot(arguments, rcsb_histogram, alpha_histogram)
    elif algorithm == "bar":
        bars.create_bar_plots()
    elif algorithm == "rcsb" or algorithm == "alpha":
        dimensionality_start = arguments.start_dimensionality
        dimensionality_end = arguments.end_dimensionality
        exponent_start = arguments.start_exponent
        exponent_end = arguments.end_exponent
        pdb_histogram = plot_comparison.get_histogram(arguments, algorithm)
        sim_histogram = plot_comparison.get_simulation_histogram(arguments.length_range)
        if dimensionality_start == dimensionality_end and exponent_start == exponent_end:
            plot_comparisont.create_plots(arguments, exponent_start, dimensionality_start, pdb_histogram, sim_histogram)
        elif dimensionality_start != dimensionality_end and exponent_start != exponent_end:
            plot_comparison.create_grid_plots(arguments, pdb_histogram, sim_histogram)
    elif algorithm == "adj":
        adjacency.plot_adjacency_matrix(arguments.file, arguments.data_type)
    elif algorithm == "2d-sim":

        # if dimensionality_start == dimensionality_end and exponent_start == exponent_end:
        # simulation.create_2d_grid_plots(arguments)
        simulation.create_2d_plots(arguments)
        # else:
    elif algorithm == "3d-sim":
        simulation.create_3d_simulation_plot(arguments)

def main():
    plot()


if __name__ == "__main__":
    main()
