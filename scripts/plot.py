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


def plot() -> None:
    """

    @return:
    """
    arguments = plot_functions.parse_command_line_arguments()
    length_range = arguments.length_range
    algorithm = arguments.algorithm

    if algorithm == "comp":
        comparison.create_comparison_plot(arguments)
    elif algorithm == "bar":
        bars.create_bar_plots()
    elif algorithm == "rcsb" or algorithm == "alpha":
        dimensionality_start = arguments.start_dimensionality
        dimensionality_end = arguments.end_dimensionality
        exponent_start = arguments.start_exponent
        exponent_end = arguments.end_exponent
        pdb_histogram = distances.get_histogram(arguments)
        sim_dataframe = pd.read_csv(f"../data/simulations/3d/simulation_{length_range}_stats.csv")
        if dimensionality_start == dimensionality_end and exponent_start == exponent_end:
            distances.create_plots(arguments, exponent_start, dimensionality_start, pdb_histogram, sim_dataframe)
        elif dimensionality_start != dimensionality_end and exponent_start != exponent_end:
            distances.create_grid_plots(arguments, pdb_histogram, sim_dataframe)
    elif algorithm == "adj":
        adjacency.plot_adjacency_matrix(arguments.file, arguments.data_type)
    elif algorithm == "2d-sim":
        dimensionality_start = arguments.start_dimensionality
        dimensionality_end = arguments.end_dimensionality
        exponent_start = arguments.start_exponent
        exponent_end = arguments.end_exponent
        if dimensionality_start == dimensionality_end and exponent_start == exponent_end:
            simulation.create_2d_plots(exponent=exponent_start, dimensionality=dimensionality_start)
        else:
            simulation.create_2d_grid_plots(arguments)


def main():
    plot()


if __name__ == "__main__":
    main()
