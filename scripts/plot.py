"""
Handle all plotting functions and actually make the plots
"""
import adjacency
import bars
import comparison
import contour
import distances
import plot_functions
import simulation


def plot() -> None:
    """
    Handle command line arguments and call required plotting function
    @return: None
    """
    arguments = plot_functions.parse_command_line_arguments()
    algorithm = arguments.algorithm
    if algorithm == "comp":
        rcsb_histogram = plot_functions.get_histogram(arguments.length_range, "rcsb", arguments.low)
        alpha_histogram = plot_functions.get_histogram(arguments.length_range, "alpha", arguments.low)
        comparison.create_comparison_plot(arguments, rcsb_histogram, alpha_histogram)
    elif algorithm == "bar":
        bars.create_bar_plots()
    elif algorithm == "rcsb" or algorithm == "alpha":
        distances.create_plots(arguments)
    elif algorithm == "adj":
        adjacency.plot_adjacency_matrix(arguments.file, arguments.data_type)
    elif algorithm == "2d-sim":
        simulation.create_2d_plots(arguments)
    elif algorithm == "3d-sim":
        simulation.create_3d_simulation_plot(arguments)
    elif algorithm == "cont":
        contour.create_contour_plots()


def main():
    plot()


if __name__ == "__main__":
    main()
