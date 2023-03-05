"""
Common functions used by all plotting codes, such as for handling command line arguments
"""
import argparse
import numpy as np
import theory_functions


def get_histogram(length_range: str, algorithm: str, low_confidence: bool) -> np.ndarray:
    """
    Get distribution of amino acid distances from histogram.
    @param length_range: chain length range
    @param algorithm: rcsb or alpha from command line arguments
    @param low_confidence: True use low confidence structures for alphafold 2, False use high confidence
    @return: histogram of RCSB/AlphaFold 2 data
    """
    histogram = None
    if algorithm == "rcsb":
        histogram = np.load(f"../data/rcsb/histogram_{length_range}_not_normed.npy", allow_pickle=True)
    elif algorithm == "alpha" and low_confidence:
        print("low")
        histogram = np.load(f"../data/alphafold/histogram_low_conf_{length_range}_not_normed.npy", allow_pickle=True)
    elif algorithm == "alpha" and not low_confidence:
        histogram = np.load(f"../data/alphafold/histogram_{length_range}_not_normed.npy", allow_pickle=True)
    else: 
        histogram = np.load(algorithm, allow_pickle=True)
    return histogram


def get_data_for_plotting(histogram: np.ndarray, arguments: argparse.Namespace, length_range: str) -> tuple:
    """
    Get protein data for plotting from dataframe
    @param histogram: numpy array that contains amino acid distance_bins
    @param arguments: command line arguments
    @param length_range: chain length range
    @return: tuple of chain length, distances, normalised measure of central tendency and confidence level bounds
    """
    distance_bins = np.linspace(start=1, stop=350, num=350)[:-1]
    lower_bound, upper_bound = theory_functions.get_confidence_interval(histogram, arguments.quantile)
    measure = theory_functions.get_measure_of_central_tendency(histogram, arguments.measure)
    normalised_measure = measure / np.sum(measure)
    normalised_lower_bound = lower_bound / np.sum(measure)
    normalised_upper_bound = upper_bound / np.sum(measure)
    chain_length = int(length_range)
    return chain_length, distance_bins, normalised_measure, normalised_lower_bound, normalised_upper_bound


def parse_command_line_arguments() -> argparse.Namespace:
    """
    Parse CL arguments for plotting.
    @return: Namespace containing CL arguments
    """
    parser = argparse.ArgumentParser(description="Plot amino acid distances and residuals.")
    sim_group = parser.add_argument_group("2d-sim, 3d-sim")
    compare_group = parser.add_argument_group("comp")
    adjacency_group = parser.add_argument_group("adj")
    parser.add_argument("-r", "--range", dest="length_range", type=str, choices=[None, "100", "200", "300"],
                        help="chain length range to be plotted; required for comp, 2d-sim, 3d-sim")
    parser.add_argument("algorithm", type=str, choices=["rcsb",
                                                        "alpha",
                                                        "comp",
                                                        "2d-sim",
                                                        "3d-sim",
                                                        "adj",
                                                        "bar",
                                                        "cont"],
                        help="See github readme for usages.")
    adjacency_group.add_argument("-t", "--type", dest="data_type", type=str, choices=["pdb", "sim"],
                                 help="data type for adjacency matrix")
    adjacency_group.add_argument("-f", "--file", dest="file", type=str,
                                 help="PDB file or simulation matrix for adjacency matrix")
    sim_group.add_argument("-startd", "--start-dimensionality",
                           dest="start_dimensionality", type=float,
                           help="starting value for dimensionality constant (A)")
    sim_group.add_argument("-endd", "--end-dimensionality", dest="end_dimensionality", type=float,
                           help="last value for dimensionality constant (A)")
    sim_group.add_argument("-starte", "--start-exponent", dest="start_exponent", type=int,
                           help="starting value for exponent (constant a)")
    sim_group.add_argument("-ende", "--end-exponent", dest="end_exponent", type=int,
                           help="last value for exponent (constant a)")
    parser.add_argument("-stepe", "--step-exponent", dest="step_exponent", type=int,
                        help="step size for constant a range")
    parser.add_argument("-stepd", "--step-dimensionality", dest="step_dimensionality", type=float, nargs="?",
                        const=0.0001, default=0.0001, help="step size for constant A range, default=0.0001")
    parser.add_argument("-m", "--measure", dest="measure", type=str, choices=["median", "mean"], nargs="?",
                        const="mean", default="mean", help="measure of central tendency to use")
    parser.add_argument("-q", "--quantile", dest="quantile", type=int, choices=[1, 2, 3], nargs="?", const=1,
                        default=1, help="confidence interval to use in multiples of standard deviations (std), "
                                        "i.e. 1 std, 2 std or 3 std")
    parser.add_argument("-startp", "--start-point", dest="start_point", type=int, nargs="?", const=1, default=1,
                        help="point from which to start fitting")
    parser.add_argument("-endp", "--end-point", dest="end_point", type=int, help="point at which to end fitting")
    parser.add_argument("-l", "--low",
                    help="get low-confidence AlphaFold 2 structures, default behaviour is to only get structures with confidence above 90",
                    action="store_true")

    compare_group.add_argument("--rcsb-startd", dest="rcsb_startd", type=float,
                               help="starting value for dimensionality constant (A) for rcsb")
    compare_group.add_argument("--rcsb-endd", dest="rcsb_endd", type=float,
                               help="last value for dimensionality constant (A) for rcsb")
    compare_group.add_argument("--rcsb-starte", dest="rcsb_starte", type=float,
                               help="starting value for exponent constant (a) for rcsb")
    compare_group.add_argument("--rcsb-ende", dest="rcsb_ende", type=float,
                               help="last value for exponent constant (a) for rcsb")

    

    arguments = parser.parse_args()
    if arguments.algorithm != "rcsb" or arguments.algorithm != "alpha" or arguments.algorithm != "bar":
        check_required_arguments(parser, arguments.algorithm, arguments.length_range, arguments.start_dimensionality,
                                 arguments.end_dimensionality, arguments.start_exponent, arguments.end_exponent,
                                 arguments.rcsb_startd, arguments.rcsb_endd, arguments.rcsb_starte, arguments.rcsb_ende,
                                 arguments.data_type, arguments.file,
                                 arguments.start_point, arguments.end_point)
    return arguments


def check_required_arguments(argument_parser: argparse.ArgumentParser, given_algorithm: str, given_length: str,
                             start_dimensionality: str, end_dimensionality: str,
                             start_exponent: str, end_exponent: str,
                             rcsb_startd: str, rcsb_endd: str, rcsb_starte: str, rcsb_ende: str,
                             data_type: str, file: str,
                             start_point: int, end_point: int) -> None:
    """
    Check given CL arguments so that all requirements are met
    @param file: PDB file or simulation matrix
    @param data_type: SIM or PDB
    @param given_length: given chain length range
    @param argument_parser: parser for CL arguments
    @param given_algorithm: either rcsb, alpha or comp
    @param start_dimensionality: starting value for dimensionality (constant A)
    @param end_dimensionality: ending value for dimensionality (constant A)
    @param start_exponent: starting value for exponent (constant a)
    @param end_exponent: ending value for exponent (constant a)
    @param rcsb_starte: starting value for exponent (constant a) for rcsb
    @param rcsb_ende: last value for exponent constant (a) for rcsb
    @param rcsb_startd: starting value for dimensionality constant (A) for rcsb
    @param rcsb_endd: last value for exponent constant (a) for rcsb
    @param start_point: point from which to start plotting
    @param end_point: point at which to end plotting
    @return: None
    """
    if (given_algorithm == "2d-sim" and start_dimensionality is None) or \
            (given_algorithm == "2d-sim" and end_dimensionality is None) or \
            (given_algorithm == "2d-sim" and start_exponent is None) or \
            (given_algorithm == "2d-sim" and end_exponent is None) or \
            (given_algorithm == "2d-sim" and given_length is not None) or \
            (given_algorithm == "2d-sim" and data_type is not None) or \
            (given_algorithm == "2d-sim" and file is not None) or \
            (given_algorithm == "2d-sim" and start_point is None) or \
            (given_algorithm == "2d-sim" and end_point is None):
        argument_parser.error("2d-sim requires -r=None -startd, -endd, -starte, -ende, -startp, -endp")

    elif (given_algorithm == "3d-sim" and start_dimensionality is None) or \
            (given_algorithm == "3d-sim" and end_dimensionality is None) or \
            (given_algorithm == "3d-sim" and start_exponent is None) or \
            (given_algorithm == "3d-sim" and end_exponent is None) or \
            (given_algorithm == "3d-sim" and given_length is None) or \
            (given_algorithm == "3d-sim" and rcsb_startd is not None) or \
            (given_algorithm == "3d-sim" and rcsb_endd is not None) or \
            (given_algorithm == "3d-sim" and rcsb_starte is not None) or \
            (given_algorithm == "3d-sim" and rcsb_ende is not None) or \
            (given_algorithm == "3d-sim" and data_type is not None) or \
            (given_algorithm == "3d-sim" and file is not None) or \
            (given_algorithm == "3d-sim" and start_point is None) or \
            (given_algorithm == "3d-sim" and end_point is None):
        argument_parser.error("3d-sim requires -r, -startd, -endd, -starte, -ende, -startp, -endp")

    elif (given_algorithm == "comp" and rcsb_startd is None) or \
            (given_algorithm == "comp" and rcsb_endd is None) or \
            (given_algorithm == "comp" and rcsb_starte is None) or \
            (given_algorithm == "comp" and rcsb_ende is None) or \
            (given_algorithm == "comp" and given_length is None) or \
            (given_algorithm == "comp" and start_dimensionality is not None) or \
            (given_algorithm == "comp" and end_dimensionality is not None) or \
            (given_algorithm == "comp" and start_exponent is not None) or \
            (given_algorithm == "comp" and end_exponent is not None) or \
            (given_algorithm == "comp" and data_type is not None) or \
            (given_algorithm == "comp" and file is not None) or \
            (given_algorithm == "comp" and end_point is None):
        argument_parser.error("comp requires -r, --rcsb-startd, --rcsb-endd, --rcsb-starte, "
                              "--rcsb-ende, -startp, -endp")

    elif (given_algorithm == "adj" and data_type is None) or \
            (given_algorithm == "adj" and file is None) or \
            (given_algorithm == "adj" and given_length is not None) or \
            (given_algorithm == "adj" and start_dimensionality is not None) or \
            (given_algorithm == "adj" and end_dimensionality is not None) or \
            (given_algorithm == "adj" and start_exponent is not None) or \
            (given_algorithm == "adj" and end_exponent is not None) or \
            (given_algorithm == "adj" and rcsb_startd is not None) or \
            (given_algorithm == "adj" and rcsb_endd is not None) or \
            (given_algorithm == "adj" and rcsb_starte is not None) or \
            (given_algorithm == "adj" and rcsb_ende is not None):
        argument_parser.error("adj requires -t, -f")
