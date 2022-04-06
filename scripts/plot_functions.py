"""
Functions used for plotting amino acid distance distributions.
"""
import argparse


def parse_command_line_arguments() -> argparse.Namespace:
    """
    Parse CL arguments for plotting.
    @return: Namespace containing CL arguments
    """
    parser = argparse.ArgumentParser(description="Plot amino acid distances and residuals.")
    pdb_group = parser.add_argument_group("rcsb, alpha, 2d")
    compare_group = parser.add_argument_group("comp")
    adjacency_group = parser.add_argument_group("adj")
    parser.add_argument("-r", "--range", dest="length_range", type=str, choices=[None, "100", "200", "300"],
                        help="chain length range to be plotted; required for rcsb, alpha, comp, 2d-sim")

    parser.add_argument("algorithm", type=str, choices=["rcsb",
                                                        "alpha",
                                                        "comp",
                                                        "2d-sim",
                                                        "adj",
                                                        "bar"],
                        help="plot distances from rcsb, [alpha]Fold, do a [comp]arison, distances from 2D "
                             "simulations (2d-sim), plot [adj]acency matrix or plot bar plots")
    adjacency_group.add_argument("-t", "--type", dest="data_type", type=str, choices=["pdb", "sim"],
                                 help="data type for adjacency matrix")
    adjacency_group.add_argument("-f", "--file", dest="file", type=str,
                                 help="PDB file or simulation matrix for adjacency matrix")
    pdb_group.add_argument("-startd", "--start-dimensionality",
                           dest="start_dimensionality", type=float,
                           help="starting value for dimensionality constant (A)")
    pdb_group.add_argument("-endd", "--end-dimensionality", dest="end_dimensionality", type=float,
                           help="last value for dimensionality constant (A)")
    pdb_group.add_argument("-starte", "--start-exponent", dest="start_exponent", type=int,
                           help="starting value for exponent (constant a)")
    pdb_group.add_argument("-ende", "--end-exponent", dest="end_exponent", type=int,
                           help="last value for exponent (constant a)")
    parser.add_argument("-stepe", "--step-exponent", dest="step_exponent", type=int, nargs="?", const=1, default=1,
                        help="step size for constant a range, default=1")
    parser.add_argument("-stepd", "--step-dimensionality", dest="step_dimensionality", type=float, nargs="?",
                        const=0.0001, default=0.0001, help="step size for constant A range, default=0.0001")
    parser.add_argument("-m", "--measure", dest="measure", type=str, choices=["median", "mean"], nargs="?",
                        const="median", default="median", help="measure of central tendency to use")
    parser.add_argument("-q", "--quantile", dest="quantile", type=int, choices=[1, 2, 3], nargs="?", const=1,
                        default=1, help="confidence interval to use in multiples of standard deviations (std), "
                                        "i.e. 1 std, 2 std or 3 std")
    parser.add_argument("-startp", "--start-point", dest="start_point", type=int, nargs="?", const=1, default=1,
                        help="point from which to start plotting")
    parser.add_argument("-endp", "--end-point", dest="end_point", type=int, help="point at which to end plotting")
    compare_group.add_argument("--rcsb-startd", dest="rcsb_startd", type=float,
                               help="starting value for dimensionality constant (A) for rcsb")
    compare_group.add_argument("--rcsb-endd", dest="rcsb_endd", type=float,
                               help="last value for dimensionality constant (A) for rcsb")
    compare_group.add_argument("--alpha-startd", dest="alpha_startd", type=float,
                               help="starting value for dimensionality constant (A) for alpha")
    compare_group.add_argument("--alpha-endd", dest="alpha_endd", type=float,
                               help="last value for dimensionality constant (A) for alpha")
    compare_group.add_argument("--rcsb-starte", dest="rcsb_starte", type=float,
                               help="starting value for exponent constant (a) for rcsb")
    compare_group.add_argument("--rcsb-ende", dest="rcsb_ende", type=float,
                               help="last value for exponent constant (a) for rcsb")
    compare_group.add_argument("--alpha-starte", dest="alpha_starte", type=float,
                               help="starting value for exponent constant (a) for alpha")
    compare_group.add_argument("--alpha-ende", dest="alpha_ende", type=float,
                               help="last value for exponent constant (a) for alpha")

    arguments = parser.parse_args()
    check_required_arguments(parser, arguments.algorithm, arguments.length_range, arguments.start_dimensionality,
                             arguments.end_dimensionality, arguments.start_exponent, arguments.end_exponent,
                             arguments.rcsb_startd, arguments.rcsb_endd, arguments.rcsb_starte, arguments.rcsb_ende,
                             arguments.alpha_startd, arguments.alpha_endd, arguments.alpha_starte, arguments.alpha_ende,
                             arguments.data_type, arguments.file, arguments.end_point)
    return arguments


def check_required_arguments(argument_parser: argparse.ArgumentParser, given_algorithm: str, given_length: str,
                             start_dimensionality: str, end_dimensionality: str,
                             start_exponent: str, end_exponent: str,
                             rcsb_startd: str, rcsb_endd: str, rcsb_starte: str, rcsb_ende: str,
                             alpha_startd: str, alpha_endd: str, alpha_starte: str, alpha_ende: str,
                             data_type: str, file: str, end_point: int) -> None:
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
    @param alpha_starte: starting value for exponent (constant a) for alpha
    @param alpha_ende: last value for exponent constant (a) for alpha
    @param alpha_startd: starting value for dimensionality constant (A) for alpha
    @param alpha_endd: last value for exponent constant (a) for alpha
    @param end_point: point at which to end plotting
    @return: None
    """
    if (given_algorithm == "rcsb" and start_dimensionality is None) or \
            (given_algorithm == "rcsb" and end_dimensionality is None) or \
            (given_algorithm == "rcsb" and start_exponent is None) or \
            (given_algorithm == "rcsb" and end_exponent is None) or \
            (given_algorithm == "rcsb" and given_length is None) or \
            (given_algorithm == "rcsb" and rcsb_startd is not None) or \
            (given_algorithm == "rcsb" and rcsb_endd is not None) or \
            (given_algorithm == "rcsb" and rcsb_starte is not None) or \
            (given_algorithm == "rcsb" and rcsb_ende is not None) or \
            (given_algorithm == "rcsb" and alpha_startd is not None) or \
            (given_algorithm == "rcsb" and alpha_endd is not None) or \
            (given_algorithm == "rcsb" and alpha_starte is not None) or \
            (given_algorithm == "rcsb" and alpha_ende is not None) or \
            (given_algorithm == "rcsb" and data_type is not None) or \
            (given_algorithm == "rcsb" and file is not None) or \
            (given_algorithm == "rcsb" and end_point is None):
        argument_parser.error("rcsb requires -r, -startd, -endd, -starte, -ende, -endp")

    elif (given_algorithm == "alpha" and start_dimensionality is None) or \
            (given_algorithm == "alpha" and end_dimensionality is None) or \
            (given_algorithm == "alpha" and start_exponent is None) or \
            (given_algorithm == "alpha" and end_exponent is None) or \
            (given_algorithm == "alpha" and given_length is None) or \
            (given_algorithm == "alpha" and rcsb_startd is not None) or \
            (given_algorithm == "alpha" and rcsb_endd is not None) or \
            (given_algorithm == "alpha" and rcsb_starte is not None) or \
            (given_algorithm == "alpha" and rcsb_ende is not None) or \
            (given_algorithm == "alpha" and alpha_startd is not None) or \
            (given_algorithm == "alpha" and alpha_endd is not None) or \
            (given_algorithm == "alpha" and alpha_starte is not None) or \
            (given_algorithm == "alpha" and alpha_ende is not None) or \
            (given_algorithm == "alpha" and data_type is not None) or \
            (given_algorithm == "alpha" and file is not None) or \
            (given_algorithm == "alpha" and end_point is None):
        argument_parser.error("rcsb requires -r, -startd, -endd, -starte, -ende, -endp")

    elif (given_algorithm == "comp" and rcsb_startd is None) or \
            (given_algorithm == "comp" and rcsb_endd is None) or \
            (given_algorithm == "comp" and rcsb_starte is None) or \
            (given_algorithm == "comp" and rcsb_ende is None) or \
            (given_algorithm == "comp" and alpha_startd is None) or \
            (given_algorithm == "comp" and alpha_endd is None) or \
            (given_algorithm == "comp" and alpha_starte is None) or \
            (given_algorithm == "comp" and alpha_ende is None) or \
            (given_algorithm == "comp" and given_length is None) or \
            (given_algorithm == "comp" and start_dimensionality is not None) or \
            (given_algorithm == "comp" and end_dimensionality is not None) or \
            (given_algorithm == "comp" and start_exponent is not None) or \
            (given_algorithm == "comp" and end_exponent is not None) or \
            (given_algorithm == "comp" and data_type is not None) or \
            (given_algorithm == "comp" and file is not None) or \
            (given_algorithm == "comp" and end_point is None):
        argument_parser.error("comp requires -r, -rcsbd, -alphad, -rcsbe, -alphae, -endp")

    elif (given_algorithm == "bar" and given_length is not None) or \
            (given_algorithm == "bar" and start_dimensionality is not None) or \
            (given_algorithm == "bar" and end_dimensionality is not None) or \
            (given_algorithm == "bar" and start_exponent is not None) or \
            (given_algorithm == "bar" and end_exponent is not None) or \
            (given_algorithm == "bar" and rcsb_startd is not None) or \
            (given_algorithm == "bar" and rcsb_endd is not None) or \
            (given_algorithm == "bar" and rcsb_starte is not None) or \
            (given_algorithm == "bar" and rcsb_ende is not None) or \
            (given_algorithm == "bar" and alpha_startd is not None) or \
            (given_algorithm == "bar" and alpha_endd is not None) or \
            (given_algorithm == "bar" and alpha_starte is not None) or \
            (given_algorithm == "bar" and alpha_ende is not None) or \
            (given_algorithm == "bar" and data_type is not None) or \
            (given_algorithm == "bar" and file is not None):
        argument_parser.error("bar does not take any arguments")

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
            (given_algorithm == "adj" and rcsb_ende is not None) or \
            (given_algorithm == "adj" and alpha_startd is not None) or \
            (given_algorithm == "adj" and alpha_endd is not None) or \
            (given_algorithm == "adj" and alpha_starte is not None) or \
            (given_algorithm == "adj" and alpha_ende is not None):
        argument_parser.error("adj requires -t, -f")
