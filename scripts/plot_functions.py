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
                           help="starting value for exponent (constant a)")
    pdb_group.add_argument("-stepe", "--step-exponent", dest="step_exponent", type=int, nargs="?", const=1, default=1,
                           help="step size for constant a range, default=1")
    pdb_group.add_argument("-stepd", "--step-dimensionality", dest="step_dimensionality", type=float, nargs="?",
                           const=0.0001, default=0.0001, help="step size for constant A range, default=0.0001")
    compare_group.add_argument("-rcsbd", "--rcsb-dimensionality", dest="d_rcsb", type=float,
                               help="constant A for rcsb algorithm in comparison plot")
    compare_group.add_argument("-alphad", "--alpha-dimensionality", dest="d_alpha", type=float,
                               help="constant A for alpha algorithm in comparison plot")
    compare_group.add_argument("-rcsbe", "--rcsbe-exponent", dest="e_rcsb", type=float,
                               help="constant a for rcsb algorithm in comparison plot")
    compare_group.add_argument("-alphae", "--alpha-exponent", dest="e_alpha", type=float,
                               help="constant a for alpha algorithm in comparison plot")
    arguments = parser.parse_args()
    check_required_arguments(parser, arguments.algorithm, arguments.length_range, arguments.start_dimensionality,
                             arguments.end_dimensionality, arguments.start_exponent, arguments.end_exponent,
                             arguments.d_rcsb, arguments.e_rcsb, arguments.d_alpha, arguments.e_alpha,
                             arguments.data_type, arguments.file)
    return arguments


def check_required_arguments(argument_parser: argparse.ArgumentParser, given_algorithm: str, given_length: str,
                             starting_d: str, ending_d: str, starting_e: str, ending_e: str, bootstrap_d: str,
                             bootstrap_e: str, chunk_d: str, chunk_e: str,
                             data_type: str, file: str) -> None:
    """
    Check given CL arguments so that all requirements are met
    @param file: PDB file or simulation matrix
    @param data_type: SIM or PDB
    @param given_length: given chain length range
    @param argument_parser: parser for CL arguments
    @param given_algorithm: either rcsb, alpha or comp
    @param starting_d: starting value for dimensionality (constant A)
    @param ending_d: ending value for dimensionality (constant A)
    @param starting_e: starting value for exponent (constant a)
    @param ending_e: ending value for exponent (constant a)
    @param bootstrap_d: dimensionality constant (A) for bootstrapping algorithm
    @param bootstrap_e: exponent constant (a) for bootstrapping algorithm
    @param chunk_d: dimensionality constant (A) for chunking algorithm
    @param chunk_e: exponent constant (a) for chunking algorithm
    @return: None
    """
    if (given_algorithm == "rcsb" and starting_d is None) or \
            (given_algorithm == "rcsb" and ending_d is None) or \
            (given_algorithm == "rcsb" and starting_e is None) or \
            (given_algorithm == "rcsb" and ending_e is None) or \
            (given_algorithm == "rcsb" and given_length is None) or \
            (given_algorithm == "rcsb" and bootstrap_d is not None) or \
            (given_algorithm == "rcsb" and bootstrap_e is not None) or \
            (given_algorithm == "rcsb" and chunk_e is not None) or \
            (given_algorithm == "rcsb" and chunk_d is not None) or \
            (given_algorithm == "rcsb" and data_type is not None) or \
            (given_algorithm == "rcsb" and file is not None):
        argument_parser.error("rcsb requires -r, -startd, -endd, -starte, -ende")

    elif (given_algorithm == "alpha" and starting_d is None) or \
            (given_algorithm == "alpha" and ending_d is None) or \
            (given_algorithm == "alpha" and starting_e is None) or \
            (given_algorithm == "alpha" and ending_e is None) or \
            (given_algorithm == "alpha" and given_length is None) or \
            (given_algorithm == "alpha" and bootstrap_d is not None) or \
            (given_algorithm == "alpha" and bootstrap_e is not None) or \
            (given_algorithm == "alpha" and chunk_e is not None) or \
            (given_algorithm == "alpha" and chunk_d is not None) or \
            (given_algorithm == "alpha" and data_type is not None) or \
            (given_algorithm == "alpha" and file is not None):
        argument_parser.error("rcsb requires -r, -startd, -endd, -starte, -ende")

    elif (given_algorithm == "comp" and bootstrap_d is None) or \
            (given_algorithm == "comp" and chunk_d is None) or \
            (given_algorithm == "comp" and bootstrap_e is None) or \
            (given_algorithm == "comp" and chunk_e is None) or \
            (given_algorithm == "comp" and given_length is None) or \
            (given_algorithm == "comp" and starting_d is not None) or \
            (given_algorithm == "comp" and ending_d is not None) or \
            (given_algorithm == "comp" and starting_e is not None) or \
            (given_algorithm == "comp" and ending_e is not None) or \
            (given_algorithm == "comp" and data_type is not None) or \
            (given_algorithm == "comp" and file is not None):
        argument_parser.error("comp requires -r, -rcsbd, -alphad, -rcsbe, -alphae")

    elif (given_algorithm == "bar" and given_length is not None) or \
            (given_algorithm == "bar" and starting_d is not None) or \
            (given_algorithm == "bar" and ending_d is not None) or \
            (given_algorithm == "bar" and starting_e is not None) or \
            (given_algorithm == "bar" and ending_e is not None) or \
            (given_algorithm == "bar" and bootstrap_d is not None) or \
            (given_algorithm == "bar" and bootstrap_e is not None) or \
            (given_algorithm == "bar" and chunk_e is not None) or \
            (given_algorithm == "bar" and chunk_d is not None) or \
            (given_algorithm == "bar" and data_type is not None) or \
            (given_algorithm == "bar" and file is not None):
        argument_parser.error("bar does not take any arguments")

    elif (given_algorithm == "adj" and data_type is None) or \
            (given_algorithm == "adj" and file is None) or \
            (given_algorithm == "adj" and given_length is not None) or \
            (given_algorithm == "adj" and starting_d is not None) or \
            (given_algorithm == "adj" and ending_d is not None) or \
            (given_algorithm == "adj" and starting_e is not None) or \
            (given_algorithm == "adj" and ending_e is not None) or \
            (given_algorithm == "adj" and bootstrap_d is not None) or \
            (given_algorithm == "adj" and bootstrap_e is not None) or \
            (given_algorithm == "adj" and chunk_e is not None) or \
            (given_algorithm == "adj" and chunk_d is not None):
        argument_parser.error("adj requires -t, -f")
