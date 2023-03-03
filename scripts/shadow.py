"""
Compute amino acid distances using Shadow map [https://doi.org/10.1021/jp300852d]
"""
import argparse
import functions as fn


# def handle_commandline_options(cl_arguments: argparse.Namespace) -> None:
#     """
#     Read in CL arguments and make a choice depending on given option
#     @param cl_arguments from user
#     @return: None
#     """
#     given_algorithm = cl_arguments.algorithm
#     given_path = cl_arguments.path_to_pdbs
#     if given_algorithm == "rcsb":
#         pass
#     elif given_algorithm == "alpha":


def main():
    arguments = fn.commandline_arguments()
    algorithm = arguments.algorithm
    length_range = arguments.length_range
    path = arguments.path_to_pdbs

    print(f"{algorithm} {length_range} {path}")


if __name__ == "__main__":
    main()
