"""
Compute amino acid distances or bootstrap.
"""
import argparse
import functions


def handle_commandline_options(cl_arguments: argparse.Namespace) -> None:
    """
    Read in CL arguments and make a choice depending on given option
    @param cl_arguments from user
    @return: None
    """
    given_algorithm = cl_arguments.algorithm
    given_inputfile = cl_arguments.inputfile
    given_linklength = cl_arguments.length_range

    if given_algorithm == "PDB":
        functions.pdb_to_pcm(log_file="log.txt", given_algorithm=given_algorithm, length_range=given_linklength)
    elif given_algorithm == "aF":
        functions.pdb_to_pcm(log_file="log.txt", given_algorithm=given_algorithm, length_range=given_linklength)
    elif given_algorithm == "SIM":
        functions.compute_simulation_distribution(length_range=given_linklength)
    elif given_algorithm == "BS":
        functions.bootstrap(inputfile=given_inputfile, sample_replacement=True, length_range=given_linklength)
    elif given_algorithm == "C":
        functions.bootstrap(inputfile=given_inputfile, sample_replacement=False, length_range=given_linklength)


def main():
    arguments = functions.commandline_arguments()
    handle_commandline_options(arguments)


if __name__ == "__main__":
    main()
