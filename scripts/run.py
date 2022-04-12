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
    given_path = cl_arguments.path_to_pdbs
    if given_algorithm == "rcsb":
        # functions.pdb_to_pcm(log_file="log.txt",
        #                      given_algorithm=given_algorithm,
        #                      length_range=given_linklength,
        #                      path_to_pdbs=given_path)
        functions.return_distance_histogram(log_file="log.txt",
                                            given_algorithm=given_algorithm,
                                            length_range=given_linklength,
                                            path_to_csvs=given_path)
    elif given_algorithm == "alpha":
        # functions.pdb_to_pcm(log_file="log.txt",
        #                      given_algorithm=given_algorithm,
        #                      length_range=given_linklength,
        #                      path_to_pdbs=given_path)
        functions.return_distance_histogram(log_file="log.txt",
                                            given_algorithm=given_algorithm,
                                            length_range=given_linklength,
                                            path_to_csvs=given_path)
    elif given_algorithm == "3d-sim":
        functions.return_3d_simulation_distance_histogram(length_range=given_linklength)
    elif given_algorithm == "2d-sim":
        functions.return_2d_simulation_distance_histogram()
    # elif given_algorithm == "boots":
    #     functions.bootstrap(inputfile=given_inputfile,
    #                         sample_replacement=True,
    #                         length_range=given_linklength)
    # elif given_algorithm == "chunk ":
    #     functions.bootstrap(inputfile=given_inputfile,
    #                         sample_replacement=False,
    #                         length_range=given_linklength)


def main():
    arguments = functions.commandline_arguments()
    handle_commandline_options(arguments)


if __name__ == "__main__":
    main()
