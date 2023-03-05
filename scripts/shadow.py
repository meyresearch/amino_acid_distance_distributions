"""
Compute amino acid distances using Shadow map [https://doi.org/10.1021/jp300852d]
"""
import argparse
import functions as fn
import plot_functions as pf


def shadow_commandline_options() -> None:
    """
    Read in CL arguments and make a choice depending on given option
    @param cl_arguments from user
    @return: None
    """
    parser = argparse.ArgumentParser(description="Get amino acid distance distributions from Shadow maps.")
    parser.add_argument("-p", "--path", type=str, help="Full path to file contianing PDB IDs and paths")
    parser.add_argument("-d", "--distances", action="store_true", help="Compute distances from Shadow map")
    parser.add_argument("-pd", "--plot-distances", type=str, help="npy file for plotting Shadow map distances")
    return parser.parse_args()


def main():
    arguments = shadow_commandline_options()
    path = arguments.path
    distances = arguments.distances
    plot_file = arguments.plot_distances

    if distances:
        print("Computing distances from Shadow map")
        fn.get_shadow_distances(path)
    if plot_file:
        print("Plotting distances from Shadow map")
        histogram = pf.get_histogram(length_range="", algorithm=plot_file, low_confidence=False)
        


if __name__ == "__main__":
    main()
