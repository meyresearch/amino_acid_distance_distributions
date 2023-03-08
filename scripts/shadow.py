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
    parser.add_argument("-p", "--path", type=str, help="Full path to input file")
    parser.add_argument("-c", "--cutoff", type=str, help="Cutoff distance in Angstrom")
    parser.add_argument("-s", "--shadow-radius", type=str, help="Shadowing radius")
    parser.add_argument("-d", "--distances", action="store_true", help="Compute distances from Shadow map")
    parser.add_argument("-pd", "--plot-distances", action="store_true", help="Plot Shadow map distances")
    parser.add_argument("-sg", "--smog", action="store_true", help="Compute Shadow maps from SMOG2")
    return parser.parse_args()


def main():
    arguments = shadow_commandline_options()
    path = arguments.path
    cutoff_radius = arguments.cutoff
    shadow_radius = arguments.shadow_radius
    distances = arguments.distances
    plot = arguments.plot_distances
    smog = arguments.smog
    if smog:
        print("Computing Shadow maps from SMOG2")
        fn.run_smog(path, cutoff_radius, shadow_radius)
    if distances:
        print("Computing distances from Shadow map")
        fn.get_shadow_distance_histograms(path, cutoff_radius, shadow_radius)
    if plot:
        print("Plotting distances from Shadow map")
        histogram = pf.get_histogram(length_range="", algorithm=path, low_confidence=False)
        


if __name__ == "__main__":
    main()
