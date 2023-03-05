"""
Get amino acid distance distributions for a user-specified value of C
"""
import argparse
import functions


def commandline_arguments():
    parser = argparse.ArgumentParser(description="Get amino acid distance distributions for a user-specified value of C")
    parser.add_argument("-a", "--algorithm", type=str, choices=["rcsb", "alpha"], default="alpha")
    parser.add_argument("-c", "--cutoff", type=str, help="Cutoff value")
    parser.add_argument("-r", "--range", type=str, choices=["100", "200", "300"], default="300", help="Chain length range")
    parser.add_argument("-p", "--path", type=str, help="Path to csv file containing PDB data")
    return parser.parse_args()
    

def main():
    arguments = commandline_arguments()
    algorithm = arguments.algorithm
    cutoff = arguments.cutoff
    length_range = arguments.range
    path = arguments.path
    functions.get_distances_with_different_cutoff(algorithm, length_range, path, cutoff)


if __name__ == "__main__":
    main()