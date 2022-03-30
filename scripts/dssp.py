# Load data from AlphaFold for chain-lengths 100, 200, 300

import dssp_functions
import glob
import numpy as np
import pandas as pd
import argparse


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(description="Get secondary structure information from AlphaFold 2 or RCSB PDBs.")
    parser.add_argument("algorithm", type=str, choices=["alpha", "rcsb"])
    return parser.parse_args()


algorithm = parse_arguments().algorithm

if algorithm == "alpha":
    dssp_functions.save_alpha_secondary_info()

elif algorithm == "rcsb":
    dssp_functions.save_rcsb_secondary_info()


                
            
