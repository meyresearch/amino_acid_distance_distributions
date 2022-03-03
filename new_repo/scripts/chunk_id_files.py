"""
Read in combined_ids_tr.csv and save files containing chunks of 1000 IDs.
"""

import pandas as pd
import numpy as np

infile = "../data/combined_ids_tr.csv"
chunk = 1000
filename_counter = 0

for chunk in pd.read_csv(infile, chunksize=chunk, sep="\n", header=None):
    chunk_array = chunk.to_numpy()
    np.savetxt(f"../data/rcsb/ids/id_file_{filename_counter}.txt", chunk_array, fmt="%s")
    filename_counter += 1
