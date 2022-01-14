import pandas as pd
import numpy as np

filename = '../data/combined_ids_tr.csv'
chunksize = 1000
filename_counter = 0

for chunk in pd.read_csv(filename, chunksize=chunksize, sep='\n', header=None):

    np_array = chunk.to_numpy()
    np.savetxt('../data/ids/id_file_'+str(filename_counter)+'.txt', np_array, fmt='%s')

    filename_counter += 1

