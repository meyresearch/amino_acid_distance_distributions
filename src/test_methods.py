import random
import numpy as np
import pandas as pd

sample_size = 1000

sample = np.arange(0,21,0.1)

DF = pd.read_csv('../data/ids/id_file_0.txt', header = None)
DF = DF.transpose()
print(DF)
