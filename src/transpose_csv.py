import pandas as pd
import numpy as np

csv_table = np.genfromtxt("../data/combined_ids.csv",delimiter=',',dtype=str)
print(csv_table[-1])
transposed = csv_table.T
np.savetxt("../data/combined_ids_tr.csv", transposed, fmt="%s")