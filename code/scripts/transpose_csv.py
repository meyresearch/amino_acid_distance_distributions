"""
Read in combined_ids.txt and save as a transposed column in a csv file
"""
import numpy as np

with open("../data/rcsb/combined_ids.txt", "r") as infile:
    data_list = infile.read().split(",")
    remove_last_bracket_list = [string.replace("]", "") for string in data_list]
    remove_first_bracket_list = [string.replace("[", "") for string in remove_last_bracket_list]
    clean_array = np.array([string.replace("\"", "") for string in remove_first_bracket_list])

transposed_array = np.transpose(clean_array)
np.savetxt("../data/rcsb/combined_ids_tr.csv", transposed_array, delimiter=",", fmt="%s")
