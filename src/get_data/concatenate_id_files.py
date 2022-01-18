"""
concat_id_files
"""
import os
import pandas as pd
import numpy as np


def concat_id_files(inpath: str) -> np.ndarray:
    """
    @param inpath: input path for id files, default "../../data/ids/"
    @return: concatenated array of pdb ids
    """
    temporary_dataframes = []
    for filename in os.listdir(inpath):
        if filename.endswith(".txt"):
            temporary_dataframe = pd.read_csv(inpath + str(filename), header=None)
            temporary_dataframes.append(temporary_dataframe)

    concatenated_dataframe = pd.concat(temporary_dataframes).reset_index()
    concatenated_dataframe_no_index = concatenated_dataframe.drop(columns="index")
    pdb_id_array = concatenated_dataframe_no_index[0].to_numpy()

    return pdb_id_array

