import glob
import pandas as pd
import numpy as np


def get_confidences(pdb_files: list) -> list:
    """
    Read AlphaFold PDB file contents to find confidence info
    @param: pdbs glob list of AlphaFold PDB files
    @return: confidence_info list containing the dictionary of confidence info
    """
    counter = 1
    confidence_info = []
    for pdb_file in pdb_files:
        print(f"At {counter}/{len(pdb_files)}")
        file = open(pdb_file)
        lines = file.readlines()
        file.close()
        confidence_scores = []
        confidence_dict = {}
        for line in lines:
            if line.startswith("ATOM"):
                current_line = line.split()
                confidence_score = float(current_line[-2])
                confidence_scores.append(confidence_score)
        confidence_dict["mean"] = np.mean(confidence_scores)
        confidence_dict["std"] = np.std(confidence_scores)
        confidence_dict["filename"] = pdb_file
        confidence_info.append(confidence_dict)
        counter += 1
    return confidence_info
                                          
                                      
def filter_confidences(pdb_files: list) -> pd.DataFrame:
    """
    Filter out AlphaFold structures with less than 90% confidence
    @param pdb_files: glob list of AlphaFold PDBs
    @return good_confidence_df: dataframe containing PDBs with 90% confidence or above
    """
    confidence_list = get_confidences(pdb_files)
    confidence_dataframe = pd.DataFrame.from_dict(confidence_list)
    confidence_filter = confidence_dataframe["mean"] > 90.0
    good_confidence_df = confidence_dataframe.where(confidence_filter).dropna()
    return good_confidence_df


pdbs = glob.glob("../data/alphafold/pdb_files/*.pdb")
confidence_df = dssp_functions.filter_confidences(pdbs)
confidence_df.to_csv("../data/alphafold/confidences.csv")
