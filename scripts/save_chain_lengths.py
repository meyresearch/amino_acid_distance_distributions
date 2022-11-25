# save chain lengths
import subprocess
import numpy as np
import os
import protein_contact_map
import glob 
import pandas as pd
import argparse


def save_alpha_chain_lengths(confidence_dataframe: pd.DataFrame, low_confidence: bool) -> None:
    """
    Open the confidence dataframe and sort files into length ranges
    @param confidence_dataframe: dataframe containing confidence info
    @param low_confidence: True use low confidence structures, False use high confidence structures
    return: None
    """
    counter = 1
    pdb_files = confidence_dataframe["filename"].to_numpy()
    mean_confidences = confidence_dataframe["mean"].to_numpy()
    std_confidences = confidence_dataframe["std"].to_numpy()
    proteins_100 = {"filename": [], "mean_conf": [], "std_conf": []}
    proteins_200 = {"filename": [], "mean_conf": [], "std_conf": []}
    proteins_300 = {"filename": [], "mean_conf": [], "std_conf": []}
    for i in range(len(pdb_files)):
        print(f"At {counter}/{len(pdb_files)}")
        pcm = protein_contact_map.ProteinContactMap(pdb_files[i])
        alpha_carbons = pcm.get_alpha_carbons
        chain_length = protein_contact_map.get_chain_length(alpha_carbons)
        if chain_length in range(85, 116):
            print("In 100s")
            proteins_100["filename"].append(pdb_files[i])
            proteins_100["mean_conf"].append(mean_confidences[i])
            proteins_100["std_conf"].append(std_confidences[i])
        elif chain_length in range(185, 216):
            print("In 200s")
            proteins_200["filename"].append(pdb_files[i])
            proteins_200["mean_conf"].append(mean_confidences[i])
            proteins_200["std_conf"].append(std_confidences[i])
        elif chain_length in range(285, 316):
            print("In 300s")
            proteins_300["filename"].append(pdb_files[i])
            proteins_300["mean_conf"].append(mean_confidences[i])
            proteins_300["std_conf"].append(std_confidences[i])
        counter += 1
    df_100 = pd.DataFrame.from_dict(proteins_100)
    df_200 = pd.DataFrame.from_dict(proteins_200)
    df_300 = pd.DataFrame.from_dict(proteins_300)

    if low_confidence: 
        df_100.to_csv("../data/alphafold/low_confidences_100.csv")
        df_200.to_csv("../data/alphafold/low_confidences_200.csv")
        df_300.to_csv("../data/alphafold/low_confidences_300.csv")
    elif not low_confidence:
        df_100.to_csv("../data/alphafold/confidences_100.csv")
        df_200.to_csv("../data/alphafold/confidences_200.csv")
        df_300.to_csv("../data/alphafold/confidences_300.csv")
    
    
def save_rcsb_chain_lengths(pdb_files: list) -> None:
    """
    Sort RCSB PDBs into length ranges.
    @param: pdbs: glob list of PDB files
    @return: None
    """
    proteins_100 = {"filename": []}
    proteins_200 = {"filename": []}
    proteins_300 = {"filename": []}
    counter = 1
    for pdb_file in pdb_files:
        print(f"Process: {counter}/{len(pdb_files)}")
        pcm = protein_contact_map.ProteinContactMap(pdb_file)
        alpha_carbons = pcm.get_alpha_carbons
        chain_length = protein_contact_map.get_chain_length(alpha_carbons)
        
        if chain_length in range(85, 116):
            print("In 100s")
            proteins_100["filename"].append(pdb_file)
        elif chain_length in range(185, 216):
            print("In 200s")
            proteins_200["filename"].append(pdb_file)
        elif chain_length in range(285, 316):
            print("In 300s")
            proteins_300["filename"].append(pdb_file)
        counter += 1
    df_100 = pd.DataFrame.from_dict(proteins_100)
    df_200 = pd.DataFrame.from_dict(proteins_200)
    df_300 = pd.DataFrame.from_dict(proteins_300)
    
    df_100.to_csv("../data/rcsb/chains_100.csv")
    df_200.to_csv("../data/rcsb/chains_200.csv")
    df_300.to_csv("../data/rcsb/chains_300.csv")
    
    
def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(description="Sort PDBs to chain length ranges.")
    parser.add_argument("algorithm", type=str, choices=["alpha", "rcsb"])
    parser.add_argument("-l",
                        "--low",
                        help="get low-confidence AlphaFold 2 structures, default behaviour is to only get structures with confidence above 90",
                        action="store_true")
    return parser.parse_args()

arguments = parse_arguments()
algorithm = arguments.algorithm
low_confidence = arguments.low

if algorithm == "alpha":
    print("Getting AlphaFold length ranges.")
    if low_confidence:
        confidence_df = pd.read_csv("../data/alphafold/low_confidences.csv")
        save_alpha_chain_lengths(confidence_df, low_confidence)
    elif not low_confidence:
        confidence_df = pd.read_csv("../data/alphafold/confidences.csv")
        save_alpha_chain_lengths(confidence_df, low_confidence)
elif algorithm == "rcsb":
    print("Getting RCSB length ranges.")
    pdbs = glob.glob("../data/rcsb/pdb_files/*.ent")
    save_rcsb_chain_lengths(pdbs)
else:
    print("An error occurred. Check arguments.")
    