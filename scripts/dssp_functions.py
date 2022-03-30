# Functions for DSSP code
import subprocess
import numpy as np
import os
import protein_contact_map
import pandas as pd
import collections
import traceback


def get_secondary_structure(filename: str, dssp_path: str, log_file: str) -> np.ndarray:
    """
    H = α-helix
    B = residue in isolated β-bridge
    E = extended strand, participates in β ladder
    G = 3-helix (310 helix)
    I = 5 helix (π-helix)
    T = hydrogen bonded turn
    S = bend

    @param filename: pdb file + path
    @param dssp_path: full path to dssp executable
    @param log_file: text file to record any exceptions
    @return: numpy array of secondary structures in pdb file
    """
    with open(log_file, "w") as log_file:
        secondary_structures = []
        try:
            subprocess.check_call(f"{dssp_path} {filename} -o result.dssp", shell=True)
            infile = open("result.dssp", "r")
            readit = False
            for line in infile:
                if readit:
                    try:
                        if line[13:15] == '!*' or line[13] == '!':
                            continue
                        else:
                            secondary_structure = line[16]
                            if line[16] == " ":
                                secondary_structure = "-"
                            secondary_structures.append(secondary_structure)
                    except:
                        continue
                if "#" in line:
                    readit = True
            infile.close()
            os.remove("result.dssp")
        except subprocess.CalledProcessError as e:
            traceback.print_exc(file=log_file)
        except Exception as e:
            traceback.print_exc(file=log_file)
    return np.array(secondary_structures)


def save_alpha_secondary_info():
    """
    Get secondary structure and confidence information for AlphaFold 2 from DSSP and save into csv
    @return:
    """
    length_ranges = ["100", "200", "300"]
    secondary_structures = []
    chain_counter = 1
    for length in length_ranges:
        print(f"Chains {chain_counter}/{len(length_ranges)}")
        pdb_counter = 1
        confidence_data = pd.read_csv(f"../data/alphafold/confidences_{length}.csv")
        pdbs = confidence_data["filename"].to_numpy()
        for pdb in pdbs:
            print(f"At {pdb_counter}/{len(pdbs)}")
            secondary_structure = get_secondary_structure(pdb, "/usr/bin/dssp", "log.txt")
            counter = collections.Counter(secondary_structure)
            chain_length = len(secondary_structure)
            for key in counter:
                counter[key] = counter[key] / chain_length
            secondary_structures.append(counter)
            pdb_counter += 1
        secondary_df = pd.DataFrame.from_dict(secondary_structures).fillna(0)
        secondary_df.to_csv(f"../data/alphafold/structures_{length}_raw.csv")
        chain_counter += 1
        secondary_df["filename"] = confidence_data["filename"]
        secondary_df["mean"] = confidence_data["mean_conf"]
        secondary_df["std"] = confidence_data["std_conf"]
        secondary_df["beta"] = secondary_df["B"] + secondary_df["E"]
        clean_df = secondary_df[["filename", "mean", "std", "H", "beta"]]
        renamed_cols_df = clean_df.rename(columns={"H": "alpha"})
        renamed_cols_df.to_csv(f"../data/alphafold/secondary_structures_{length}.csv")


def save_rcsb_secondary_info():
    """
    Get secondary structure information for RCSB from DSSP and save into csv
    @return None:
    """
    length_ranges = ["100", "200", "300"]
    secondary_structures = []
    chain_counter = 1
    for length in length_ranges:
        print(f"Chains {chain_counter}/{len(length_ranges)}")
        pdb_counter = 1
        chain_data = pd.read_csv(f"../data/rcsb/chains_{length}.csv")
        pdbs = chain_data["filename"].to_numpy()
        for pdb in pdbs:
            print(f"At {pdb_counter}/{len(pdbs)}")
            secondary_structure = get_secondary_structure(pdb, "/usr/bin/dssp", "log.txt")
            counter = collections.Counter(secondary_structure)
            chain_length = len(secondary_structure)
            for key in counter:
                counter[key] = counter[key] / chain_length
            secondary_structures.append(counter)
            pdb_counter += 1
        secondary_df = pd.DataFrame.from_dict(secondary_structures).fillna(0)
        secondary_df["filename"] = chain_data["filename"]
        secondary_df.to_csv(f"../data/rcsb/structures_{length}_raw.csv")
        chain_counter += 1
        secondary_df["beta"] = secondary_df["B"] + secondary_df["E"]
        clean_df = secondary_df[["filename", "H", "beta", "B", "E"]]
        renamed_cols_df = clean_df.rename(columns={"H": "alpha"})
        renamed_cols_df.to_csv(f"../data/rcsb/secondary_structures_{length}.csv")
