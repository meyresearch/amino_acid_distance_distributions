# Functions for DSSP code
import subprocess
import numpy as np
import os
import protein_contact_map
import pandas as pd
import collections

    
def get_secondary_structure(filename: str, dssp_path: str) -> np.ndarray:
    
    """
    H = α-helix
    B = residue in isolated β-bridge
    E = extended strand, participates in β ladder
    G = 3-helix (310 helix)
    I = 5 helix (π-helix)
    T = hydrogen bonded turn
    S = bend
    """
    
    #call DSSP
    try:
        subprocess.check_call(f"{dssp_path} {filename} -o result.dssp",shell=True)
        infile = open("result.dssp","r")
    except Exception as e:
        raise Exception(f"Could not calculate secondary structure! {e}")

    #parse output
    readit = False
    secondary_structures = []
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

    # clean temporary files
    os.remove("result.dssp")

    return np.array(secondary_structures)


def save_secondary_info():
    
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
            secondary_structure = get_secondary_structure(pdb, "/usr/bin/dssp")
            counter = collections.Counter(secondary_structure)
            chain_length = len(secondary_structure)
            for key in counter:
                counter[key] = counter[key]/chain_length
            secondary_structures.append(counter)
            pdb_counter += 1
              
        secondary_df = pd.DataFrame.from_dict(secondary_structures)
        secondary_df.to_csv(f"../data/alphafold/structures_{length}_raw.csv")
        chain_counter += 1
        
        secondary_df["filename"] = confidence_data["filename"]
        secondary_df["mean"] = confidence_data["mean_conf"]
        secondary_df["std"] = confidence_data["std_conf"]
        secondary_df["beta"] = secondary_df["B"] + secondary_df["E"]
        secondary_df = secondary_df[["filename","mean","std","H","beta"]]
        
        secondary_df.to_csv(f"../data/alphafold/secondary_structures_{length}.csv")
        

