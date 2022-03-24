# Functions for DSSP code
import subprocess
import numpy as np
import os

def get_secondary_structure(filename, dssp_path):
    
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
