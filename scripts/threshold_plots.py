""" Supplemental plots for investigating the effect of the threshold value on amino acid distance
    distributions.
"""

import functions
import numpy as np
import pandas as pd
import protein_contact_map
import matplotlib.pyplot as plt

def pdb_to_adjacency(pdb_file: str, threshold: float) -> tuple:
    """
    Convert given PDB file to an adjacency matrix
    @param pdb_file: PDB file from RCSB or AlphaFold in each length range
    @param threshold: threshold value to use for counting contacts
    @return: tuple of the adjacency and distance matrices as numpy arrays
    """
    pcm = protein_contact_map.ProteinContactMap(pdb_file, default_threshold=threshold)
    alpha_carbons = pcm.get_alpha_carbons
    distance_array = protein_contact_map.get_distance_array(alpha_carbons)
    distance_matrix = protein_contact_map.get_distance_matrix(alpha_carbons, distance_array)
    adjacency_matrix = pcm.get_adjacency_matrix(alpha_carbons, distance_array)
    return distance_matrix, adjacency_matrix


   
        
