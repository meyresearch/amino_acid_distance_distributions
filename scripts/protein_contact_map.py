"""
Base class for PCMs.

"""
import networkx as nx
import MDAnalysis as mdA
import numpy as np
import MDAnalysis.analysis.contacts
import MDAnalysis.analysis.distances


class ProteinContactMap:
    """
    Base class for ProteinContactNetworks.
    """

    def __init__(self, pdb_file, default_threshold=8.0):
        """
        Construct the ProteinContactMap and initialise default threshold.
        @param pdb_file: str
        @param default_threshold: float
        """
        self.universe = mdA.Universe(pdb_file)
        self.threshold = default_threshold

    @property
    def threshold(self):
        """
        Get threshold
        @return:
        """
        return self._threshold

    @threshold.setter
    def threshold(self, threshold_value: float) -> None:
        """
        Set threshold value
        @param threshold_value:
        @return: None
        """
        if threshold_value != float(threshold_value):
            raise TypeError('The threshold value must be a float.')

        if threshold_value >= 0:
            self._threshold = threshold_value

        else:
            raise ValueError('The threshold value must be larger than 0.')

    @threshold.deleter
    def threshold(self) -> None:
        """
        Deleter for threshold
        @return: None
        """
        raise AttributeError("The threshold value cannot be deleted. You can set it to 0.")

    @property
    def get_alpha_carbons(self) -> None:
        """
        Take a universe, get segment A and return the C-alphas for the PDB.
        @return: None
        """
        segments = self.universe.residues.segments
        n_segments = len(segments)
        alpha_carbons = None
        for i in range(n_segments):
            segment_id = self.universe.residues.segments.segids[i]
            try:
                protein_atoms = self.universe.select_atoms('protein')
                mda_altlocs = protein_atoms.altLocs
                alternative_locations = []
                for loc in mda_altlocs:
                    if loc != '':
                        alternative_locations.append(loc)
                # Remove duplicates and sort alphabetically
                alternative_locations = np.sort(list(set(alternative_locations)))
                n_altlocs = len(alternative_locations)
                if n_altlocs > 1:
                    residue_ids = protein_atoms.residues.resids
                    alpha_carbons = self.universe.select_atoms(
                        f"resid {residue_ids[0]}:{residue_ids[-1]} and name CA and segid {segment_id}")
                    for j in range(1, n_altlocs):
                        exclude_atoms = self.universe.select_atoms(
                            f"resid {residue_ids[0]}:{residue_ids[-1]} and name CA and segid {segment_id} and altLoc "
                            f"{alternative_locations[j]}")
                        alpha_carbons -= exclude_atoms
                else:
                    residue_ids = protein_atoms.residues.resids
                    alpha_carbons = self.universe.select_atoms(
                        f"resid {residue_ids[0]}:{residue_ids[-1]} and name CA and segid {segment_id}")
                break

            except mdA.SelectionError:
                print("Error in first segment id. Trying the next one.")
                continue
        return alpha_carbons


    def get_contacts_array(self, distance_array: np.ndarray) -> np.ndarray:
        """
        Use the distance array to get contacts array
        @return: contacts array
        """
        return MDAnalysis.analysis.contacts.contact_matrix(distance_array, radius=self.threshold)


    def get_adjacency_matrix(self, alpha_carbons: mdA.AtomGroup, distance_array: np.ndarray) -> np.ndarray:
        """
        Get adjacency matrix from PDB
        @param distance_array: distance array from C-alphas
        @param alpha_carbons: C-alphas from PDB
        @return: adjacency matrix as numpy array
        """
        contacts_array = self.get_contacts_array(distance_array)
        adjacency_array = np.zeros(len(contacts_array))
        for i in range(len(contacts_array)):
            if contacts_array[i]:
                adjacency_array[i] = 1
        n_alpha_carbons = get_chain_length(alpha_carbons)
        pcm_matrix = np.zeros((n_alpha_carbons, n_alpha_carbons))
        n_rows = pcm_matrix.shape[0]
        diagonal = 1
        upper_triangle = np.triu_indices(n_rows, k=diagonal)
        pcm_matrix[upper_triangle] = adjacency_array
        np.fill_diagonal(pcm_matrix, 0)
        return np.where(pcm_matrix, pcm_matrix, pcm_matrix.T)



def get_chain_length(alpha_carbons: mdA.AtomGroup) -> int:
    """
    Calculate the chain length of protein from alpha carbons.
    @param alpha_carbons
    @return: length of chain
    """
    chain_length = 0
    if alpha_carbons:
        chain_length = len(alpha_carbons)
    return chain_length


def get_distance_array(alpha_carbons: mdA.AtomGroup) -> np.ndarray:
    """
    Use the C-alphas to get the distance array
    @return: distance array
    """
    return MDAnalysis.analysis.distances.self_distance_array(alpha_carbons.positions)

