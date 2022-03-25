"""
Base class for PCMs.

"""
import networkx as nx
import MDAnalysis as mdA
import numpy as np
import warnings

warnings.filterwarnings("ignore")


def create_connected_component_subgraphs(graph: nx.Graph) -> None:
    """
    Create connected component subgraphs within protein graph.
    @param graph: protein graph
    @return:
    """
    for components in nx.connected_components(graph):
        yield graph.subgraph(components)


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


def get_distance_matrix(alpha_carbons: mdA.AtomGroup, distance_array: np.ndarray) -> np.ndarray:
    """
    Get distance matrix for PDBs from distance array and C-alphas
    @param: alpha_carbons C-alphas from PDB
    @param: distance_array distance array from C-alphas
    @return: distance_matrix as a numpy array
    """
    n_alpha_carbons = get_chain_length(alpha_carbons)
    distance_matrix = np.zeros((n_alpha_carbons, n_alpha_carbons))
    n_rows = distance_matrix.shape[0]
    diagonal = 1
    upper_triangle = np.triu_indices(n_rows, k=diagonal)
    distance_matrix[upper_triangle] = distance_array
    np.fill_diagonal(distance_matrix, 0)
    for i in range(len(distance_matrix)):
        for j in range(i):
            distance_matrix[i, j] = distance_matrix[j, i]
    return distance_matrix


def get_contacts_array(distance_array: np.ndarray) -> np.ndarray:
    """
    Use the distance array to get contacts array
    @return: contacts array
    """
    return mdA.contacts.self_distance_array(distance_array)


def get_distance_array(alpha_carbons: mdA.AtomGroup) -> np.ndarray:
    """
    Use the C-alphas to get the distance array
    @return: distance array
    """
    return mdA.analysis.distances.self_distance_array(alpha_carbons.positions)


def get_adjacency_matrix(alpha_carbons: mdA.AtomGroup, distance_array: np.ndarray) -> np.ndarray:
    """
    Get adjacency matrix from PDB
    @param distance_array: distance array from C-alphas
    @param alpha_carbons: C-alphas from PDB
    @return: adjacency matrix as numpy array
    """
    contacts_array = get_contacts_array(distance_array)
    adjacency_array = np.zeros(len(contacts_array))
    adjacency_array[contacts_array == True] = 1
    n_alpha_carbons = get_chain_length(alpha_carbons)
    pcm_matrix = np.zeros((n_alpha_carbons, n_alpha_carbons))
    n_rows = pcm_matrix.shape[0]
    diagonal = 1
    upper_triangle = np.triu_indices(n_rows, k=diagonal)
    pcm_matrix[upper_triangle] = adjacency_array
    np.fill_diagonal(pcm_matrix, 1)
    return np.where(pcm_matrix, pcm_matrix, pcm_matrix.T)


def get_discreet_distances(adjacency_matrix: np.ndarray) -> np.ndarray:
    """
    Get the amino acid distances from the adjacency matrix
    @param: adjacency_matrix from PDB
    @return list of discrete distances from adjacency matrix
    """
    distances = []
    rows = len(adjacency_matrix)
    columns = len(adjacency_matrix)
    for row in rows:
        for col in columns:
            if adjacency_matrix[row][col] == 1:
                distance = np.abs(col - row)
                distances.append(distance)
    return np.asarray(distances)


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

    def get_link_lengths(self, alpha_carbons: mdA.AtomGroup) -> np.ndarray:
        """
        Use the C-alphas to calculate the amino acid distances.

        @param alpha_carbons
        @return: array of amino acid distances
        """
        link_lengths = []
        parent_graph = self.get_protein_graph(alpha_carbons)
        n_parent_graph_edges = len(parent_graph.edges())
        n_parent_graph_nodes = len(parent_graph.nodes())
        # Add links to list
        if n_parent_graph_edges > 0:
            subgraphs = list(create_connected_component_subgraphs(parent_graph))
            n_subgraphs = len(subgraphs)
            protein_graph = subgraphs[0]  # changed [0]
            n_protein_graph_nodes = len(protein_graph.nodes())
            if n_subgraphs > 1:  # more than one subgraph
                print("More than one subgraph. This PDB will be excluded.")
            else:  # just one subgraph
                cycle_graph = nx.cycle_graph(n_parent_graph_nodes)
                links = list(set(protein_graph.edges()) - set(cycle_graph.edges()))
                for link in links:
                    link_length = abs(link[0] - link[1])
                    if link_length <= 1:
                        print("There is a link of length 1.")
                    else:
                        link_lengths.append(link_length)
        return np.array(link_lengths)

    def get_protein_graph(self, alpha_carbons: mdA.AtomGroup) -> nx.Graph:
        """
        Create protein graph from alpha carbons.
        @param alpha_carbons
        @return: protein graph
        """
        alpha_carbon_positions = alpha_carbons.positions
        n_alpha_carbons = len(alpha_carbons)
        parent_graph = nx.empty_graph(n_alpha_carbons)

        for i in range(n_alpha_carbons):
            for j in range(i):
                # Get the distance between two adjacent atoms
                distance = np.linalg.norm(alpha_carbon_positions[i] - alpha_carbon_positions[j])
                if distance < self.threshold:
                    parent_graph.add_edge(i, j)
        return parent_graph
