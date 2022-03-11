"""
Base class for PCNs. 

"""
import MDAnalysis
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


def get_chain_length(alpha_carbons: MDAnalysis.AtomGroup) -> int:
    """
    Calculate the chain length of protein from alpha carbons.
    @param alpha_carbons
    @return: length of chain
    """
    chain_length = 0
    if alpha_carbons:
        chain_length = len(alpha_carbons)
    return chain_length


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
        @return:
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

            except:
                print("Error in first segment id. Trying the next one.")
                continue
        return alpha_carbons

    def get_link_lengths(self, alpha_carbons: MDAnalysis.AtomGroup) -> np.ndarray:
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

    def get_protein_graph(self, alpha_carbons: MDAnalysis.AtomGroup) -> nx.Graph:
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
