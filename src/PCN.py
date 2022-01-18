"""
Base class for PCNs. 

"""
import networkx as nx
import MDAnalysis as mda
import numpy as np
import random
import warnings
warnings.filterwarnings("ignore")
import pandas as pd

class PCN():
    """
    Base class for PCNs. 
    
    Takes a PDB file and creates a PCN from the alpha carbons.

    """
    
    def __init__(self, PDB_file, default_threshold = 8.0):
        """
        Class constructor for PCN(). Creates an MDAnalysis Universe from supplied PDB file.
        Initialises the link length threshold for the PCN object. 
        Default threshold value is 8.0 Å.
        
        Parameters
        ----------
        PDB_file: str
            Exact location and name of the PDB file. 

        default_threshold: float
            Link length threshold value. 

        Return
        ------
        None
        """
        self.universe = mda.Universe(PDB_file)
        self.threshold = default_threshold


    @property
    def threshold(self):
        """
        Get the current link length threshold value. 

        Parameters
        ----------
        None

        Return
        ------
        self._threshold: float
            Current link length threshold.
        """
        return self._threshold
    
    @threshold.setter
    def threshold(self, threshold_value):
        """
        Set the link length threshold value. 

        Parameters
        ----------
        threshold_value: float
            Link length threshold value in Å.

        Return
        ------
        None
        """
        if threshold_value != float(threshold_value):
            raise TypeError('The threshold value must be a float.')

        if threshold_value >= 0:
            self._threshold = threshold_value

        else:
            raise ValueError('The threshold value must be larger than 0.')
            
    @threshold.deleter
    def threshold(self):
        """
        'Deleter' for link length threshold. Will raise an error if deleted. 

        Parameters
        ----------
        None

        Return
        ------
        None
        """
        raise AttributeError('The threshold value cannot be deleted. You can set it to 0.')

    def get_alpha_carbons(self):
        
        """
        Take a universe, get segment A and return the C-alphas for the PDB. 

        Parameters
        ----------
        universe: MDAnalysis.Universe(PDB_file)

        Return
        ------
        alpha_carbons: MDAnalysis AtomGroup
            contains the C-alphas from the PDB
        """
        "------------------------MOVE-------------------------------"
        segments = self.universe.residues.segments
            
        n_segments = len(segments)

        for i in range(n_segments):
            segment_id = self.universe.residues.segments.segids[i]

            try:
                protein_atoms = self.universe.select_atoms('protein')
                MDA_altlocs = protein_atoms.altLocs
                alternative_locations = []
                for loc in MDA_altlocs:
                    if loc != '':
                        alternative_locations.append(loc)
                # Remove duplicates and sort alphabetically
                alternative_locations = np.sort(list(set(alternative_locations)))
                n_altlocs = len(alternative_locations)

                if n_altlocs > 1:
                    first_altloc = alternative_locations[0]
                    residue_ids = protein_atoms.residues.resids
                    alpha_carbons = self.universe.select_atoms(f'resid {residue_ids[0]}:{residue_ids[-1]} and name CA and segid {segment_id}')

                    for i in range(1, n_altlocs):
                        exclude_atoms = self.universe.select_atoms(f'resid {residue_ids[0]}:{residue_ids[-1]} and name CA and segid {segment_id} and altLoc {alternative_locations[i]}')
                        alpha_carbons -= exclude_atoms
                else:
                    residue_ids = protein_atoms.residues.resids
                    alpha_carbons = self.universe.select_atoms(f'resid {residue_ids[0]}:{residue_ids[-1]} and name CA and segid {segment_id}')

                break
            
            except:
                print('Error in first segment id. Trying the next one.')
                continue

        "-----------------------------------------------------------"
        return alpha_carbons

    def get_chain_length(self, alpha_carbons):
        """
        Take C-alphas and use that to return the length of the chain.

        Parameters
        ----------
        alpha_carbons: MDAnalysis AtomGroup
            contains the C-alphas from the PDB

        Return
        ------
        chain_length: int
            length of the chain
        """
        chain_length = len(alpha_carbons)
        
        return chain_length

    def create_connected_component_subgraphs(self, graph):
        """
        Take a protein graph and create connected component subgraphs.
        
        Parameters
        ----------
        graph: nx.graph()
            protein graph containing the links

        Return
        ------
        cc_subgraph: nx.subgraph() generator
            subgraph of connected components in the protein graph
        """
        for components in nx.connected_components(graph):
            yield graph.subgraph(components)

    def get_link_lengths(self, alpha_carbons, PDB_ID):
        """
        Use the C-alphas to calculate the link lengths between adjacent atoms.
        
        Get the positions of C-alphas, create a proteingraph and compute link 
        engths.

        The link lengths are calculated from an adjacency matrix, where links 
        are added if the distance between adjacent atoms is less than the 
        threshold.

        Parameters
        ----------
        alpha_carbons: MDAnalysis AtomGroup
            contains the C-alphas from the PDB

        Return
        ------
        link_lengths: list
            list containing lengths of links 
        """

        link_lengths = []
        alpha_carbon_positions = alpha_carbons.positions
        n_alpha_carbons = len(alpha_carbons)
        parent_graph = nx.empty_graph(n_alpha_carbons)

        for i in range(n_alpha_carbons):
            for j in range(i):
                # Get the distance between two adjacent atoms
                distance = np.linalg.norm(alpha_carbon_positions[i] - alpha_carbon_positions[j])
                if distance < self.threshold:
                    parent_graph.add_edge(i,j)
        n_parent_graph_edges = len(parent_graph.edges())
        n_parent_graph_nodes = len(parent_graph.nodes())
        # Add links to list
        if n_parent_graph_edges > 0:
            subgraphs = list(self.create_connected_component_subgraphs(parent_graph))
            n_subgraphs = len(subgraphs)
            protein_graph = subgraphs[0] # changed [0]
            n_protein_graph_nodes = len(protein_graph.nodes())
            if n_subgraphs > 1: # more than one subgraph
                print('More than one subgraph. This PDB will be excluded.')
            else: # just one subgraph
                cycle_graph = nx.cycle_graph(n_parent_graph_nodes)
                links = list(set(protein_graph.edges()) - set(cycle_graph.edges()))
                for link in links:
                    link_length = abs(link[0] - link[1])
                    if link_length <= 1:
                        print('There is a link of length 1.')
                    else:
                        link_lengths.append(link_length)   
        return link_lengths







