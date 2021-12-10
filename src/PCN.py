"""
Base class for PCNs. 

"""
import networkx as nx
import MDAnalysis as mda
import numpy as np
import random

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

    def get_C_alphas(self):
        
        """
        Take a universe, get segment A and return the C-alphas for the PDB. 

        Parameters
        ----------
        universe: MDAnalysis.Universe(PDB_file)

        Return
        ------
        C_alphas: MDAnalysis AtomGroup
            contains the C-alphas from the PDB
        """
        "------------------------MOVE-------------------------------"
        segments = self.universe.residues.segments
            
        number_of_segments = len(segments)

        for index in range(number_of_segments):
            segment_id = self.universe.residues.segments.segids[index]

            try:
                C_alphas = self.universe.select_atoms(f'name CA and segid {segment_id}')
                #print(segment_id)
                break
            
            except:
                print('Error in first segment id. Trying the next one.')
                continue

        "-----------------------------------------------------------"
        return C_alphas

    def get_chain_length(self, C_alphas):
        """
        Take C-alphas and use that to return the length of the chain.

        Parameters
        ----------
        C_alphas: MDAnalysis AtomGroup
            contains the C-alphas from the PDB

        Return
        ------
        chain_length: int
            length of the chain
        """
        chain_length = len(C_alphas)
        
        return chain_length

    def create_connected_component_subgraphs(self, protein_graph):
        """
        Take a protein graph and create connected component subgraphs.
        
        Parameters
        ----------
        protein_graph: nx.graph()
            protein graph containing the links

        Return
        ------
        cc_subgraph: nx.subgraph() generator
            subgraph of connected components in the protein graph
        """
        for components in nx.connected_components(protein_graph):
            yield protein_graph.subgraph(components)

    def get_link_lengths(self, C_alphas):
        """
        Use the C-alphas to calculate the link lengths between adjacent atoms.
        
        Get the positions of C-alphas, create a proteingraph and compute link 
        engths.

        The link lengths are calculated from an adjacency matrix, where links 
        are added if the distance between adjacent atoms is less than the 
        threshold.

        Parameters
        ----------
        C_alphas: MDAnalysis AtomGroup
            contains the C-alphas from the PDB

        Return
        ------
        link_lengths: list
            list containing lengths of links 
        """

        link_lengths = []
        C_alpha_positions = C_alphas.positions
        protein_graph = nx.empty_graph(len(C_alpha_positions))

        for i in range(len(C_alpha_positions)):
            for j in range(i):
                # Get the distance between two adjacent atoms
                distance = np.linalg.norm(C_alpha_positions[i] - C_alpha_positions[j])
                if distance < self.link_length_threshold:
                    protein_graph.add_edge(i,j)
        # Add links to list
        if len(protein_graph.edges()) > 0:
            protein_graph = list(self.create_connected_component_subgraphs(protein_graph))[0]
            number_protein_nodes = len(protein_graph.nodes())
            links = list(set(protein_graph.edges()) - set(nx.cycle_graph(number_protein_nodes).edges()))
            for link in links:
                link_lengths.append(abs(link[0] - link[1]))

        return link_lengths







