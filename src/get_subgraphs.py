import PCN
import functions
import numpy as np
import networkx as nx
import os 


PDB_IDs = functions.concat_ID_files()
number_PDBs = len(PDB_IDs)

counter = 0
exclude_counter = 0

exclude_IDs = np.load('../data/exclude.npy')

subgraphs_100 = []
subgraphs_200 = []
subgraphs_300 = []
ids_100 = []
ids_200 = []
ids_300 = []
for i in range(number_PDBs):

    if PDB_IDs[i] not in exclude_IDs:
        print('--------------------------------------------------------------------------------------------------')
        print(f'At entry {counter+1} out of {number_PDBs+1} entries.')
        print(f'The PDB ID is {PDB_IDs[i]}')

        try:
            # PDB_tarfile = tarfile.open(f'../data/PDBs.tar.gz', 'r:gz')
            PDB_IDs_lower = str(PDB_IDs[i]).lower()
            # PDB_file = PDB_tarfile.extractfile(f'../data/PDBs/pdb{PDB_IDs_lower}.ent')
            PDB_file = f'/Volumes/Seagate_Extension_Plus/PDBs/pdb{PDB_IDs_lower}.ent'
            PDB_check = open(PDB_file, 'r')
            print('Successfully opened file.')
        except:
            print('Could not open file.') # Log the exception
            continue # Jump back to the top of the loop and try another PDB

        if os.path.isfile(PDB_file):
            # Create the PCN instance
            protein_contact_network = PCN.PCN(PDB_file)
            alpha_carbons = protein_contact_network.get_alpha_carbons()
            chain_length = protein_contact_network.get_chain_length(alpha_carbons)
            # Split into chain length ranges
            if chain_length in range(85,116):
                print('This PDB is in the length range 85-115.')
                alpha_carbon_positions = alpha_carbons.positions
                n_alpha_carbons = len(alpha_carbons)
                parent_graph = nx.empty_graph(n_alpha_carbons)

                for i in range(n_alpha_carbons):
                    for j in range(i):
                        # Get the distance between two adjacent atoms
                        distance = np.linalg.norm(alpha_carbon_positions[i] - alpha_carbon_positions[j])
                        if distance < protein_contact_network.threshold:
                            parent_graph.add_edge(i,j)
                n_parent_graph_edges = len(parent_graph.edges())
                if n_parent_graph_edges > 0:
                    subgraphs = list(protein_contact_network.create_connected_component_subgraphs(parent_graph))
                    n_subgraphs = len(subgraphs)
                    
                    if n_subgraphs > 1: # more than one subgraph
                        print('More than one subgraph')
                        subgraphs_100.append(subgraphs)
                        ids_100.append(PDB_IDs[i])

            if chain_length in range(185,216):
                print('This PDB is in the length range 185-215.')
                alpha_carbon_positions = alpha_carbons.positions
                n_alpha_carbons = len(alpha_carbons)
                parent_graph = nx.empty_graph(n_alpha_carbons)

                for i in range(n_alpha_carbons):
                    for j in range(i):
                        # Get the distance between two adjacent atoms
                        distance = np.linalg.norm(alpha_carbon_positions[i] - alpha_carbon_positions[j])
                        if distance < protein_contact_network.threshold:
                            parent_graph.add_edge(i,j)
                n_parent_graph_edges = len(parent_graph.edges())
                if n_parent_graph_edges > 0:
                    subgraphs = list(protein_contact_network.create_connected_component_subgraphs(parent_graph))
                    n_subgraphs = len(subgraphs)
                    
                    if n_subgraphs > 1: # more than one subgraph
                        print('More than one subgraph')
                        subgraphs_200.append(subgraphs)
                        ids_200.append(PDB_IDs[i])
           
            if chain_length in range(285,316):
                print('This PDB is in the length range 285-315.')
                alpha_carbon_positions = alpha_carbons.positions
                n_alpha_carbons = len(alpha_carbons)
                parent_graph = nx.empty_graph(n_alpha_carbons)

                for i in range(n_alpha_carbons):
                    for j in range(i):
                        # Get the distance between two adjacent atoms
                        distance = np.linalg.norm(alpha_carbon_positions[i] - alpha_carbon_positions[j])
                        if distance < protein_contact_network.threshold:
                            parent_graph.add_edge(i,j)
                n_parent_graph_edges = len(parent_graph.edges())
                if n_parent_graph_edges > 0:
                    subgraphs = list(protein_contact_network.create_connected_component_subgraphs(parent_graph))
                    n_subgraphs = len(subgraphs)
                    
                    if n_subgraphs > 1: # more than one subgraph
                        print('More than one subgraph')
                        subgraphs_300.append(subgraphs)
                        ids_300.append(PDB_IDs[i])
            else:
                print('PDB is outside chain length ranges.')
                counter += 1

np.save(f'../data/subgraph_data/subgraphs_100.npy', subgraphs_100)
np.save(f'../data/subgraph_data/subgraphs_200.npy', subgraphs_200)
np.save(f'../data/subgraph_data/subgraphs_300.npy', subgraphs_300)

np.save(f'../data/subgraph_data/ids_100.npy', ids_100)
np.save(f'../data/subgraph_data/ids_200.npy', ids_200)
np.save(f'../data/subgraph_data/ids_300.npy', ids_300)



