# Imports
import numpy as np
import networkx as nx
import traceback
import os
from Bio.PDB import *
from MDAnalysis import *
import warnings
import urllib
import pandas as pd
import socket
import random
import operator
warnings.filterwarnings('ignore')

def concat_ID_files():
	''' Opens id_file_N.txt as a DF and concats them all into one DF '''
	frames = []
	for filename in os.listdir('../data/ids/'):
		if filename.endswith(".txt"):
			temp_df = pd.read_csv('../data/ids/'+str(filename), header=None)
			frames.append(temp_df)

	DF = pd.concat(frames).reset_index()
	DF = DF.drop(columns='index')
	PDB_list = DF[0].to_numpy()
	return PDB_list

def make_connected_component_subgraphs(graph):
	'''Take a graph and get subgraphs.'''
	for components in nx.connected_components(graph):
		yield graph.subgraph(components)


def get_link_lengths(CA):
	'''Get link lengths and return a list.'''
	# Define the threshold for connections
	threshold = 8.0 #Å

	link_length_list = []
	linklengths = []

	# Get positions of C alphas
	r = CA.positions
	# Create an empty graph
	graph = nx.empty_graph(len(r))

	# Define max and average distances
	max_dist = 0
	avg_dist = 0
	'''-------------------------------???-------------------------------'''
	# Create the adjacency matrix
	for i in range(len(r)):
		for j in range(i):
			# Get the norm, i.e. distnces between adjacent atoms
			distance = np.linalg.norm(r[i] - r[j])
			max_dist = max(distance, max_dist)
			avg_dist += distance ** 2
			# If closer than threshold, add link
			if distance < threshold:
				graph.add_edge(i,j)
				linklengths.append(distance)
	'''-------------------------------???-------------------------------'''
	# Get the average distance
	avg_dist = np.sqrt(avg_dist) / (len(r))
	# Take the largest component
	if len(graph.edges()) > 0:
		# Get the connected component subgraphs
		graph = list(make_connected_component_subgraphs(graph))[0]
		# Get number of nodes
		nodes = len(graph.nodes())
		# Make a list of links ?
		links = list(set(graph.edges()) - set(nx.cycle_graph(nodes).edges()))
		# Loop over links
		for link in links: 
			# Get link length
			link_length_list.append(abs(link[0] - link[1]))

	return link_length_list


def bootsrap(ll_list, link_length_id_list):
	'''Bootstrap LL data.'''
	# Specify a sample size to use for bootstrapping
	sample_size = 1000
	data_histogram_list = []
	mean_list = []
	std_list = []
	# Loop over number of samples
	for n in range(sample_size):
		# Randomly choose a list of observations from the sample
		bts_sample = [random.choice(ll_list) for _ in ll_list]
		# Sample with replacement
		bts_sample = [item for items in bts_sample for item in items]
		# Make a histogram
		data_histogram, data_histogram_edges = np.histogram(bts_sample, bins = range(3, 1000), density = True)
		# Normalise the histogram? But doesn't density already do this?
		data_histogram = data_histogram / len(link_length_id_list)
		data_histogram_list.append(data_histogram)
	data_histogram_list = np.array(data_histogram_list)
	mean_list = data_histogram_list.mean(axis = 0)
	std_list = data_histogram_list.std(axis = 0)





# Get pdbs with residue lengths in ranges [85-115], [185-215] and [285-315]
# based on number of C-alpha in PDBs, from entire PDB
if __name__ == '__main__':
	
	# Get PDB ids 
	PDBs = concat_ID_files()

	# Create lists for residue length ranges
	hundred = [] # 85-115
	twohundered = [] # 185-215
	threehundred = [] # 285-315
	# Create lists for link length ranges
	ll_100 = []
	ll_200 = []
	ll_300 = []

	counter = 0

	# Number of PDB ids we want to retrieve 
	#n_pdbs = len(query_DF)
	
	# print(n_pdb)
	#n_pdb = 1000 
	#PDBs = ['3UQI']

	# Open log file & write to it
	with open("log.txt", "w") as log:
	# 	# Loop over PDB ids
		for pdb in PDBs:
			print(f'We\'re at PDB ID {pdb}.')
			# Make url request
			download = f'https://files.rcsb.org/download/{pdb}.pdb' 
			
			try:
				PDB_file = f'../data/temp/{pdb}.pdb'
				urllib.request.urlretrieve(download,PDB_file)
			
			except:
				print(f'Failed at PDB ID: {pdb}')
				traceback.print_exc(file=log) # Log the exception
				#continue # Jump back to the top of the loop and try another PDB

			# MDAnalysis: 
			if os.path.isfile(PDB_file):
				# Create a new Universe
				u = Universe(PDB_file) 
				# Select SEG ID A 
				segid_A = u.residues.segments.segids[0] 
				# Select C-alphas
				CA = u.select_atoms(f'name CA and segid {segid_A}')
				# Get chain length from number of C-alphas
				chain_length = len(CA)
				# Split by chain lengths
				if chain_length in range(85, 116): 
					# Save PDB ID in 100s file
					hundred.append(pdb)
					# Get link lengths
					ll_100.append(get_link_lengths(CA))
					#print(f'this is ll_100: {ll_100}')
				elif chain_length in range(185, 216):
					# Save PDB ID in 100s file
					twohundered.append(pdb)
					# Get link lengths
					ll_200.append(get_link_lengths(CA))
					#print(f'this is ll_200: {ll_200}')
				elif chain_length in range(285, 316):
					# Save PDB ID in 100s file
					threehundred.append(pdb)
					# Get link lengths
					ll_300.append(get_link_lengths(CA))
					#print(f'this is ll_300: {ll_300}')
					
				# Delete the file after we're done
				os.remove(PDB_file)
				counter += 1
			else:
				continue

			'''----------------------------------------------------------'''
			# print('Inside first FOR loop')
			# if counter % 5000 == 0:
			# 	print(f"Counter after mod: {counter}")
			# 	print(f"")
			# 	# print('Passed second IF at %i' %counter)
			# 	# print("We are at entry %d/%d!" % (counter, n_pdbs), flush=True)
			# 	# save_bootstrapped_current_data('100', ll_100, hundred)
			# 	counter += 1
			'''----------------------------------------------------------'''






				 