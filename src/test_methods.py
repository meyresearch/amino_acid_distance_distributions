import random
import numpy as np
import pandas as pd
import PCN
import os
import MDAnalysis as mda
PDB_IDs = ['5EUB']

for PDB_ID in PDB_IDs:
	PDB_ID_lower = str(PDB_ID).lower()
	PDB_file = f'/Volumes/Seagate_Extension_Plus/PDBs/pdb{PDB_ID_lower}.ent'
	if os.path.isfile(PDB_file):
		u = mda.Universe(PDB_file)


		# for seg_index in range(number_of_segments):

		segments = u.universe.residues.segments
			
		number_of_segments = len(segments)

		for index in range(number_of_segments):
			segment_id = u.universe.residues.segments.segids[index]

			try:
				C_alphas = u.universe.select_atoms(f'name CA and segid {segment_id}')
				
				print(len(C_alphas))
				break
			
			except:
				print('Error in first segment id. Trying the next one.')
				continue




		# if segment_id:
		# 	print(segment_id)
		# segments = self.universe.residues.segments
		# number_of_segments = len(segments)
		# for seg_index in range(number_of_segments):
		# 	str_seg_ids = str(segments[seg_index])
			
		# 	if '>' != str_seg_ids[9]:
		# 	 	segment_id = segments[seg_index]


