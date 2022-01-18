
import numpy as np
import pandas as pd


filepath = "/Users/jasminguven/Desktop/"
sequence_lengths = ["100", "200", "300"]
uniprot_id_list = []
for sequence_length in sequence_lengths:
    uniprot_ids_dataframe = pd.read_excel(f"{filepath}/test_{sequence_length}.xlsx")
    uniprot_id_list.append(uniprot_ids_dataframe["Entry"].values)
uniprot_ids = np.asarray(uniprot_id_list)
print(uniprot_ids.flatten())

# sequence_length = '100'
# uniprot_IDs = functions.get_uniprot_ids(sequence_length)
#
# download = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_IDs[0]}-F1-model_v2.pdb'
#
# try:
#     file_name = f'../data/alphafold_data/temp/AF-{uniprot_IDs[0]}-F1-model_v2.pdb'
#     #bio.retrieve_pdb_file(pdb,pdir='../data/tmep/', file_format='pbd')
#     urllib.request.urlretrieve(download, file_name)
# except Exception:
#     print("failed at %s" % pdb, flush=True)
#     traceback.print_exc(file=log)
#
# u = mda.Universe(file_name)
# print(u)

# PDB_IDs = ['5EUB']

# for PDB_ID in PDB_IDs:
#     PDB_ID_lower = str(PDB_ID).lower()
#     PDB_file = f'/Volumes/Seagate_Extension_Plus/PDBs/pdb{PDB_ID_lower}.ent'
#     if os.path.isfile(PDB_file):
#         u = mda.Universe(PDB_file)


#         # for seg_index in range(number_of_segments):

#         segments = u.universe.residues.segments
            
#         number_of_segments = len(segments)

#         for index in range(number_of_segments):
#             segment_id = u.universe.residues.segments.segids[index]

#             try:
#                 alpha_carbons = u.universe.select_atoms(f'name CA and segid {segment_id}')
                
#                 print(len(alpha_carbons))
#                 break
            
#             except:
#                 print('Error in first segment id. Trying the next one.')
#                 continue



        # if segment_id:
        #     print(segment_id)
        # segments = self.universe.residues.segments
        # number_of_segments = len(segments)
        # for seg_index in range(number_of_segments):
        #     str_seg_ids = str(segments[seg_index])
            
        #     if '>' != str_seg_ids[9]:
        #          segment_id = segments[seg_index]


