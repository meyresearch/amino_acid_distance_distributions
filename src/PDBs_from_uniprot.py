
import functions
import urllib

sequence_lengths = ['100', '200', '300']

for sequence_length in sequence_lengths:
    counter = 0

    uniprot_IDs = functions.get_uniprot_IDs(sequence_length)
    for uniprot_ID in uniprot_IDs:
        print(f'At entry {counter}/{len(uniprot_IDs)}')
        download = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_ID}-F1-model_v2.pdb'

        try:
            file_name = f'../data/alphafold_data/AF_PDBs/AF-{uniprot_ID}-F1-model_v2.pdb'
            #bio.retrieve_pdb_file(pdb,pdir='../data/tmep/', file_format='pbd')
            urllib.request.urlretrieve(download, file_name)
        except Exception:
            
            'error'
        counter += 1
         
