"""
Get UniProt ids and loop through them to save AlphaFold PDBs in ../data/alphafold/pdb_files/
"""
import urllib.request
from sequence_distance_distribution.scripts import functions

sequence_lengths = ["100", "200", "300"]

for sequence_length in sequence_lengths:
    counter = 0

    uniprot_IDs = functions.get_uniprot_ids(sequence_length)
    for uniprot_ID in uniprot_IDs:
        print(f"At entry {counter}/{len(uniprot_IDs)}")
        download = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_ID}-F1-model_v2.pdb"

        try:
            file_name = f"../data/alphafold/pdb_files/AF-{uniprot_ID}-F1-model_v2.pdb"
            urllib.request.urlretrieve(download, file_name)
        except ValueError:
            print("Error downloading file.")

        counter += 1
         
