"""
Get UniProt ids and loop through them to save AlphaFold PDBs in ../data/alphafold/pdbs/
"""
import urllib.request
import functions
import argparse


def commandline_arguments() -> argparse.Namespace:
    """
    Parser for commandline arguments
    @return: commandline arguments from user
    """

    parser = argparse.ArgumentParser(description="Get AlphaFold 2 PDBs.")
    parser.add_argument("range", type=str, choices=["100", "200", "300"], help="Chain length range")
    cl_arguments = parser.parse_args()
    return cl_arguments


sequence_length = commandline_arguments().range

counter = 0

uniprot_IDs = functions.get_uniprot_ids(sequence_length)
for uniprot_ID in uniprot_IDs:
    print(f"At entry {counter}/{len(uniprot_IDs)}")
    print(f"ID: {uniprot_ID}")
    download = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_ID}-F1-model_v2.pdb"
        
    try:
        file_name = f"../data/alphafold/pdb_files/AF-{uniprot_ID}-F1-model_v2.pdb"
        urllib.request.urlretrieve(download, file_name)
    except urllib.error.HTTPError:
        print("No such file.")

        counter += 1
         
