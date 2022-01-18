"""
bio_retrieve_pdbs and get_alphafold_pdbs

"""
import numpy as np
import pandas as pd
from Bio.PDB.PDBList import PDBList
from concatenate_id_files import concat_id_files


def get_rcsb_pdbs(outpath: str) -> None:
    """
    @param outpath: destination directory to save pdbs, default "../../data/rcsb_data/"
    @return: None
    """
    filepath = "../../data/ids"
    pdbs = concat_id_files(filepath)
    counter = 1
    pdblist = PDBList()
    for pdb in pdbs:
        print(f"At entry {counter}")
        pdblist.retrieve_pdb_file(pdb_code=pdb,
                                  file_format="pdb",
                                  pdir=outpath)
        counter += 1


def get_uniprot_ids() -> np.ndarray:
    """
    @return: array of uniprot ids
    """
    filepath = "../data/alphafold_data"
    sequence_lengths = ["100", "200", "300"]
    uniprot_id_list = []
    for sequence_length in sequence_lengths:
        uniprot_ids_dataframe = pd.read_excel(f"{filepath}/uniprot_{sequence_length}s.xlsx")
        uniprot_id_list.append(uniprot_ids_dataframe["Entry"].values)
    uniprot_ids = np.asarray(uniprot_id_list).flatten()
    return uniprot_ids
