"""
Opens each chunked ID file and retrieves PDB files from RCSB to ../data/rcsb/pdb_files/.
"""

from Bio.PDB.PDBList import PDBList
import functions

pdblist = PDBList()
pdb_ids = functions.concatenate_rcsb_id_files()
counter = 1
for pdb in pdb_ids:
    print(f"At entry {counter}")
    pdblist.retrieve_pdb_file(pdb_code=pdb,
                              file_format="pdb",
                              pdir="../data/rcsb/pdb_files/")
    counter += 1

