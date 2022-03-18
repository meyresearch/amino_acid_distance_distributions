# Getting RCSB PDB data

## 1. Get all current PDB IDs 

Run `source get_current_entry_ids.sh` to save all current PDB entry IDs into `data/rcsb/combined_ids.txt`.

Alternatively, go to https://data.rcsb.org/rest/v1/holdings/current/entry_ids and manually save the files. 

## 2. Read in and chunk ID files

First, run `scripts/transpose_csv.py` to transpose the `combined_ids.txt` into one column and save as `data/rcsb/combined_ids_tr.csv`.

Then run `/scripts/chunk_id_files.py` to read in the transposed ID file and chunk the IDs into files of 1000 IDs. 
These are saved in `/data/rcsb/ids/`. 

## 3. Retrieve PDB files from RCSB

Run `scripts/retrieve_rcsb_pdbs.py` to retrieve the PDB files in `data/rcsb/pdb_files/`. 

**Note:** This may take a long time to run. Please see the `trivial_parallelisation/` directory for an option to use 
trivial parallelisation to retrieve the PDB files. 