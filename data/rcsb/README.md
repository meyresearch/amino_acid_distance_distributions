# RCSB PDB Data

## 1. Get all current PDB IDs 

Run `source get_current_entry_ids.sh` to save all current PDB entry IDs into `data/rcsb/combined_ids.txt`.

Alternatively, go to https://data.rcsb.org/rest/v1/holdings/current/entry_ids and manually save the files. 

## 2. Read in and chunk ID files

First, run `scripts/transpose_csv.py` to transpose the `combined_ids.txt` into one column and save as 
`data/rcsb/combined_ids_tr.csv`.

Then run `/scripts/chunk_id_files.py` to read in the transposed ID file and chunk the IDs into files of 1000 IDs. 
These are saved in `/data/rcsb/ids/`. 

## 3. Retrieve PDB files from RCSB

Run `scripts/retrieve_rcsb_pdbs.py` to retrieve the PDB files in `data/rcsb/pdb_files/`. 

**Note:** This may take a long time to run. Please see the 
[`trivial_parallelisation/`](https://github.com/meyresearch/sequence_distance_distribution/tree/readmes/scripts/trivial_parallelisation) 
directory for an option to use 
trivial parallelisation to retrieve the PDB files. 

## 4. Get amino acid distance distributions

Amino acid distance distributions are computed by running `python run.py rcsb --range=<chain-length> 
--path=<path-to-pdb-files>` where `<chain-length>` is one of `100, 200` or `300`. This creates Protein Contact Maps 
from each PDB in the given chain length
range and calculates the amino acid distances. The amino acid distances will be saved in 
`data/rcsb/lls_<chain-length>.npy`.

## 5. Bootstrap RCSB amino acid distances

To bootstrap the distance data, run `python run.py boots <chain-length> --file=<data/rcsb/lls_<chain-length>.npy>`. 
This will save the "raw" bootstrapped data in `data/rcsb/bootstrap_<chain-length>_raw.csv` and the mean frequencies
of distances, as well as the confidence level bounds, in `data/rcsb/bootstrap_<chain-length>_stats.csv`. 