# RCSB PDB Data

This directory contains the RCSB PDB sequence distance distribution data. All code is available at [`../scripts`](https://github.com/meyresearch/sequence_distance_distribution/tree/36a848c5b15a62dc646de583f03c824f680873da/scripts).

## Table of contents
### [1. Download PDB IDs](#1-download-pdb-ids-1)
### [2. Retrieve RCSB PDB structures](#2-retrieve-rcsb-pdb-structures-1)
### [3. Sort files into chain length ranges](#3-sort-files-into-chain-length-ranges-1)
### [4. Get secondary structure information with DSSP](#4-get-secondary-structure-information-with-dssp-1)
### [5. Compute amino acid distances](#5-compute-amino-acid-distances-1)
### [6. Plot](#6-plot-1)

## How to get the data

### 1. Download and handle PDB IDs

Run `get_current_entry_ids.sh`. This will save all the current RCSB PDB entry IDs in `combined_ids.txt`. 

Next, run `python transpose_csv.py`. This reads in the `combined_ids.txt` and saves the IDs in a single column format
in a csv file called `combined_ids_tr.csv`.

Finally, run `python chunk_id_files.py` which reads in the `combined_ids_tr.csv` and separates the IDs into file chunks
of 1000 IDs each. The files are saved in `ids/id_file_<i>.txt`. 

### 2. Retrieve RCSB PDB structures

Go into [`../scripts/trivial_parallelisation/`](https://github.com/meyresearch/sequence_distance_distribution/tree/bump/scripts/trivial_parallelisation)
and run `make_directories.sh`, which will create 15 subdirectories within `ids/` called `id_file_<i>`. You can then divide the chunked id files 
into the `id_file_<i>` directories. The code will also copy `retrieve_pdbs_parallel.py` into each subdirectory.

Finally, run `embarrassing_parallel.sh` to simultaneously run `retrieve_pdbs_parallel.py` in each subdirectory. The
code retrieves the RCSB PDB structures and saves them into `pdb_files/`.

Optionally, you can directly run `retrieve_rcsb_pdbs.py`, but this will take significantly longer to run.

### 3. Sort files into chain length ranges

Run `python save_chain_lengths.py rcsb` to Sort the PDBs into the three chain length ranges. The data is saved in 
`chains_<chain-length-range>.csv`. 

### 4. Get secondary structure information with DSSP

> ### ❗️Note❗️
> This step requires DSSP to be installed on your machine. 
> For more info and instructions, please see [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/).

Run `python dssp.py rcsb` to  output secondary structure information on each PDB in each range. 
The data is saved in `secondary_structures_<chain-length-range>.csv` where `<chain-length-range>` is `100`, `200` or `300`. 
The code loops through all three ranges automatically, so you don't have to input this yourself.

### 5. Compute amino acid distances

Finally, run `python run.py rcsb -r <chain-length-range> -p <../data/rcsb/secondary_structures_{range}.csv>` for each
range to compute the amino acid distances. This returns a histogram of the amino distance distribution in each range. 
These are saved in `histogram_<chain-length-range>_not_normed.npy`. 

### 6. Plot

For plotting see [`../plots/`](https://github.com/meyresearch/sequence_distance_distribution/tree/bump/plots). 
