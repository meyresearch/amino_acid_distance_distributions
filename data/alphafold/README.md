# AlphaFold 2 Data

This directory contains the AlphaFold 2 sequence distance distribution data. All code is available at [`scripts`](https://github.com/meyresearch/sequence_distance_distribution/tree/36a848c5b15a62dc646de583f03c824f680873da/scripts).

## Table of contents
### [1. Download Uniprot IDs](#1-download-uniprot-ids-1)
### [2. Retrieve AlphaFold 2 predicted structures](#2-retrieve-alphafold-2-predicted-structures-1)
### [3. Get AlphaFold 2 per-residue confidence scores](#3-get-alphafold-2-per-residue-confidence-scores-1)
### [4. Sort files into chain length ranges](#4-sort-files-into-chain-length-ranges-1)
### [5. Get secondary structure information with DSSP](#5-get-secondary-structure-information-with-dssp-1)
### [6. Remove structures describing the same protein from different organisms](#6-remove-structures-describing-the-same-protein-from-different-organisms-1)
### [7. Compute amino acid distances](#7-compute-amino-acid-distances-1)
### [8. Plot](#8-plot-1)

## How to get the data

### 1. Download Uniprot IDs

The Uniprot ID files (accessed March 2022) are available in the [uniprot_ids](https://github.com/meyresearch/sequence_distance_distribution/tree/36a848c5b15a62dc646de583f03c824f680873da/data/alphafold/uniprot_ids) 
directory.

To download the current IDs, go to [uniprot.org](https://www.uniprot.org) and select Advanced search. 
In the first selection term, click the "All" drop-down menu and from "Sequence" select "Sequence length".
For each of the length ranges, use the following:

* chain lengths of approx. 100 amino acids: From 85 To 115,
* approx. 200 amino acids: From 185 To 215 and
* approx. 300 amino acids: From 285 To 315.

Then, from the left-hand side, choose "Filter by: Reviewed (SwissProt)". Next, you can remove any unnecessary 
columns by clicking the "Columns" option at the top of the search results and leaving only "Entry name" and 
"Length" selected. Click "Save" and then "Download". From the "Download" option, choose Excel as the format and click "Go". 

### 2. Retrieve AlphaFold 2 predicted structures

From [`scripts`](https://github.com/meyresearch/sequence_distance_distribution/tree/36a848c5b15a62dc646de583f03c824f680873da/scripts), 
run `python retrieve_alphafold_pdbs.py <chain-length-range>`, where `<chain-length-range>` is `100`, `200` or `300`. Optionally,
run `get_alphafold_pdbs.sh` which will call `nohup` and run the above Python code simultaneously for all three ranges.

The code will save the AlphaFold 2 PDB files in [pdb_files/](https://github.com/meyresearch/sequence_distance_distribution/tree/36a848c5b15a62dc646de583f03c824f680873da/data/alphafold/pdb_files).

### 3. Get AlphaFold 2 per-residue confidence scores

Per-residue scores of above 90 were used in this paper. To filter for the confidence scores, run `python confidences.py`. This reads through the PDB files and returns a `.csv` file with PDBs with confidences above 90. The confidence data is saved in `confidences.csv`.

### 4. Sort files into chain length ranges

Run `python save_chain_lengths.py alpha` to read in the confidence data and sort the PDBs into the three chain length ranges. This also saves the mean and standard deviation of the confidence scores in `confidences_<chain-length-range>.csv`. 

### 5. Get secondary structure information with DSSP

> ### ❗️Note❗️
> This step requires DSSP to be installed on your machine. 
> For more info and instructions, please see [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/).

Run `python dssp.py alpha` to again read in the confidence data and output secondary structure information on each PDB in each range. The data is saved in `secondary_structures_<chain-length-range>.csv` where `<chain-length-range>` is `100`, `200` or `300`. The code loops through all three ranges automatically, so you don't have to input this yourself.

### 6. Remove structures describing the same protein from different organisms

Next, run `python remove_duplicates.py` to filter out PDBs that describe the same protein but from different organisms. The *unique* structures are then saved in `unique_secondary_structures_<chain-length-range>.csv` where `<chain-length-range>` is `100`, `200` or `300`. The code loops through all three ranges automatically, so you don't have to input this yourself.

### 7. Compute amino acid distances

Finally, run `python run.py alpha -r <chain-length-range> -p <../data/alphafold/unique_secondary_structures_{range}.csv>` for each range to compute the amino acid distances. This returns a histogram of the amino distance distribution in each range. These are saved in `histogram_<chain-length-range>_not_normed.npy`. 

### 8. Plot

For plotting see [`../plots/`](https://github.com/meyresearch/sequence_distance_distribution/tree/36a848c5b15a62dc646de583f03c824f680873da/plots). 



