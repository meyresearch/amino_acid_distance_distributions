# AlphaFold 2 PDB Data

## 1. Get UniProt IDs 

Go to https://www.uniprot.org and enter Advanced Search. In the first Term field, click on All -> Sequence -> Sequence 
length. For each chain length range, enter the range as

- **85 to 115**,
- **185 to 215** or
- **285 to 315**

Click Search. In the Filter by option on the left, select Reviewed. Next click on Columns and remove ticks on non-essential columns. 
We used Entry name and Length. Click Save.

Alternatively, use the below links:

- [**100 amino acids**](https://www.uniprot.org/uniprot/?query=length:[85%20TO%20115]&fil=reviewed%3Ayes&sort=score#): From: 85, To: 115,
- [**200 amino acids**](https://www.uniprot.org/uniprot/?query=length%3A%5B185+TO+215%5D+reviewed%3Ayes&sort=score): From: 185, To: 215 and
- [**300 amino acids**](https://www.uniprot.org/uniprot/?query=length%3A%5B285+TO+315%5D+reviewed%3Ayes&sort=score): From: 285, To: 315.

Select all and click Download. In the Download tab, change Format to Excel and click Go.

Save files as `uniprot_<chain-length>s.xlsx` in `data/alphafold/uniprot_ids/`.

## 2. Retrieve PDB files from AlphaFold 2

Run `scripts/retrieve_alphafold_pdbs.py`. This loops over each of the Uniprot ID files and retrieves the PDB files 
from the AlphaFold 2 [database](https://alphafold.ebi.ac.uk). Files are saved in `/data/alphafold/pdb_files/` as
`AF-<uniprot_ID>-F1-model_v2.pdb`.

## 3. Get amino acid distance distributions

Amino acid distance distributions are computed by running `python run.py af --range=<chain-length> --path=<path-to-pdb-files>`
where `<chain-length>` is one of 100, 200 or 300. This creates Protein Contact Maps from each PDB in the given chain 
length range and calculates the amino acid distances. The amino acid distances will be saved in 
`data/alphafold/lls_<chain-length>.npy`.

## 4. 'Chunk' AlphaFold distances

To bootstrap the distance data, run python `run.py chunk <chain-length> --file=<data/alphafold/lls_<chain-length>.npy>`. 
This will save the "raw" "chunked" data in `data/alphafold/chunk_<chain-length>_raw.csv` and the mean frequencies of 
distances, as well as the confidence level bounds, in `data/alphafold/chunk_<chain-length>_stats.csv`.