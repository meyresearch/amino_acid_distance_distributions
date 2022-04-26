# Plotting

This directory contains all the figures used in the main paper as well as in the supplementary
material. All plotting code is available in [`../scripts/`](https://github.com/meyresearch/sequence_distance_distribution/tree/bump/scripts).
In general, you should only need to interact with `plot.py` and use the different command line arguments
to create the different plots. Included in this directory are also the PowerPoint files that can be used to create the final plots.

The different command line options can be viewed by calling `python plot.py --help`.

The bash script `../scripts/plot.sh` can also be run to plot all the plots described below.

## Figure 1

Fig. 1. is available both as a pdf and a svg file.

## Figure 2

### Fig. 2a. 

To plot Fig.2a. run 
```
python plot.py bar
```
This calls the `bars.py` script and creates the bar plots with data from 
[`../data/pdb_statistics/`](https://github.com/meyresearch/sequence_distance_distribution/tree/bump/data/pdb_statistics).

This option does not take any arguments.

### Fig. 2b.

Fig. 2b. is an example of an adjacency matrix from [PDB ID: 3UFG](https://www.rcsb.org/structure/3UFG).
This structure is included in `adjacency_matrices/`. To create this figure, run 
```
python plot.py adj --type=pdb --file=../plots/adjacency_matrices/3ufg.pdb
```
The plot is automatically saved as a jpeg in
`../plots/adjacency_matrices/` as `pdb_matrix.jpg`.

## Figure 3

### Fig. 3a.
To plot Fig. 3a. for the 2D simulation amino acid distance distribution, run 
```
python plot.py 2d-sim --measure=mean --quantile=2 --start-dimensionality=0.0002 --end-dimensionality=0.0004 --start-exponent=6 --end-exponent=8 --start-point=4 --end-point=250
```

The code will also print out the fit parameters on the command line. 

### Fig. 3b.

Similarly, for the 3D simulation plot, run
```
python plot.py 3d-sim --range=300 --measure=mean --quantile=2 --start-dimensionality 0.0006 --end-dimensionality 0.0008 --start-exponent 0 --end-exponent 2 --start-point 4 --end-point 150
```

and similarly for the `100`:
```
python plot.py 3d-sim -r=100 -m mean -q 2 -startd 0.003 -endd 0.0032 -starte 2 -ende 4 -startp 4 -endp 50
```
and `200`:
```
python plot.py 3d-sim -r=200 -m mean -q 2 -startd 0.001 -endd 0.0013 -starte 0 -ende 2 -startp 4 -endp 100

```
 
### Fig. 3c.

Fig. 3c. is an example of an adjacency matrix from a 3D simulation run. This matrix was created from 
`../data/simulations/3d/matrices/matrix_300_14.txt`. To create this figure, run 
```
python plot.py adj --type=sim --file=../data/simulations/3d/matrices/matrix_300_14.txt
```
The plot is automatically
saved as a jpeg in `../plots/adjacency_matrices/` as `sim_matrix.jpg`.

## Figure 4

### Fig. 4a. 

To plot the comparison of all three ranges for RCSB, run
```
python plot.py rcsb
```
and similarly for AlphaFold:
```
python plot.py alpha
```

### Fig. 4b.

To plot the theoretical fit to `300` data, run:
```
python plot.py comp --range=300 --rcsb-startd 0.0007 --rcsb-endd 0.0009 --rcsb-starte 4 --rcsb-ende 6 --quantile 2 --measure-mean --start-point 20 --end-point 150
```
For `100`:
```
python plot.py comp -r=100 --rcsb-startd 0.0015 --rcsb-endd 0.0016 --rcsb-starte 1 --rcsb-ende 3 -q 2 -m mean -startp 4 -endp 50
```
and for `200`:
```
python plot.py comp -r=200 --rcsb-startd 0.0008 --rcsb-endd 0.0010 --rcsb-starte 2 --rcsb-ende 4 -q 2 -m mean -startp 4 -endp 100
```

### Fig. 4c.

To plot all the contour plots for all ranges, run
```
python plot.py cont
```