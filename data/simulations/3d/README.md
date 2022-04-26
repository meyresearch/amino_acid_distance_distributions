# 3D Simulation Data

The adjacency matrices from 3D simulations are available in `matrices/` and more information is available in 
[this paper](https://doi.org/10.1371/journal.pone.0229230). 

To compute the amino acid distance distributions from these matrices, for `100`, run

```
python run.py 3d-sim --range=100
```
and for `200`:
```
python run.py 3d-sim --range=200
```
and `300`:
```
python run.py 3d-sim --range=300
```

This will create histograms and save the data in `histogram_<chain-length-range>_not_normed.npy`.