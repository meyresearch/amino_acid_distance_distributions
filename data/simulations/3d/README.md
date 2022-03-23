# 3D Simulation Data

To get the amino acid distance distributions from the 3D simulated data, run 
`run.py 3d-sim --range=<chain-length>` where `<chain-length>` is one of `100, 200` or `300`.

This will save the "raw" distance data as well as the mean frequencies and 
confidence level bounds in `data/simulations/3d/simulation_<chain-length>_raw.csv`
and `data/simulations/3d/simulation_<chain-length>_stats.csv`, respectively.
