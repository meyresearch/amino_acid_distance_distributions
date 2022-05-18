#!/bin/bash

# RCSB ranges 
python plot.py rcsb

# AlphaFold 2 ranges
python plot.py alpha

# RCSB 100 comparison
python plot.py comp -r=100 --rcsb-startd 0.0015 --rcsb-endd 0.0016 --rcsb-starte 1 --rcsb-ende 3 -q 2 -m mean -startp 10 -endp 50

# RCSB 200 comparison
python plot.py comp -r=200 --rcsb-startd 0.0007 --rcsb-endd 0.0009 --rcsb-starte 2 --rcsb-ende 4 -q 2 -m mean -startp 20 -endp 100

# RCSB 300 comparison
python plot.py comp -r=300 --rcsb-startd 0.0007 --rcsb-endd 0.0009 --rcsb-starte 4 --rcsb-ende 6 -q 2 -m mean -startp 20 -endp 150

# Contour plots
python plot.py cont

# 3D simulation 100
python plot.py 3d-sim -r=100 -m mean -q 2 -startd 0.003 -endd 0.0032 -starte 2 -ende 4 -startp 4 -endp 50

# 3D simulation 200
python plot.py 3d-sim -r=200 -m mean -q 2 -startd 0.001 -endd 0.0013 -starte 0 -ende 2 -startp 4 -endp 100

# 3D simulation 300
python plot.py 3d-sim -r=300 -m mean -q 2 -startd 0.0006 -endd 0.0008 -starte 0 -ende 2 -startp 4 -endp 150

# 2D simulation 
python plot.py 2d-sim -m mean -q 2 -startd 0.0002 -endd 0.0004 -starte 6 -ende 8 -startp 4 -endp 250





