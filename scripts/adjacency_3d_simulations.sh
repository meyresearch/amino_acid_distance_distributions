#!/bin/bash
# Bash script to plot different 3D adjacency matrix plots

for i in `seq 0 1 27`
do
	python plot.py adj -t sim -f ../data/simulations/3d/matrices/matrix_300_$i.txt &
done
