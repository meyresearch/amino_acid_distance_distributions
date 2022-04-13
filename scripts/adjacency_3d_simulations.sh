#!/bin/bash

for i in `seq 0 1 27`
do
	python plot.py adj -t sim -f ../data/simulations/3d/matrices/matrix_300_$i.txt &
done
