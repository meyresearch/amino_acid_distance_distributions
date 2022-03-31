#!/bin/bash

python plot.py -r=100 rcsb -startd 0.0016 -endd 0.0022 -starte 0 -ende 4 -m mean 
python plot.py -r=200 rcsb -startd 0.0011 -endd 0.0016 -starte 0 -ende 4 -m mean 
python plot.py -r=300 rcsb -startd 0.0006 -endd 0.0012 -starte 2 -ende 7 -m mean 
python plot.py -r=100 alpha -startd 0.0018 -endd 0.0023 -starte 1 -ende 5 -m mean 
python plot.py -r=200 alpha -startd 0.0009 -endd 0.0014 -starte 1 -ende 5 -m mean 
python plot.py -r=300 alpha -startd 0.0009 -endd 0.0014 -starte 4 -ende 8 -m mean 
