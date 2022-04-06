#!/bin/bash

python plot.py -r=100 comp -rcsbd 0.0019 -rcsbe 2 -alphad 0.0021 -alphae 3 -m=mean -q=1
python plot.py -r=200 comp -rcsbd 0.0012 -rcsbe 3 -alphad 0.0011 -alphae 3 -m=mean -q=1
python plot.py -r=300 comp -rcsbd 0.0009 -rcsbe 4 -alphad 0.0010 -alphae 5 -m=mean -q=1
