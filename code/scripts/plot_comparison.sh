#!/bin/bash

python plot.py -r=100 BOTH --d-BS 0.0019 --e-BS 2 --d-C 0.0021 --e-C 3
python plot.py -r=200 BOTH --d-BS 0.0012 --e-BS 3 --d-C 0.0011 --e-C 3
python plot.py -r=300 BOTH --d-BS 0.0009 --e-BS 4 --d-C 0.0010 --e-C 5
