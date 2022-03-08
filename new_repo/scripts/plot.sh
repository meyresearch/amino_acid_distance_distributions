#!/bin/bash

python plot.py 100 BS --d-begin 0.0016 --d-end 0.0022 --e-begin 0 --e-end 4
python plot.py 200 BS --d-begin 0.0011 --d-end 0.0016 --e-begin 0 --e-end 4
python plot.py 300 BS --d-begin 0.0006 --d-end 0.0012 --e-begin 2 --e-end 7
python plot.py 100 C --d-begin 0.0018 --d-end 0.0023 --e-begin 1 --e-end 5
python plot.py 200 C --d-begin 0.0009 --d-end 0.0014 --e-begin 1 --e-end 5
python plot.py 300 C --d-begin 0.0009 --d-end 0.0014 --e-begin 4 --e-end 8
 
