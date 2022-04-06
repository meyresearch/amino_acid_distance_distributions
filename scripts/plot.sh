#!/bin/bash

#python plot.py -r=100 rcsb -startd 0.0060 -endd 0.0065 -stepd 0.0001 -starte 0 -ende 4 -startp 4 -endp 75
#python plot.py -r=200 rcsb -startd 0.0060 -endd 0.0063 -starte 2 -ende 4 -startp 4 -endp 100 
#python plot.py -r=300 rcsb -startd 0.0045 -endd 0.0051 -stepd 0.0001 -starte 4 -ende 7 -startp 10 -endp 150

# python plot.py -r=100 alpha -startd 0.0060 -endd 0.0065 -starte 0 -ende 6 -startp 4 -endp 75 
# python plot.py -r=200 alpha -startd 0.0060 -endd 0.0063 -starte 2 -ende 4 -startp 4 -endp 100
# python plot.py -r=300 alpha -startd 0.013 -endd 0.015 -stepd 0.001 -starte 4 -ende 6 -startp 10 -endp 100

python plot.py -r=100 comp --rcsb-startd 0.0060 --rcsb-endd 0.0065 --rcsb-starte 0 --rcsb-ende 6 --alpha-startd 0.0060 --alpha-endd 0.0065 --alpha-starte 0 --alpha-ende 6 -startp 5 -endp 100 

python plot.py -r=200 comp --rcsb-startd 0.0070 --rcsb-endd 0.0071 --rcsb-starte 0 --rcsb-ende 3 --alpha-startd 0.0060 --alpha-endd 0.0063 --alpha-starte 2 --alpha-ende 4 -startp 1 -endp 200 

python plot.py -r=300 comp --rcsb-startd 0.0060 --rcsb-endd 0.0066 --rcsb-starte 4 --rcsb-ende 7 --alpha-startd 0.0056 --alpha-endd 0.0060 --alpha-starte 4 --alpha-ende 7 -startp 10 -endp 150


python plot.py -r=200 rcsb -startd 0.0066 -endd 0.0070 -starte 0 -ende 4 -startp 5 -endp 200
