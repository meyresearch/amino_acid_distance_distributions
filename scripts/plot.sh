#!/bin/bash

# RCSB 100 comparison
python plot.py comp -r=100 --rcsb-startd 0.0015 --rcsb-endd 0.0016 --rcsb-starte 1 --rcsb-ende 3 -q 2 -m mean -startp 4 -endp 50

#-----------------------RCSB THEORY-----------------------
#N: 99.000000 +/- 13.409797 
#H-N/2: 4.499205 +/- 1.259433
#a: 1.532645 +/- 0.785944
#A: 0.001500 +/- 0.000288
#----------------------RCSB POWER LAW----------------------
#gamma: 0.8288237555616077 +/- 0.02177892605304199
#constant: 0.0723936032968393+/- 0.004595302419742699

# RCSB 200 comparison
python plot.py comp -r=200 --rcsb-startd 0.0008 --rcsb-endd 0.0010 --rcsb-starte 2 --rcsb-ende 4 -q 2 -m mean -startp 4 -endp 100

#-----------------------RCSB THEORY-----------------------
#N: 199.000000 +/- 27.828499 
#H-N/2: 5.187378 +/- 1.252461
#a: 2.077775 +/- 0.647283
#A: 0.000800 +/- 0.000130
#----------------------RCSB POWER LAW----------------------
#gamma: 0.8392897880610493 +/- 0.017390772327924024
#constant: 0.06652861846208974+/- 0.003796831547281125

# RCSB 300 comparison
python plot.py comp -r=300 --rcsb-startd 0.0007 --rcsb-endd 0.0009 --rcsb-starte 4 --rcsb-ende 6 -q 2 -m mean -startp 20 -endp 150

#-----------------------RCSB THEORY-----------------------
#N: 300.000000 +/- 15.412474 
#H-N/2: 5.591171 +/- 1.192838
#a: 4.194035 +/- 0.474760
#A: 0.000740 +/- 0.000075
#----------------------RCSB POWER LAW----------------------
#gamma: 1.4455442099013252 +/- 0.02052068532309389
#constant: 0.667044825116215+/- 0.057087884295147495

# 3D simulation 100
python plot.py 3d-sim -r=100 -m mean -q 2 -startd 0.003 -endd 0.0032 -starte 2 -ende 4 -startp 4 -endp 50

#----------------Parameters----------------
#Chain length: 99.99999989179496 +/- 19.023065141692637
#N/2 Harmonic: 4.499195338329424 +/- 1.7627263419474972
#Exponent: 2.000000000035432 +/- 0.7144553859269824
#Dimensionality: 0.0031000000006269424 +/- 0.0006242577889197593

# 3D simulation 200
python plot.py 3d-sim -r=200 -m mean -q 2 -startd 0.001 -endd 0.0013 -starte 0 -ende 2 -startp 4 -endp 100

#----------------Parameters----------------
#Chain length: 199.00000000000003 +/- 56.94275888261524
#N/2 Harmonic: 5.187367517639622 +/- 0.9150485958887088
#Exponent: 0.8589551144924402 +/- 0.949948240647899
#Dimensionality: 0.0011660604925413437 +/- 0.0003817539242172025

# 3D simulation 300
python plot.py 3d-sim -r=300 -m mean -q 2 -startd 0.0006 -endd 0.0008 -starte 0 -ende 2 -startp 4 -endp 150

#----------------Parameters----------------
#Chain length: 299.99999999999994 +/- 88.68405340932338
#N/2 Harmonic: 5.591170588643882 +/- 1.0065599894983752
#Exponent: 0.912573466567533 +/- 1.033185137292476
#Dimensionality: 0.0007000000005955203 +/- 0.00024325295701141397

# 2D simulation 
python plot.py 2d-sim -m mean -q 2 -startd 0.0002 -endd 0.0004 -starte 6 -ende 8 -startp 4 -endp 250

#----------------Parameters----------------
#Chain length: 497.9999999999989 +/- 105.72821437106455
#N/2 Harmonic: 6.096675249432579 +/- 1.1357693277760001
#Exponent: 6.000000000000923 +/- 1.132664567346109
#Dimensionality: 0.00029999999999996295 +/- 5.611169236115122e-05



