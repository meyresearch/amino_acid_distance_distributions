#!/bin/bash
for i in `seq 100 100 300`
do
	nohup python run.py rcsb -r=$i -p=../data/rcsb/secondary_structures_$i.csv > ../data/rcsb/nohup_run_$i.out &
done
