#!/bin/bash
for i in `seq 100 100 300`
do
	nohup python run.py alpha -r=$i -p=../data/alphafold/unique_secondary_structures_$i.csv > nohup_run_not_normed_$i.out &
done
