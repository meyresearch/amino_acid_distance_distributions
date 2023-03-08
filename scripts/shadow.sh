#!/bin/bash

path_to_pdb=$1
cutoff=$2
shadow=$3
IFS="/" read -ra PATH <<< $path_to_pdb

pdb_file=${PATH[-1]}
filename=${pdb_file//.pdb}

smog_dir="/home/jguven/Software/smog-2.4.5/bin/"
shadow_dir="/home/jguven/projects/amino_acid_distance_distributions/data/alphafold/shadow_maps_s_${shadow}_c_${cutoff}_A"
$smog_dir/smog_adjustPDB -i $path_to_pdb -o $shadow_dir/${filename}_adjusted.pdb -altLocs

$smog_dir/smog2 -i $shadow_dir/${filename}_adjusted.pdb -AA -o $shadow_dir/${filename}_adjusted.top -g $shadow_dir/${filename}_adjusted.gro -s $shadow_dir/${filename}_adjusted.contacts -backup no

scm_dir="/home/jguven/Software/smog-2.4.5/src/tools"
/usr/bin/java -jar $scm_dir/SCM.jar -g $shadow_dir/${filename}_adjusted.gro -t $shadow_dir/${filename}_adjusted.top -o $shadow_dir/${filename}_contacts --default --coarse CA -s ${shadow} -c ${cutoff} --distance -po $shadow_dir/${filename}_scm.pdb

