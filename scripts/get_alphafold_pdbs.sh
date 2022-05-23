for i in `seq 100 100 300`
do
nohup python retrieve_alphafold_pdbs.py $i &
done
