for i in ../data/rcsb/ids/ids_* 
do
cd $i
nohup python retrieve_rcsb_pdbs.py
cd ../
done
