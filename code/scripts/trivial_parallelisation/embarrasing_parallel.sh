for i in {1..14} 
do 
cd ids/id_file_$i
nohup python retrieve_pdbs_parallel_$i.py &
cd ../../
done 
