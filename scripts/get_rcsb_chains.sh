for i in `seq 100 100 300`
do
nohup python run.py rcsb -r=$i -p=../data/rcsb/chains_$i.csv &
done
