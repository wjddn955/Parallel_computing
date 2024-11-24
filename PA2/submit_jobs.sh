mkdir results

numproc=8

for p in 1 2 4 8 16 32
do
  n=`echo "($p-1)/$numproc+1" | bc`
  sge_options="-V -cwd -S /bin/bash -l h_rt=1:00:00 -q leopard-short.q -pe leopard-short $(($numproc*$n)) -N mpi.$p -o ./results -e ./results -v p=$p"

  qsub $sge_options ./run.sh 
done
