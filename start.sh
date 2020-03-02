#!/bin/bash
   
#export OMP_NUM_THREADS=4
MPIRUN=mpirun

rm -r result
rm log.log
mkdir result




#$MPIRUN
./nbody.out -N 3001  -T 200 -D 10 -n 3005  -s 1e-5 > log.log &
#./nbody.out  -s 1e-9 > log.log &
#T:timeend 


