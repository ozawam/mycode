#!/bin/bash
   
#export OMP_NUM_THREADS=4
MPIRUN=mpirun

rm -r result
rm log.log
mkdir result




#$MPIRUN
./nbody.out -T 100 -D 0.1 -s 1e-5 > log.log &
#T:timeend 


