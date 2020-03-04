#!/bin/bash
   
#export OMP_NUM_THREADS=4
MPIRUN=mpirun

rm -r result
rm log.log
mkdir result




#$MPIRUN
./nbody.out no_collision_IDA1992 -N 401 -T 2000 -D 1 -n 401  -s 1e-3 > log.log &
#./nbody.out  -s 1e-9 > log.log &
#T:timeend 


