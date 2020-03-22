#!/bin/bash
   
#export OMP_NUM_THREADS=4
MPIRUN=mpirun

rm -r result
rm log.log
mkdir result
mkdir result/original
mkdir result/snap



#$MPIRUN
./nbody.out B_ten_particle -N 12  -T 12000 -D 10 -n 15  -s 1e-3  > log.log &
#./nbody.out  -s 1e-9 > log.log &
#T:timeend 


