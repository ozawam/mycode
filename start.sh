#!/bin/bash
   
#export OMP_NUM_THREADS=4
MPIRUN=mpirun

rm -r result
rm log.log
mkdir result
mkdir result/original
mkdir result/snap



#$MPIRUN
./nbody.out cpu_test_Beauge1994 -N 1002  -T 0.2 -D 0.1 -n 1005  -s 1e-3  > log.log &
#./nbody.out  -s 1e-9 > log.log &
#T:timeend 


