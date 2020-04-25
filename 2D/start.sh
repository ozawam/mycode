#!/bin/bash

CUDIR=$(pwd)

cd ../
DIR1=$(basename `pwd`)

cd ${CUDIR}
DIR2=$(basename `pwd`)
name=1AU_${DIR1}_${DIR2}


rm -r output
rm log.log
rm collision_after.dat  energy.dat *.dat
mkdir output
mkdir output/snap01
mkdir output/xyz_snap01
mkdir output/restart

#./rebound &
nohup ./REBOUND $name  > log.log  &
