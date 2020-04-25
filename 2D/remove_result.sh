#!/bin/sh
start=335
end=406
 
#for i in `seq -f %05g $1 $2`
for i in `seq -f %05g ${start} ${end}`
do
rm -r output/snap01/snap"${i}.dat"
rm -r output/snap01/snap"${i}.dat_Sun"
rm -r output/restart/restart"${i}.bin"
rm -r output/xyz_snap01/snap"${i}.dat"
rm -r output/xyz_snap01/snap"${i}.dat_sp"
done
 

