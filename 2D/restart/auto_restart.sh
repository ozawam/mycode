
#!/bin/bash

name=M5_2d_restart

tmp_num=`ls -l ../output/snap01/*_Sun | wc -l`
tmp_num=`expr ${tmp_num} - 1`
snap_num=`printf %05d ${tmp_num}`
echo ${snap_num}
  
rm log.log
rm collision_after.dat  energy.dat *.dat

nohup ./restart ${snap_num} $name  > log.log  &
