#!/bin/bash

rm -r output
rm nohup.out
rm collision_after.dat collision_before.dat collision.dat energy.dat *.dat
mkdir output
mkdir output/snap01
mkdir output/xyz_snap01
mkdir output/restart

./restart_3D_dt_1e_5_M8_reproduction   &
