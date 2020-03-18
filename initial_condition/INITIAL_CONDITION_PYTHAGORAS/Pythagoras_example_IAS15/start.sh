#PBS -N IAS15_test
#PBS -l nodes=1
#PBS -q debug
#PBS -j oe
#PBS -M ozawa_m@eps.s.u-tokyo.ac.jp
#PBS -me
#PBS -r n
   
cd $PBS_O_WORKDIR 
  
#export OMP_NUM_THREADS=4
#export OMP_STACKSIZE=512M
rm -r output
rm energy.dat
#rm nohup.out
#rm collision_after.dat collision_before.dat collision.dat energy.dat *.dat
mkdir output
#mkdir output/snap01
#mkdir output/xyz_snap01
#mkdir output/restart

aprun -n 1 -N 1 ./REBOUND

