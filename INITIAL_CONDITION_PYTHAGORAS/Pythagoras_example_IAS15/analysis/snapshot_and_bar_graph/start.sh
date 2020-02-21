#PBS -N A_REBOUND
#PBS -l nodes=1
#PBS -q test-b
#PBS -j oe
#PBS -M ozawa_m@eps.s.u-tokyo.ac.jp
#PBS -me
#PBS -r n
   
cd $PBS_O_WORKDIR 
#export OMP_NUM_THREADS=4
#export OMP_STACKSIZE=512M

./remove_dat.sh
# start directory, end directory, f gas, start au
aprun -n 1 -N 1 ./snapshot 1 1 1.0 10.75


