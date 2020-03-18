#PBS -N 2D_M5_REBOUND_restart
#PBS -l nodes=1
#PBS -q bulk-b
#PBS -j oe
#PBS -M ozawa_m@eps.s.u-tokyo.ac.jp
#PBS -me
#PBS -r n
   
cd $PBS_O_WORKDIR 

tmp_num=`ls -l ../output/snap01/*_Sun | wc -l`
tmp_num=`expr ${tmp_num} - 1`
snap_num=`printf %05d ${tmp_num}`
echo ${snap_num}
  
#export OMP_NUM_THREADS=4
#export OMP_STACKSIZE=512M


aprun -n 1 -N 1 ./restart ${snap_num}

