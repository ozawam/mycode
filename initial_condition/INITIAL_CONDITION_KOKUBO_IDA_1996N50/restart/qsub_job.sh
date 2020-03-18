#!/bin/sh

job_name=$2.sdb
echo ${job_name}

for job_num in `seq 1 $1`
do
  job_id=`qsub -W depend=afterany:${job_name} auto_restart.sh`
  job_name=${job_id}
  echo ${job_num} ${job_name}
done
