#!/bin/bash
while read line
do
  col1=`echo ${line} | cut -d ':' -f 1`
  col2=`echo ${line} | cut -d ':' -f 2`
  echo "Column1: ${col1} Column2: ${col2}"
done<./seed.dat

#  echo "Column1: ${col1} Column2: ${col2}"
kill ${col2}

