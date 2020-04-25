#!/bin/sh

for i in `seq  $1 $2`
do
  qdel "${i}"
done

