#!/bin/bash

i=0
thread=$(nproc --all)
cd work
for com in *.com; do
  echo $com
  (g16 $com) &
  i=$((i+1))
  if [[ $i -ge $thread ]];then
    wait
    i=0
  fi
done
cd ..

