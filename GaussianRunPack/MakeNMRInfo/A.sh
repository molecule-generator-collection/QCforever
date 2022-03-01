#!/bin/sh

#$ -S /bin/bash
#$ -cwd
#$ -pe impi 12
##$ -jc pcc-normal
#$ -jc pcc-large
##$ -jc pcc-skl

. /fefs/opt/x86_64/Gaussian/envset.sh
source ~/.bashrc

ulimit -s unlimited

#python /home/sumita/GaussianRun_0.3/ParaChem/main.py  ch2o.sdf > log
./Execute.sh >& log 
