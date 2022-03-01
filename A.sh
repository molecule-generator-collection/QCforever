#!/bin/sh

#$ -S /bin/bash
#$ -cwd
#$ -pe OpenMP 20
##$ -jc pcc-skl
#$ -jc pcc-large

. /fefs/opt/x86_64/Gaussian/envset.sh
source ~/.bashrc

ulimit -s unlimited
export OMP_NUM_THREADS=$NSLOTS

#conda activate py37
python -V

#python /home/sumita/GaussianRun_2.0/main.py ch2o.sdf  > log
python /home/sumita/GaussianRun_2.0/main.py TENDI-X3LYP-631Gd-TD.chk > TENDI.log
#python /home/sumita/GaussianRun_2.0/main.py CheckMolopt11.sdf  > log
#./test_unix.sh
