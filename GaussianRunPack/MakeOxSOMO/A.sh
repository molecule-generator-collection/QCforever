#!/bin/sh

#$ -S /bin/bash
#$ -cwd
#$ -pe impi 4
##$ -jc pcc-normal
#$ -jc pcc-large
##$ -jc pcc-skl

. /fefs/opt/x86_64/Gaussian/envset.sh
source ~/.bashrc

ulimit -s unlimited

conda activate py37
python -V

#/home/sumita/GaussianRun_1.3/GaussianRunPack/MakeOxSOMO/Execute.sh >& log 
./Execute.sh >& log 
