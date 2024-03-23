#!/bin/bash


## this program requires 6GB RAM to run, output files ~40GB, and may take 10 min. 
## outputfile1  will be the "left moving"  transform
## outputfile2  will be the "right-moving" transform


export OMP_NUM_THREADS=4
export OMP_STACKSIZE=10000000

./wavelet2d.exe olr_sym_0 outputfile 16071 144 35 33 1 1
#~/SDE/sde-external-9.14.0-2022-10-25-lin/sde64 -icl -- ./wavelet2d.exe olr_sym_0 outputfile 16071 144 35 33 1 1



## some verifications
./wavelet2d_stat.exe
