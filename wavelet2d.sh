#!/bin/bash


## this program requires 6GB RAM to run, output files ~40GB, and may take more than 20 min. 
## outputfile1  will be the "left moving"  transform
## outputfile2  will be the "right-moving" transform


export OMP_NUM_THREADS=4
./wavelet2d.exe olr_sym_0 outputfile 16071 144 35 33 1 1

