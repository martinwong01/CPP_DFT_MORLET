#!/bin/bash

optim="-mavx -fopenmp -Ofast"

g++ $optim -o table.exe table.cpp
./table.exe > table.h

g++ $optim -o testdft.exe testdft.cpp
g++ $optim -o testmorlet.exe testmorlet.cpp
g++ $optim -o wavelet2d.exe wavelet2d.cpp


#checkavx=$(grep "avx" /proc/cpuinfo|wc -l)
#checkavx=0
#g++ -Ofast -mavx -o fft_bench.exe fft_bench.cpp -D"AVX=${checkavx}"

