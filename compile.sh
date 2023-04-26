#!/bin/bash


tryavx=1                       # 1: use avx/avx512, 0: not used
compiler="intel"               # "gnu" or "intel"
maxn=131072                    # set to 2^
maxs=64                        # maximum number of wavelet scales

### you may not need to modify beyond this



if [ $tryavx -eq 1 ]; then
    avx=$(grep " avx " /proc/cpuinfo|wc -l)
    avx512=$(grep " avx512f " /proc/cpuinfo|wc -l) 
else
    avx=0
    avx512=0
fi
macros="-D AVX512=${avx512} -D AVX=${avx} -D MAXN=${maxn} -D MAXS=${maxs}"
avxflags=""
if [ $avx512 -gt 0 ]; then avxflags="-mavx512f"; fi 
if [ $avx -gt 0 ]; then avxflags="$avxflags -mavx"; fi

if [ $compiler == "gnu" ]; then
    #command="g++ $avxflags -fopenmp -Ofast -march=skylake-avx512 $macros"
    command="g++ $avxflags -fopenmp -Ofast $macros"
elif [ $compiler == "intel" ]; then
    command="icc $avxflags -diag-disable=10441 -qopenmp -Ofast $macros"
fi


echo $command


$command -o table.exe table.cpp
./table.exe > table.h
$command -o fft_bench.exe fft_bench.cpp
$command -o testdft.exe testdft.cpp
$command -o testmorlet.exe testmorlet.cpp
$command -o wavelet2d.exe wavelet2d.cpp
$command -S wavelet2d.cpp
$command -o wavelet2d_stat.exe wavelet2d_stat.cpp
### put your programs here
