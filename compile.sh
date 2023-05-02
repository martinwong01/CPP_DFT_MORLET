#!/bin/bash


tryvector=1                    # 1: use avx/avx512/fma, 0: not used
compiler="intel"               # "clang" or "intel" or "gnu"
maxn=131072                    # set to 2^
maxs=64                        # maximum number of wavelet scales

### you may not need to modify beyond this



avx=0
avx512=0
fma=0
if [ $tryvector -eq 1 ]; then
    fma=$(grep " fma " /proc/cpuinfo|wc -l)
    if [ $fma -gt 0 ]; then
        avx=$(grep " avx " /proc/cpuinfo|wc -l)
        avx512=$(grep " avx512f " /proc/cpuinfo|wc -l)
    fi
fi
macros="-D AVX512=${avx512} -D AVX=${avx} -D FMA=${fma} -D MAXN=${maxn} -D MAXS=${maxs}"
flags=""
if [ $avx512 -gt 0 ]; then flags="-mavx512f"; fi 
if [ $avx -gt 0 ]; then flags="$flags -mavx"; fi
if [ $fma -gt 0 ]; then flags="$flags -mfma"; fi

if [ $compiler == "gnu" ]; then
    #command="g++ $flags -fopenmp -Ofast -march=skylake-avx512 $macros"
    command="g++ $flags -fopenmp -Ofast $macros"
elif [ $compiler == "intel" ]; then
    command="icc $flags -diag-disable=10441 -qopenmp -Ofast $macros"
elif [ $compiler == "clang" ]; then
    command="clang++ $flags -fopenmp=libomp -Ofast $macros"
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
