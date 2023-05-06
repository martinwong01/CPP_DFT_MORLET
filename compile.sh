#!/bin/bash


vector=1                       # 1: use avx/avx512/fma, 0: not used
compiler="clang"               # "clang" or "intel" or "gnu" or "nvidia"
maxn=131072                    # set to 2^
maxs=64                        # maximum number of wavelet scales

### you may not need to modify beyond this



avx=0
avx512=0
fma=0
if [ $vector -eq 1 ]; then
    fma=$(grep " fma " /proc/cpuinfo|wc -l)
    if [ $fma -gt 0 ]; then
        avx=$(grep " avx " /proc/cpuinfo|wc -l)
        avx512=$(grep " avx512f " /proc/cpuinfo|wc -l)
    fi
fi
macros="-D AVX512=${avx512} -D AVX=${avx} -D FMA=${fma} -D MAXN=${maxn} -D MAXS=${maxs}"
vector_flags=""
if [ $avx512 -gt 0 ]; then vector_flags="-mavx512f"; fi 
if [ $avx -gt 0 ]; then vector_flags="$vector_flags -mavx"; fi
if [ $fma -gt 0 ]; then vector_flags="$vector_flags -mfma"; fi

if [ $compiler == "gnu" ]; then
    #command="g++ $vector_flags -std=c++17 -fopenmp -Ofast -march=skylake-avx512 $macros"
    command="g++ $vector_flags -std=c++17 -fopenmp -Ofast $macros"
elif [ $compiler == "intel" ]; then
    command="icc $vector_flags -std=c++17 -qopenmp -Ofast -diag-disable=10441 $macros"
elif [ $compiler == "clang" ]; then
    #command="clang++ $vector_flags -std=c++17 -fopenmp=libomp -Ofast -march=skylake-avx512 $macros"
    command="clang++ $vector_flags -std=c++17 -fopenmp=libomp -Ofast $macros"
elif [ $compiler == "nvidia" ]; then
    #command="nvc++ $vector_flags -std=c++17 -fopenmp -Ofast -march=skylake-avx512 $macros --diag_suppress declared_but_not_referenced --diag_suppress set_but_not_used"
    command="nvc++ $vector_flags -std=c++17 -fopenmp -Ofast $macros --diag_suppress declared_but_not_referenced --diag_suppress set_but_not_used"
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
