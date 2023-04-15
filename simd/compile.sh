#!/bin/bash


g++ -mavx -o testsimd_g.exe testsimd.cpp
g++ -o test_g.exe test.cpp
g++ -mavx -o complexsimd_g.exe complexsimd.cpp
g++ -fopenmp -o testpragma_g.exe testpragma.cpp

icc -diag-disable=10441 -o testsimd_i.exe testsimd.cpp
icc -diag-disable=10441 -o test_i.exe test.cpp
icc -diag-disable=10441 -o complexsimd_i.exe complexsimd.cpp
icc -diag-disable=10441 -qopenmp -o testpragma_i.exe testpragma.cpp

#g++ -S test.cpp
#g++ -mavx -S complexsimd.cpp
#g++ -fopenmp -S testpragma.cpp


if [ 1 -eq 0 ]; then
g++ -mavx -S testsimd.cpp
mv testsimd.s testsimd_noopt.s
g++ -Ofast -mavx -S testsimd.cpp
mv testsimd.s testsimd_Ofast.s
g++ -S test.cpp
mv test.s test_noopt.s
g++ -Ofast -S test.cpp
mv test.s test_Ofast.s
g++ -O1 -S test.cpp
mv test.s test_O1.s
g++ -O2 -S test.cpp
mv test.s test_O2.s
g++ -O3 -S test.cpp
mv test.s test_O3.s

fi

