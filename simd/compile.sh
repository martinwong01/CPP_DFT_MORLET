#!/bin/bash


g++ -mavx -o testsimd.exe testsimd.cpp
g++ -o test.exe test.cpp
g++ -o complexsimd.exe complexsimd.cpp

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
