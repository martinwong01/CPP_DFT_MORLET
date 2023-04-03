#!/bin/bash

optim="-fopenmp -Ofast -fopt-info-vec-all"

g++ $optim -o testdft.exe testdft.cpp
g++ $optim -o testmorlet.exe testmorlet.cpp
g++ $optim -o wavelet2d.exe wavelet2d.cpp
