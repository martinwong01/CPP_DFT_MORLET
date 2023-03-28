#!/bin/bash

optim="fast"

g++ -fopenmp -O$optim -o testdft.exe testdft.cpp
g++ -fopenmp -O$optim -o testmorlet.exe testmorlet.cpp
g++ -fopenmp -O$optim -o wavelet2d.exe wavelet2d.cpp
