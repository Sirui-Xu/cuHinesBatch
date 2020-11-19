#!/bin/sh
echo "compile serial program..."
g++ -std=c++11 -fopenmp -lm ../src/serial.cc -o ../src/serial
echo "compile parallel program..."
g++ -std=c++11 -fopenmp -lm ../src/parallel.cc -o ../src/parallel