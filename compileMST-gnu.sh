#!/bin/bash
module purge
module load gcc/9.2.0
g++ -fopenmp -O3 graphMST.cpp -o graphMST.X
