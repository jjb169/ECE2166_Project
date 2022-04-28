#!/bin/bash
module purge
module load gcc/9.2.0
g++ -fopenmp -O3 graphMST_baseline.cpp -o graphMST_baseline.X
