#!/bin/bash
module purge
module load gcc/9.2.0
g++ -fopenmp -O3 graphBFS_baseline.cpp -o graphBFS_baseline.X
