###
 # @Author: shen fan
 # @Date: 2022-04-26 00:53:21
 # @LastEditors: shen fan
 # @LastEditTime: 2022-04-26 00:53:21
 # @Description:       
 # @FilePath: /pj/graphBFS_omp.sh
### 
#!/bin/bash
module purge
module load gcc/9.2.0
g++ -fopenmp -O3 graphBFS_omp.cpp -o graphBFS_omp.X
