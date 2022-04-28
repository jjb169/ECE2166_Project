###
 # @Author: shen fan
 # @Date: 2022-04-26 00:56:25
 # @LastEditors: shen fan
 # @LastEditTime: 2022-04-26 00:56:25
 # @Description:       
 # @FilePath: /pj/graphBFS_baseline.sh
### 
#!/bin/bash
module purge
module load gcc/9.2.0
g++ -fopenmp -O3 graphBFS_base.cpp -o graphBFS_base.X
