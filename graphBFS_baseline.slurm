#!/bin/bash
#SBATCH --job-name=graphBFS
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --time=00:05:00
#SBATCH --account=ece2166_2022s

#purge and load needed modules
module purge
module load gcc/9.2.0

# Give details of CPU running program
echo -n "Running on: "
cat /proc/cpuinfo | grep "model name" | head -n 1 | cut -c 13-


#execute in typical fashion
./graphBFS_base.X graphs/mawi_201512020130.mt 1 &> graphBFS_base_130.out

#conclude with stats
crc-job-stats.py
