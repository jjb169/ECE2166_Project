# ECE2166_Project
## Parallel Computer Architecute (ECE2166) course project exploring shared memory multiprocessors utilized for MST and BFS graph operations 
#### To compile the code locally: g++ -fopenmp -O3 graphXXX.cpp -o graphXXX.X
#### To compile on Pitt CRC or Linux based computing resources with .sh file: ./compileXXX-gnu.sh
####
#### Serial Execution: ./graphXXX "graph_name" DELIMETER
#### Parallel Execution: ./graphXXX "graph_name" DELIMETER #CORES
#### DELIMETER = 0 if using csv file (ie an edge looks like: 1,2,3)
#### DELIMETER = 1 if using txt file with spaces (ie an dedge looks like: 1 2 3)
#### DELIMETER = 2 if using txt file with tabs (ie an dedge looks like: 1 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 2 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; 3)
#### #CORES is the number of cores/threads you would like to use within OpenMP 
####
#### Added in .slurm files for submission to CRC- submit with 'sbatch graphXXX.slurm'
