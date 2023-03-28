#!/bin/sh
#BSUB -J Global
#BSUB -q hpc
#BSUB -R "rusage[mem=100GB]"
#BSUB -M 100GB
#BSUB -B
#BSUB -N
#BSUB -u kha@aqua.dtu.dk
#BSUB -o Output_%J.txt
#BSUB -e Error_%J.txt
#BSUB -W 4:00 
#BSUB -n 1

module load matlab/R2021b

#cd /zhome/09/7/104700/local/cluster_runs/
matlab -nodisplay -r runGlobalCluster -logfile cluster.log
