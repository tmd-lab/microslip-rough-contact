#!/bin/bash
# Single node 
#SBATCH --nodes=1
#SBATCH --exclusive
### Run with: 
### > export X=5 
### > sbatch --job-name="EPMC"$X run_epmcX.slurm
#SBATCH --export=ALL
# SMP Execution only one task:
#SBATCH --ntasks 1
# Use physical cores not hyper threads:
#SBATCH --threads-per-core=1
# Request multiple cores:
#SBATCH --cpus-per-task=16
# Other details of specific runs:
#SBATCH --job-name=EPMC
# scavenge or commons
#SBATCH --partition=scavenge
#SBATCH --time=2:55:00
### #SBATCH --mail-user=
### #SBATCH --mail-type=ALL
#
echo "I ran on:"
cd $SLURM_SUBMIT_DIR
echo $SLURM_NODELIST
#
grep MemTotal /proc/meminfo
lspci
#
module load MATLAB/2020a
#
matlab -nodisplay -r "main_epmc($X); quit"

