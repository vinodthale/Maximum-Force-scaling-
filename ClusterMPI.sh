#!/bin/bash
######################################################################
# sample SLURM batch script
# to run do something like "sbatch run.sh"
######################################################################
# set name of job
#SBATCH --job-name=P1
# set working dirrectory
#SBATCH --chdir=/storage/Vat/Vinod/WaterVat/feb2023/ForceAirliquiddrop/Silicone500Cst/V0.10
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err

# set the number of nodes 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --threads-per-core=1   #thread per node use

# Load all the dependant modules
#module purge            #clear all modules from the environment
#module load bullxmpi
#module load intel

qcc -source -D_MPI=1 Bdropimpact.c
mpicc -O2 -Wall -std=c99 -D_MPI=1 -D_FORTIFY_SOURCE=0 _Bdropimpact.c -o Bdropimpact -lm
srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS ./Bdropimpact  
