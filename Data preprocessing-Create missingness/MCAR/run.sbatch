#!/bin/bash
#----------------------------------------------------
# Sample Slurm job script
#   for TACC Stampede2
#
#   *** MPI Job on Normal Queue ***
# 
# Last revised: 4 June 2021
#
# Notes:
#
#   -- If you're running out of memory, try running
#      fewer tasks per node to give each task more memory.
#
#----------------------------------------------------

#SBATCH -J MNAR       # Job name
#SBATCH -o MNAR.o%j   # Name of the output file
#SBATCH -e MNAR.e%j   # Name of the error file
#SBATCH -p cpu           # Queue (partition) name
#SBATCH -N 1             # Total # of nodes 
#SBATCH -n 20             # Total # of MPI tasks
#SBATCH -t 24:30:00      # Run time (hh:mm:ss)



# Launch MPI code
export I_MPI_FABRICS=shm:ofi FI_PROVIDER=verbs;mpirun -np 20 ./main_MPI

# ---------------------------------------------------
