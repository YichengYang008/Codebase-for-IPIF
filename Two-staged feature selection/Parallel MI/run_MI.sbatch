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

#SBATCH -J P-MI          # Job name
#SBATCH -o P-MI.o%j      # Name of the output file
#SBATCH -e P-MI.e%j      # Name of the error file
#SBATCH -p normal        # Queue (partition) name
#SBATCH -N 1             # Total # of nodes 
#SBATCH -n 6             # Total # of MPI tasks
#SBATCH -t 00:30:00      # Run time (hh:mm:ss)

# Load the Intel compiler
module load intel/18.0.2

# Launch MPI code
ibrun ./main_MPI

# ---------------------------------------------------
