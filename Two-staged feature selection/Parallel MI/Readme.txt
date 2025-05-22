In this example for running parallel MI, we need three input files with the exact names:

(1) input.txt
(2) daty_column_binary.bin (Note that the data distribution in this file must be column-oriented)
(3) run.sbatch

The input dataset has 2000 instances with 1000 variables, and every eight variables are internally dependent. 
Suppose we need to select 12 important variables, including the target variable.

The steps to run parallel MI with the example dataset are summarized as follows:
1. Copy input files to the directory of the parallel MI program. Change the name of the input file with the exact name "input.txt".
   Change the name of the raw data file with the exact name "daty_column_binary.bin". Change the name of job script with the exact
   name "run.sbatch"

2. Compile the main function using the Intel compiler:
   mpicxx -o main_MPI main_MPI.cpp

2. Run the job script:
   sbatch run.sbatch

See "Appendix to Ulta Data-Oriented Parallel Fractional Hot-Deck Imputation with Efficient Linearized Variance Estimation" for more details. 