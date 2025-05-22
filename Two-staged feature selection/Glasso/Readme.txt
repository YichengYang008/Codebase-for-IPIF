In this example for running R package glasso, we need two input files with the exact names:

(1) selected_daty.txt (dataset with 2000 instances and 12 variables)
(2) selected_id.txt (global indices of 12 variables)

The original input dataset has 2000 instances with 1000 variables, and every eight variables are internally dependent. 
We adopted parallel mutual information to extract 12 important features among 1000 variables. Note that the variable indexed by 1000 is the target. 
In this example, we will adopt the graphical lasso to build a graphic structure of these 12 features and further extract 8 important features.

See "Appendix to Ulta Data-Oriented Parallel Fractional Hot-Deck Imputation with Efficient Linearized Variance Estimation" for more details. 