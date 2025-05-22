#-----------------------------------------
# Example R code of using graphical lasso
# Developer: Yicheng Yang (Iowa State Univeristy)
# Modified on June 15, 2021
#-----------------------------------------

#Load required packages
library(network)
library(glasso)

#Import input files
i_row = 2000 #number of instances
i_col = 12 #number of variables

#Read selected features by parallel MI. Note that users must define their own directory to input files
selected_daty <- scan('C:/Users/yicheng/Box Sync/IEEE UP-FHDI/Glasso/selected_daty.txt')
selected_daty <- matrix(selected_daty, ncol = i_col, byrow = TRUE)

#Read global ID of selected features 
selected_id <- scan('C:/Users/yicheng/Box Sync/IEEE UP-FHDI/Glasso/selected_id.txt')
selected_id <- matrix(selected_id, ncol = i_col, byrow = TRUE)

#Rename indices of seletced features with global ID
colnames(selected_daty) = as.character(selected_id)  

#--------------
#Perform glasso
#--------------
S <- cov(selected_daty) #covariance matrix of selected features
a <- glasso(S, 0.6) #apply lasso penalty
P <- a$wi #estimated precision matrix

#Visualize the graphic model
A <- ifelse(P!=0 & row(P)!=col(P),1,0)
g <- network(A)
plot(g,label=selected_id,main="Graphical model with rho = 0.6")

