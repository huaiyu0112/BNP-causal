rm(list=ls())  

library(MASS) # for function mvrnorm (drawing from Multivariate Normal Distribution)

rep_seq = c(1:1000) #1:40 # 1:40
n_sample = 500
DATA = NULL

mean_L3_rep = mean_L4_rep = rep(NA, 1000)

# Given Hy = 1
# Coefficients of L4
beta_int_L4_Hy1 = 1.5 
beta_L1_2_L4_Hy1 = 1.5
beta_L2_2_L4_Hy1 = 1.5
beta_L2_3_L4_Hy1 = 1.5
# Coefficients of L3
beta_int_L3_Hy1 = 2
beta_L1_2_L3_Hy1 = 2
beta_L2_2_L3_Hy1 = 2
beta_L2_3_L3_Hy1 = 2
Beta_L4_Hy1 = c(beta_int_L4_Hy1, beta_L1_2_L4_Hy1, beta_L2_2_L4_Hy1, beta_L2_3_L4_Hy1)
Beta_L3_Hy1 = c(beta_int_L3_Hy1, beta_L1_2_L3_Hy1, beta_L2_2_L3_Hy1, beta_L2_3_L3_Hy1)

# Given Hy = 2
# Coefficients of L4
beta_int_L4_Hy2 = -1 
beta_L1_2_L4_Hy2 = -1
beta_L2_2_L4_Hy2 = -1
beta_L2_3_L4_Hy2 = -1
# Coefficients of L3
beta_int_L3_Hy2 = -2
beta_L1_2_L3_Hy2 = -2
beta_L2_2_L3_Hy2 = -2
beta_L2_3_L3_Hy2 = -2
Beta_L4_Hy2 = c(beta_int_L4_Hy2, beta_L1_2_L4_Hy2, beta_L2_2_L4_Hy2, beta_L2_3_L4_Hy2)
Beta_L3_Hy2 = c(beta_int_L3_Hy2, beta_L1_2_L3_Hy2, beta_L2_2_L3_Hy2, beta_L2_3_L3_Hy2)

expit_fn = function(x) { 1/(1+exp(-x)) } # inverse of the logit function

for (i_rep in rep_seq) {

  set.seed(100 + i_rep)
  
  ###############################################
  ########## Simulated data #####################
  ###############################################
  # top level components
  lambda_z = c(0.5, 0.5)
  z = sample.int(2, n_sample, replace=TRUE, prob=lambda_z)
  
  # categorical components 
  Hx = Hy = numeric(n_sample)
  Hx[which(z == 1)] = sample.int(2, sum(z == 1), replace=TRUE, prob=c(0.1, 0.9))
  Hx[which(z == 2)] = sample.int(2, sum(z == 2), replace=TRUE, prob=c(0.9, 0.1))
  # continuous components
  Hy[which(z == 1)] = sample.int(2, sum(z == 1), replace=TRUE, prob=c(0.9, 0.1))
  Hy[which(z == 2)] = sample.int(2, sum(z == 2), replace=TRUE, prob=c(0.1, 0.9))

  # Generate categorical variables
  L1 = L2 =  numeric(n_sample)
  # Given Hx = 1
  L1[which(Hx==1)] = sample.int(2, sum(Hx==1), replace=TRUE, prob=c(0.9, 0.1))
  L2[which(Hx==1)] = sample.int(3, sum(Hx==1), replace=TRUE, prob=c(0.10, 0.10, 0.80))
  # Given Hx = 2
  L1[which(Hx==2)] = sample.int(2, sum(Hx==2), replace=TRUE, prob=c(0.1, 0.9))
  L2[which(Hx==2)] = sample.int(3, sum(Hx==2), replace=TRUE, prob=c(0.40, 0.40, 0.20))

  #Making dummy variables with dummy_cols()
  L1_2 = cbind(L1, L2)
  L1_2_dummies <- fastDummies::dummy_cols(L1_2, select_columns = c("L1", "L2"))
  L1_2_dummies = L1_2_dummies[,-c(1:2)]
  L1_2_dummies = L1_2_dummies[, c("L1_2", "L2_2", "L2_3")]

  # Generate continuous variables
  X_cat_dummies = as.matrix(cbind(1, L1_2_dummies))
  
  L4_L3= array(data= NA, dim = c(n_sample, 2))
  
  for (i in 1: n_sample) {
    
    if (Hy[i] ==1){
      mu_L4_Hy1 = X_cat_dummies %*% Beta_L4_Hy1
      mu_L3_Hy1 = X_cat_dummies %*% Beta_L3_Hy1
      mu_L4_L3_Hy1 = cbind(mu_L4_Hy1, mu_L3_Hy1)
      Sigma_L4_L3_Hy1 = matrix(c(1, 0.3, 0.3, 1), nrow=2, ncol = 2, byrow = T)
      L4_L3[i,] = mvrnorm(n = 1, mu= mu_L4_L3_Hy1[i,], Sigma = Sigma_L4_L3_Hy1)
      
    } else if (Hy[i] == 2){
      mu_L4_Hy2 = X_cat_dummies %*% Beta_L4_Hy2
      mu_L3_Hy2 = X_cat_dummies %*% Beta_L3_Hy2
      mu_L4_L3_Hy2 = cbind(mu_L4_Hy2, mu_L3_Hy2)
      Sigma_L4_L3_Hy2 = matrix(c(1, 0.2, 0.2, 1), nrow=2, ncol = 2, byrow = T)
      L4_L3[i,] = mvrnorm(n = 1, mu= mu_L4_L3_Hy2[i,], Sigma = Sigma_L4_L3_Hy2)  
    } 
  }
  
  L4 = L4_L3[,1]
  L3 = L4_L3[,2]
  L1_4 = cbind(L1, L2, L3, L4)
  
  mean_L3_rep[i_rep] = mean(L3)
  mean_L4_rep[i_rep] = mean(L4)
	DATA[[i_rep]] = L1_4
	
}

mean_L3 = mean(mean_L3_rep)
mean_L4 = mean(mean_L4_rep)

save(DATA, mean_L3,  mean_L4, file = "1_DATA_BD.RData")


rm(list=ls())
load("1_DATA_BD.RData")

# Plots of continuous variables
i_rep =100
L1_4 = DATA[[i_rep]]
L4= L1_4[, 4]
L3 = L1_4[,3]
plot(L3, L4, xlab = "L3", ylab = "L4", pch = 16, col= "grey", xlim = c(-9, 9), ylim = c(-9, 9), cex= 1, cex.lab=1)



