rm(list=ls())  

library(MASS) # for function mvrnorm (drawing from Multivariate Normal Distribution)

expit_fn = function(x) { 1/(1+exp(-x)) } # inverse of the logit function

# True causal effect 
ATE_true = c(1.5, 0.5)
n_sample = 1000

rep_seq = c(1:400) #1:40 # 1:40
CompleteDATA = IncompleteDATA = NULL

for (i_rep in rep_seq) {
  
  print(i_rep)
  set.seed(100 + i_rep)
  
  # top level components
  lambda_z = c(0.05, 0.4, 0.3, 0.15, 0.1)
  z = sample.int(5, n_sample, replace=TRUE, prob=lambda_z)
  # categorical components 
  s = r = numeric(n_sample)
  s[which(z == 1)] = sample.int(5, sum(z == 1), replace=TRUE, prob=c(0.15, 0.10, 0.35, 0.15, 0.25))
  s[which(z == 2)] = sample.int(5, sum(z == 2), replace=TRUE, prob=c(0.35, 0.20, 0.15, 0.15, 0.15))
  s[which(z == 3)] = sample.int(5, sum(z == 3), replace=TRUE, prob=c(0.25, 0.15, 0.15, 0.35, 0.10))
  s[which(z == 4)] = sample.int(5, sum(z == 4), replace=TRUE, prob=c(0.15, 0.25, 0.15, 0.15, 0.30))
  s[which(z == 5)] = sample.int(5, sum(z == 5), replace=TRUE, prob=c(0.45, 0.15, 0.15, 0.10, 0.15))
  # continuous components
  r[which(z == 1)] = sample.int(3, sum(z == 1), replace=TRUE, prob=c(0.50, 0.20, 0.30))
  r[which(z == 2)] = sample.int(3, sum(z == 2), replace=TRUE, prob=c(0.70, 0.15, 0.15))
  r[which(z == 3)] = sample.int(3, sum(z == 3), replace=TRUE, prob=c(0.40, 0.20, 0.40))
  r[which(z == 4)] = sample.int(3, sum(z == 4), replace=TRUE, prob=c(0.20, 0.50, 0.30))
  r[which(z == 5)] = sample.int(3, sum(z == 5), replace=TRUE, prob=c(0.30, 0.20, 0.50))
  
  # Generate categorical variables
  L1 = L2 = L3 = L4 = L5 = numeric(n_sample)
  # Given s = 1
  L1[which(s==1)] = sample.int(2, sum(s==1), replace=TRUE, prob=c(0.60, 0.40))
  L2[which(s==1)] = sample.int(4, sum(s==1), replace=TRUE, prob=c(0.55, 0.15, 0.20, 0.10))
  L3[which(s==1)] = sample.int(3, sum(s==1), replace=TRUE, prob=c(0.20, 0.10, 0.70))
  L4[which(s==1)] = sample.int(2, sum(s==1), replace=TRUE, prob=c(0.50, 0.50))
  L5[which(s==1)] = sample.int(3, sum(s==1), replace=TRUE, prob=c(0.30, 0.50, 0.20))
  
  # Given s = 2
  L1[which(s==2)] = sample.int(2, sum(s==2), replace=TRUE, prob=c(0.80, 0.20))
  L2[which(s==2)] = sample.int(4, sum(s==2), replace=TRUE, prob=c(0.55, 0.35, 0.15, 0.05))
  L3[which(s==2)] = sample.int(3, sum(s==2), replace=TRUE, prob=c(0.20, 0.20, 0.60))
  L4[which(s==2)] = sample.int(2, sum(s==2), replace=TRUE, prob=c(0.90, 0.10))
  L5[which(s==2)] = sample.int(3, sum(s==2), replace=TRUE, prob=c(0.50, 0.10, 0.40))
  
  # Given s = 3
  L1[which(s==3)] = sample.int(2, sum(s==3), replace=TRUE, prob=c(0.55, 0.45))
  L2[which(s==3)] = sample.int(4, sum(s==3), replace=TRUE, prob=c(0.25, 0.15, 0.05, 0.55))
  L3[which(s==3)] = sample.int(3, sum(s==3), replace=TRUE, prob=c(0.20, 0.30, 0.50))
  L4[which(s==3)] = sample.int(2, sum(s==3), replace=TRUE, prob=c(0.50, 0.50))
  L5[which(s==3)] = sample.int(3, sum(s==3), replace=TRUE, prob=c(0.50, 0.30, 0.20))
  
  # Given s = 4
  L1[which(s==4)] = sample.int(2, sum(s==4), replace=TRUE, prob=c(0.55, 0.45))
  L2[which(s==4)] = sample.int(4, sum(s==4), replace=TRUE, prob=c(0.25, 0.35, 0.05, 0.35))
  L3[which(s==4)] = sample.int(3, sum(s==4), replace=TRUE, prob=c(0.30, 0.40, 0.30))
  L4[which(s==4)] = sample.int(2, sum(s==4), replace=TRUE, prob=c(0.65, 0.35))
  L5[which(s==4)] = sample.int(3, sum(s==4), replace=TRUE, prob=c(0.50, 0.20, 0.30))
  
  # Given s = 5
  L1[which(s==5)] = sample.int(2, sum(s==5), replace=TRUE, prob=c(0.25, 0.75))
  L2[which(s==5)] = sample.int(4, sum(s==5), replace=TRUE, prob=c(0.45, 0.05, 0.15, 0.35))
  L3[which(s==5)] = sample.int(3, sum(s==5), replace=TRUE, prob=c(0.45, 0.20, 0.35))
  L4[which(s==5)] = sample.int(2, sum(s==5), replace=TRUE, prob=c(0.75, 0.25))
  L5[which(s==5)] = sample.int(3, sum(s==5), replace=TRUE, prob=c(0.75, 0.15, 0.10))
  
  #Making dummy variables with dummy_cols()
  L1_5 = cbind(L1, L2, L3, L4, L5)
  L1_5_dummies <- fastDummies::dummy_cols(L1_5, select_columns = c("L1", "L2", "L3", "L4", "L5"))
  L1_5_dummies = L1_5_dummies[,-c(1:5)]
  L1_5_dummies = L1_5_dummies[, c("L1_2", "L2_2", "L2_3", "L2_4", "L3_2", "L3_3", "L4_2", "L5_2", "L5_3")]
  
  # Generate continuous variables
  # Given r = 1
  # Coefficients of Y0
  beta_int_Y0_r1 = 0.8 
  beta_L1_2_Y0_r1 = 0.5
  beta_L2_2_Y0_r1 = 1.2
  beta_L2_3_Y0_r1 = 2.5
  beta_L2_4_Y0_r1 = 0.7
  beta_L3_2_Y0_r1 = 1.3
  beta_L3_3_Y0_r1 = 0.5
  beta_L4_2_Y0_r1 = 0.8
  beta_L5_2_Y0_r1 = 1.2
  beta_L5_3_Y0_r1 = 0.5
  
  # Given r = 1
  # Coefficients of X0
  beta_int_X0_r1 = 0.6 
  beta_L1_2_X0_r1 = 1.5
  beta_L2_2_X0_r1 = 2.2
  beta_L2_3_X0_r1 = -2.5
  beta_L2_4_X0_r1 = 3.7
  beta_L3_2_X0_r1 = -1.3
  beta_L3_3_X0_r1 = -0.5
  beta_L4_2_X0_r1 = 2.8
  beta_L5_2_X0_r1 = 3.2
  beta_L5_3_X0_r1 = 0.7
  
  # Coefficients of L6
  beta_int_L6_r1 = 1.5 
  beta_L1_2_L6_r1 = 0.5
  beta_L2_2_L6_r1 = 1
  beta_L2_3_L6_r1 = 2
  beta_L2_4_L6_r1 = 1.5
  beta_L3_2_L6_r1 = 2
  beta_L3_3_L6_r1 = 0.5
  beta_L4_2_L6_r1 = 0.5
  beta_L5_2_L6_r1 = 0.3
  beta_L5_3_L6_r1 = 2.5
  
  # Coefficients of L7
  beta_int_L7_r1 = -0.3 
  beta_L1_2_L7_r1 = -1.5
  beta_L2_2_L7_r1 = -0.7
  beta_L2_3_L7_r1 = -2.5
  beta_L2_4_L7_r1 = -0.5
  beta_L3_2_L7_r1 = -2
  beta_L3_3_L7_r1 = -2.5
  beta_L4_2_L7_r1 = 2.5
  beta_L5_2_L7_r1 = 1.3
  beta_L5_3_L7_r1 = 1.5
  
  Beta_Y0_r1 = c(beta_int_Y0_r1, beta_L1_2_Y0_r1, beta_L2_2_Y0_r1, beta_L2_3_Y0_r1, beta_L2_4_Y0_r1, beta_L3_2_Y0_r1, beta_L3_3_Y0_r1, beta_L4_2_Y0_r1, beta_L5_2_Y0_r1, beta_L5_3_Y0_r1)
  Beta_X0_r1 = c(beta_int_X0_r1, beta_L1_2_X0_r1, beta_L2_2_X0_r1, beta_L2_3_X0_r1, beta_L2_4_X0_r1, beta_L3_2_X0_r1, beta_L3_3_X0_r1, beta_L4_2_X0_r1, beta_L5_2_X0_r1, beta_L5_3_X0_r1)
  Beta_L6_r1 = c(beta_int_L6_r1, beta_L1_2_L6_r1, beta_L2_2_L6_r1, beta_L2_3_L6_r1, beta_L2_4_L6_r1, beta_L3_2_L6_r1, beta_L3_3_L6_r1, beta_L4_2_L6_r1, beta_L5_2_L6_r1, beta_L5_3_L6_r1)
  Beta_L7_r1 = c(beta_int_L7_r1, beta_L1_2_L7_r1, beta_L2_2_L7_r1, beta_L2_3_L7_r1, beta_L2_4_L7_r1, beta_L3_2_L7_r1, beta_L3_3_L7_r1, beta_L4_2_L7_r1, beta_L5_2_L7_r1, beta_L5_3_L7_r1)
  
  # Given r = 2
  # Coefficients of Y0
  beta_int_Y0_r2 = 0.3 
  beta_L1_2_Y0_r2 = 1.5
  beta_L2_2_Y0_r2 = -1.2
  beta_L2_3_Y0_r2 = -1.5
  beta_L2_4_Y0_r2 = 0.7
  beta_L3_2_Y0_r2 = 1.3
  beta_L3_3_Y0_r2 = 1.5
  beta_L4_2_Y0_r2 = 2.8
  beta_L5_2_Y0_r2 = -1.2
  beta_L5_3_Y0_r2 = 1.5
  
  # Coefficients of X0
  beta_int_X0_r2 = 1.3 
  beta_L1_2_X0_r2 = -1.5
  beta_L2_2_X0_r2 = 1.2
  beta_L2_3_X0_r2 = -2.5
  beta_L2_4_X0_r2 = 0.8
  beta_L3_2_X0_r2 = 1.7
  beta_L3_3_X0_r2 = -1.5
  beta_L4_2_X0_r2 = 3.8
  beta_L5_2_X0_r2 = 1.2
  beta_L5_3_X0_r2 = 1.9
  
  # Coefficients of L6
  beta_int_L6_r2 = -2.5 
  beta_L1_2_L6_r2 = -0.1
  beta_L2_2_L6_r2 = -3.1
  beta_L2_3_L6_r2 = -2.2
  beta_L2_4_L6_r2 = 1.5
  beta_L3_2_L6_r2 = -2.0
  beta_L3_3_L6_r2 = -0.5
  beta_L4_2_L6_r2 = -2.5
  beta_L5_2_L6_r2 = -1.3
  beta_L5_3_L6_r2 = -0.5
  # Coefficients of L7
  beta_int_L7_r2 = 1.3 
  beta_L1_2_L7_r2 = 0.5
  beta_L2_2_L7_r2 = -0.7
  beta_L2_3_L7_r2 = 1.5
  beta_L2_4_L7_r2 = 0.5
  beta_L3_2_L7_r2 = -1.5
  beta_L3_3_L7_r2 = 1.5
  beta_L4_2_L7_r2 = 0.5
  beta_L5_2_L7_r2 = -1.3
  beta_L5_3_L7_r2 = 0.3
  
  Beta_Y0_r2 = c(beta_int_Y0_r2, beta_L1_2_Y0_r2, beta_L2_2_Y0_r2, beta_L2_3_Y0_r2, beta_L2_4_Y0_r2, beta_L3_2_Y0_r2, beta_L3_3_Y0_r2, beta_L4_2_Y0_r2, beta_L5_2_Y0_r2, beta_L5_3_Y0_r2)
  Beta_X0_r2 = c(beta_int_X0_r2, beta_L1_2_X0_r2, beta_L2_2_X0_r2, beta_L2_3_X0_r2, beta_L2_4_X0_r2, beta_L3_2_X0_r2, beta_L3_3_X0_r2, beta_L4_2_X0_r2, beta_L5_2_X0_r2, beta_L5_3_X0_r2)
  Beta_L6_r2 = c(beta_int_L6_r2, beta_L1_2_L6_r2, beta_L2_2_L6_r2, beta_L2_3_L6_r2, beta_L2_4_L6_r2, beta_L3_2_L6_r2, beta_L3_3_L6_r2, beta_L4_2_L6_r2, beta_L5_2_L6_r2, beta_L5_3_L6_r2)
  Beta_L7_r2 = c(beta_int_L7_r2, beta_L1_2_L7_r2, beta_L2_2_L7_r2, beta_L2_3_L7_r2, beta_L2_4_L7_r2, beta_L3_2_L7_r2, beta_L3_3_L7_r2, beta_L4_2_L7_r2, beta_L5_2_L7_r2, beta_L5_3_L7_r2)
  
  # Given r = 3
  # Coefficients of Y0
  beta_int_Y0_r3 = -1 
  beta_L1_2_Y0_r3 = -1
  beta_L2_2_Y0_r3 = -1
  beta_L2_3_Y0_r3 = -1
  beta_L2_4_Y0_r3 = -1
  beta_L3_2_Y0_r3 = -1
  beta_L3_3_Y0_r3 = -1
  beta_L4_2_Y0_r3 = -1
  beta_L5_2_Y0_r3 = -1
  beta_L5_3_Y0_r3 = -1
  # Coefficients of X0
  beta_int_X0_r3 = 1.2 
  beta_L1_2_X0_r3 = 1.2
  beta_L2_2_X0_r3 = -1
  beta_L2_3_X0_r3 = -1
  beta_L2_4_X0_r3 = -1
  beta_L3_2_X0_r3 = -1
  beta_L3_3_X0_r3 = -1
  beta_L4_2_X0_r3 = -1
  beta_L5_2_X0_r3 = 1.2
  beta_L5_3_X0_r3 = 1.2
  # Coefficients of L6
  beta_int_L6_r3 = 1 
  beta_L1_2_L6_r3 = 1
  beta_L2_2_L6_r3 = 1
  beta_L2_3_L6_r3 = 1
  beta_L2_4_L6_r3 = 1
  beta_L3_2_L6_r3 = 1
  beta_L3_3_L6_r3 = 1
  beta_L4_2_L6_r3 = 1
  beta_L5_2_L6_r3 = 1
  beta_L5_3_L6_r3 = 1
  # Coefficients of L7
  beta_int_L7_r3 = -1 
  beta_L1_2_L7_r3 = -1
  beta_L2_2_L7_r3 = -1
  beta_L2_3_L7_r3 = -1
  beta_L2_4_L7_r3 = -1
  beta_L3_2_L7_r3 = -1
  beta_L3_3_L7_r3 = -1
  beta_L4_2_L7_r3 = -1
  beta_L5_2_L7_r3 = -1
  beta_L5_3_L7_r3 = -1
  
  Beta_Y0_r3 = c(beta_int_Y0_r3, beta_L1_2_Y0_r3, beta_L2_2_Y0_r3, beta_L2_3_Y0_r3, beta_L2_4_Y0_r3, beta_L3_2_Y0_r3, beta_L3_3_Y0_r3, beta_L4_2_Y0_r3, beta_L5_2_Y0_r3, beta_L5_3_Y0_r3)
  Beta_X0_r3 = c(beta_int_X0_r3, beta_L1_2_X0_r3, beta_L2_2_X0_r3, beta_L2_3_X0_r3, beta_L2_4_X0_r3, beta_L3_2_X0_r3, beta_L3_3_X0_r3, beta_L4_2_X0_r3, beta_L5_2_X0_r3, beta_L5_3_X0_r3)
  Beta_L6_r3 = c(beta_int_L6_r3, beta_L1_2_L6_r3, beta_L2_2_L6_r3, beta_L2_3_L6_r3, beta_L2_4_L6_r3, beta_L3_2_L6_r3, beta_L3_3_L6_r3, beta_L4_2_L6_r3, beta_L5_2_L6_r3, beta_L5_3_L6_r3)
  Beta_L7_r3 = c(beta_int_L7_r3, beta_L1_2_L7_r3, beta_L2_2_L7_r3, beta_L2_3_L7_r3, beta_L2_4_L7_r3, beta_L3_2_L7_r3, beta_L3_3_L7_r3, beta_L4_2_L7_r3, beta_L5_2_L7_r3, beta_L5_3_L7_r3)
  
  X_cat_dummies = as.matrix(cbind(1, L1_5_dummies))
  
  Y0_X0_L6_7= array(data= NA, dim = c(n_sample, 4))
  for (i in 1: n_sample) {
    
    if (r[i] ==1){
      
      mu_Y0_r1 = X_cat_dummies %*% Beta_Y0_r1
      mu_X0_r1 = X_cat_dummies %*% Beta_X0_r1
      mu_L6_r1 = X_cat_dummies %*% Beta_L6_r1
      mu_L7_r1 = X_cat_dummies %*% Beta_L7_r1
      mu_Y0_X0_L6_7_r1 = cbind(mu_Y0_r1, mu_X0_r1, mu_L6_r1, mu_L6_r1)
      
      Sigma_Y0_X0_L6_7_r1 = matrix(c(1, 0.8, 0.7, 0.6, 
                                     0.8, 1, 0.8, 0.5, 
                                     0.7, 0.8, 1, 0.6,
                                     0.6, 0.5, 0.6, 1),
                                   nrow=4, ncol = 4, byrow = T)
      
      Y0_X0_L6_7[i,] = mvrnorm(n = 1, mu= mu_Y0_X0_L6_7_r1[i,], Sigma = Sigma_Y0_X0_L6_7_r1)
      
    } else if (r[i] == 2){
      
      mu_Y0_r2 = X_cat_dummies %*% Beta_Y0_r2
      mu_X0_r2 = X_cat_dummies %*% Beta_X0_r2
      mu_L6_r2 = X_cat_dummies %*% Beta_L6_r2
      mu_L7_r2 = X_cat_dummies %*% Beta_L7_r2
      mu_Y0_X0_L6_7_r2 = cbind(mu_Y0_r2, mu_X0_r2, mu_L6_r2, mu_L7_r2)
      Sigma_Y0_X0_L6_7_r2 = matrix(c(1, 0.7, 0.5, 0.5, 
                                     0.7, 1, 0.6, 0.8, 
                                     0.5, 0.6, 1, 0.5, 
                                     0.5, 0.8, 0.5, 1),
                                   nrow=4, ncol = 4, byrow = T)
      Y0_X0_L6_7[i,] = mvrnorm(n = 1, mu= mu_Y0_X0_L6_7_r2[i,], Sigma = Sigma_Y0_X0_L6_7_r2)  
    } else if (r[i] == 3){
      
      mu_Y0_r3 = X_cat_dummies %*% Beta_Y0_r3
      mu_X0_r3 = X_cat_dummies %*% Beta_X0_r3
      mu_L6_r3 = X_cat_dummies %*% Beta_L6_r3
      mu_L7_r3 = X_cat_dummies %*% Beta_L7_r3
      mu_Y0_X0_L6_7_r3 = cbind(mu_Y0_r3, mu_X0_r3, mu_L6_r3, mu_L7_r3)
      Sigma_Y0_X0_L6_7_r3 = matrix(c(1, 0.5, 0.1,0.3, 
                                     0.5, 1, 0, 0.5, 
                                     0.1, 0, 1, 0.5, 
                                     0.3, 0.5, 0.5, 1),
                                   nrow=4, ncol = 4, byrow = T)
      Y0_X0_L6_7[i,] = mvrnorm(n = 1, mu= mu_Y0_X0_L6_7_r3[i,], Sigma = Sigma_Y0_X0_L6_7_r3)
    }
  }
  
  Y0 = Y0_X0_L6_7[,1]
  X0 = Y0_X0_L6_7[,2]
  L6 = Y0_X0_L6_7[,3]
  L7 = Y0_X0_L6_7[,4]
  L6_7 = cbind(L6, L7)
  L1_7 = cbind(L1, L2, L3, L4, L5, L6, L7)
  
  # Generate treatment indicator A
  expit_fn = function(x) { 1/(1+exp(-x)) } # inverse of the logit function
  prob_temp = expit_fn(-0.3 *(L1_7[,6]+L1_7[,7]) +3*ifelse(L5==1, 1, 0)+ 1.5*ifelse(L5==3, 1, 0)- 5*ifelse(L3==2, 1, 0) - ifelse(L3==3, 1, 0))
  A = TrueA = rbinom(n=n_sample, size=1, prob=prob_temp)
  
  # Generate observed potential outcome
  Y1 = Y0 + ATE_true[1]
  X1 = X0 + ATE_true[2]
  Y_obs = A*Y1+(1-A)*Y0
  X_obs = A*X1+(1-A)*X0
  
  # MAR Missness can be expressed in terms of L5 and L7
  Outcome = Outcome_missing = cbind(Y_obs, X_obs)
  L1_7_missing= L1_7
  
  Y_missing_ind = rbinom(n_sample, 1, expit_fn(0.2 + 0.5*ifelse(L1_7[,5]==1, 1, 0) - 0.3*L1_7[,7]))
  X_missing_ind = rbinom(n_sample, 1, expit_fn(0.6 + 1.2*ifelse(L1_7[,5]==1, 1, 0) - 0.5*L1_7[,7]))
  L1_missing_ind = rbinom(n_sample, 1, expit_fn(-2.5 - 1.5*ifelse(L1_7[,5]==1, 1, 0) + 0.5*L1_7[,7]))
  L6_missing_ind = rbinom(n_sample, 1, expit_fn(-1.5 + 0.1*ifelse(L1_7[,5]==1, 1, 0) - 0.1*L1_7[,7]))
  
  Outcome_missing[Y_missing_ind==1 , 1] = NA
  Outcome_missing[X_missing_ind==1 , 2] = NA
  L1_7_missing[L1_missing_ind==1, 1]= NA
  L1_7_missing[L6_missing_ind==1, 6]= NA
  
  CompleteDATA[[i_rep]] = list(Y_obs=Y_obs, X_obs = X_obs, TrueA=TrueA, L1_7=L1_7)
  IncompleteDATA[[i_rep]] = list(Outcome=Outcome_missing, TrueA=TrueA, L1_7=L1_7_missing)
  
}

save(CompleteDATA, IncompleteDATA, rep_seq, ATE_true, file = "1_DATA.RData")
