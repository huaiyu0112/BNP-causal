rm(list=ls())  

library(MASS) # for function mvrnorm (drawing from Multivariate Normal Distribution)

# True causal effect 
ATE_true = 1.5

rep_seq = c(1:400) #1:40 # 1:40
n_sample = 500
DATA = NULL

expit_fn = function(x) { 1/(1+exp(-x)) } # inverse of the logit function

for (i_rep in rep_seq) {
  
  # i_rep = 1
  set.seed(100 + i_rep)
  
  ###############################################
  ########## Simulated data #####################
  ###############################################
  # top level components
  lambda_z = c(0.05, 0.4, 0.3, 0.15, 0.1)
  z = sample.int(5, n_sample, replace=TRUE, prob=lambda_z)
  # categorical components 
  Hx = Hy = numeric(n_sample)
  Hx[which(z == 1)] = sample.int(5, sum(z == 1), replace=TRUE, prob=c(0.15, 0.10, 0.35, 0.15, 0.25))
  Hx[which(z == 2)] = sample.int(5, sum(z == 2), replace=TRUE, prob=c(0.35, 0.20, 0.15, 0.15, 0.15))
  Hx[which(z == 3)] = sample.int(5, sum(z == 3), replace=TRUE, prob=c(0.25, 0.15, 0.15, 0.35, 0.10))
  Hx[which(z == 4)] = sample.int(5, sum(z == 4), replace=TRUE, prob=c(0.15, 0.25, 0.15, 0.15, 0.30))
  Hx[which(z == 5)] = sample.int(5, sum(z == 5), replace=TRUE, prob=c(0.45, 0.15, 0.15, 0.10, 0.15))
  # continuous components
  Hy[which(z == 1)] = sample.int(3, sum(z == 1), replace=TRUE, prob=c(0.50, 0.20, 0.30))
  Hy[which(z == 2)] = sample.int(3, sum(z == 2), replace=TRUE, prob=c(0.70, 0.15, 0.15))
  Hy[which(z == 3)] = sample.int(3, sum(z == 3), replace=TRUE, prob=c(0.40, 0.20, 0.40))
  Hy[which(z == 4)] = sample.int(3, sum(z == 4), replace=TRUE, prob=c(0.20, 0.50, 0.30))
  Hy[which(z == 5)] = sample.int(3, sum(z == 5), replace=TRUE, prob=c(0.30, 0.20, 0.50))
  
  # Generate categorical variables
  L1 = L2 = L3 = L4 = L5 = numeric(n_sample)
  # Given Hx = 1
  L1[which(Hx==1)] = sample.int(2, sum(Hx==1), replace=TRUE, prob=c(0.60, 0.40))
  L2[which(Hx==1)] = sample.int(4, sum(Hx==1), replace=TRUE, prob=c(0.55, 0.15, 0.20, 0.10))
  L3[which(Hx==1)] = sample.int(3, sum(Hx==1), replace=TRUE, prob=c(0.20, 0.10, 0.70))
  L4[which(Hx==1)] = sample.int(2, sum(Hx==1), replace=TRUE, prob=c(0.50, 0.50))
  L5[which(Hx==1)] = sample.int(3, sum(Hx==1), replace=TRUE, prob=c(0.30, 0.50, 0.20))
  
  # Given Hx = 2
  L1[which(Hx==2)] = sample.int(2, sum(Hx==2), replace=TRUE, prob=c(0.80, 0.20))
  L2[which(Hx==2)] = sample.int(4, sum(Hx==2), replace=TRUE, prob=c(0.55, 0.35, 0.15, 0.05))
  L3[which(Hx==2)] = sample.int(3, sum(Hx==2), replace=TRUE, prob=c(0.20, 0.20, 0.60))
  L4[which(Hx==2)] = sample.int(2, sum(Hx==2), replace=TRUE, prob=c(0.90, 0.10))
  L5[which(Hx==2)] = sample.int(3, sum(Hx==2), replace=TRUE, prob=c(0.50, 0.10, 0.40))
  
  # Given Hx = 3
  L1[which(Hx==3)] = sample.int(2, sum(Hx==3), replace=TRUE, prob=c(0.55, 0.45))
  L2[which(Hx==3)] = sample.int(4, sum(Hx==3), replace=TRUE, prob=c(0.25, 0.15, 0.05, 0.55))
  L3[which(Hx==3)] = sample.int(3, sum(Hx==3), replace=TRUE, prob=c(0.20, 0.30, 0.50))
  L4[which(Hx==3)] = sample.int(2, sum(Hx==3), replace=TRUE, prob=c(0.50, 0.50))
  L5[which(Hx==3)] = sample.int(3, sum(Hx==3), replace=TRUE, prob=c(0.50, 0.30, 0.20))
  
  # Given Hx = 4
  L1[which(Hx==4)] = sample.int(2, sum(Hx==4), replace=TRUE, prob=c(0.55, 0.45))
  L2[which(Hx==4)] = sample.int(4, sum(Hx==4), replace=TRUE, prob=c(0.25, 0.35, 0.05, 0.35))
  L3[which(Hx==4)] = sample.int(3, sum(Hx==4), replace=TRUE, prob=c(0.30, 0.40, 0.30))
  L4[which(Hx==4)] = sample.int(2, sum(Hx==4), replace=TRUE, prob=c(0.65, 0.35))
  L5[which(Hx==4)] = sample.int(3, sum(Hx==4), replace=TRUE, prob=c(0.50, 0.20, 0.30))
  
  # Given Hx = 5
  L1[which(Hx==5)] = sample.int(2, sum(Hx==5), replace=TRUE, prob=c(0.25, 0.75))
  L2[which(Hx==5)] = sample.int(4, sum(Hx==5), replace=TRUE, prob=c(0.45, 0.05, 0.15, 0.35))
  L3[which(Hx==5)] = sample.int(3, sum(Hx==5), replace=TRUE, prob=c(0.45, 0.20, 0.35))
  L4[which(Hx==5)] = sample.int(2, sum(Hx==5), replace=TRUE, prob=c(0.75, 0.25))
  L5[which(Hx==5)] = sample.int(3, sum(Hx==5), replace=TRUE, prob=c(0.75, 0.15, 0.10))
  
  #Making dummy variables with dummy_cols()
  L1_5 = cbind(L1, L2, L3, L4, L5)
  L1_5_dummies <- fastDummies::dummy_cols(L1_5, select_columns = c("L1", "L2", "L3", "L4", "L5"))
  L1_5_dummies = L1_5_dummies[,-c(1:5)]
  L1_5_dummies = L1_5_dummies[, c("L1_2", "L2_2", "L2_3", "L2_4", "L3_2", "L3_3", "L4_2", "L5_2", "L5_3")]
  
  # Generate continuous variables
  # Given Hy = 1
  # Coefficients of Y0
  beta_int_Y0_Hy1 = 0.8 
  beta_L1_2_Y0_Hy1 = 0.5
  beta_L2_2_Y0_Hy1 = 1.2
  beta_L2_3_Y0_Hy1 = 2.5
  beta_L2_4_Y0_Hy1 = 0.7
  beta_L3_2_Y0_Hy1 = 1.3
  beta_L3_3_Y0_Hy1 = 0.5
  beta_L4_2_Y0_Hy1 = 0.8
  beta_L5_2_Y0_Hy1 = 1.2
  beta_L5_3_Y0_Hy1 = 0.5
  # Coefficients of L6
  beta_int_L6_Hy1 = 1.5 
  beta_L1_2_L6_Hy1 = 0.5
  beta_L2_2_L6_Hy1 = 1
  beta_L2_3_L6_Hy1 = 2
  beta_L2_4_L6_Hy1 = 1.5
  beta_L3_2_L6_Hy1 = 2
  beta_L3_3_L6_Hy1 = 0.5
  beta_L4_2_L6_Hy1 = 0.5
  beta_L5_2_L6_Hy1 = 0.3
  beta_L5_3_L6_Hy1 = 2.5
  # Coefficients of L7
  beta_int_L7_Hy1 = -0.3 
  beta_L1_2_L7_Hy1 = -1.5
  beta_L2_2_L7_Hy1 = -0.7
  beta_L2_3_L7_Hy1 = -2.5
  beta_L2_4_L7_Hy1 = -0.5
  beta_L3_2_L7_Hy1 = -2
  beta_L3_3_L7_Hy1 = -2.5
  beta_L4_2_L7_Hy1 = 2.5
  beta_L5_2_L7_Hy1 = 1.3
  beta_L5_3_L7_Hy1 = 1.5
  
  Beta_Y0_Hy1 = c(beta_int_Y0_Hy1, beta_L1_2_Y0_Hy1, beta_L2_2_Y0_Hy1, beta_L2_3_Y0_Hy1, beta_L2_4_Y0_Hy1, beta_L3_2_Y0_Hy1, beta_L3_3_Y0_Hy1, beta_L4_2_Y0_Hy1, beta_L5_2_Y0_Hy1, beta_L5_3_Y0_Hy1)
  Beta_L6_Hy1 = c(beta_int_L6_Hy1, beta_L1_2_L6_Hy1, beta_L2_2_L6_Hy1, beta_L2_3_L6_Hy1, beta_L2_4_L6_Hy1, beta_L3_2_L6_Hy1, beta_L3_3_L6_Hy1, beta_L4_2_L6_Hy1, beta_L5_2_L6_Hy1, beta_L5_3_L6_Hy1)
  Beta_L7_Hy1 = c(beta_int_L7_Hy1, beta_L1_2_L7_Hy1, beta_L2_2_L7_Hy1, beta_L2_3_L7_Hy1, beta_L2_4_L7_Hy1, beta_L3_2_L7_Hy1, beta_L3_3_L7_Hy1, beta_L4_2_L7_Hy1, beta_L5_2_L7_Hy1, beta_L5_3_L7_Hy1)
  
  # Given Hy = 2
  # Coefficients of Y0
  beta_int_Y0_Hy2 = 0.3 
  beta_L1_2_Y0_Hy2 = 1.5
  beta_L2_2_Y0_Hy2 = -1.2
  beta_L2_3_Y0_Hy2 = -1.5
  beta_L2_4_Y0_Hy2 = 0.7
  beta_L3_2_Y0_Hy2 = 1.3
  beta_L3_3_Y0_Hy2 = 1.5
  beta_L4_2_Y0_Hy2 = 2.8
  beta_L5_2_Y0_Hy2 = -1.2
  beta_L5_3_Y0_Hy2 = 1.5
  # Coefficients of L6
  beta_int_L6_Hy2 = -2.5 
  beta_L1_2_L6_Hy2 = -0.1
  beta_L2_2_L6_Hy2 = -3.1
  beta_L2_3_L6_Hy2 = -2.2
  beta_L2_4_L6_Hy2 = 1.5
  beta_L3_2_L6_Hy2 = -2.0
  beta_L3_3_L6_Hy2 = -0.5
  beta_L4_2_L6_Hy2 = -2.5
  beta_L5_2_L6_Hy2 = -1.3
  beta_L5_3_L6_Hy2 = -0.5
  # Coefficients of L7
  beta_int_L7_Hy2 = 1.3 
  beta_L1_2_L7_Hy2 = 0.5
  beta_L2_2_L7_Hy2 = -0.7
  beta_L2_3_L7_Hy2 = 1.5
  beta_L2_4_L7_Hy2 = 0.5
  beta_L3_2_L7_Hy2 = -1.5
  beta_L3_3_L7_Hy2 = 1.5
  beta_L4_2_L7_Hy2 = 0.5
  beta_L5_2_L7_Hy2 = -1.3
  beta_L5_3_L7_Hy2 = 0.3
  
  Beta_Y0_Hy2 = c(beta_int_Y0_Hy2, beta_L1_2_Y0_Hy2, beta_L2_2_Y0_Hy2, beta_L2_3_Y0_Hy2, beta_L2_4_Y0_Hy2, beta_L3_2_Y0_Hy2, beta_L3_3_Y0_Hy2, beta_L4_2_Y0_Hy2, beta_L5_2_Y0_Hy2, beta_L5_3_Y0_Hy2)
  Beta_L6_Hy2 = c(beta_int_L6_Hy2, beta_L1_2_L6_Hy2, beta_L2_2_L6_Hy2, beta_L2_3_L6_Hy2, beta_L2_4_L6_Hy2, beta_L3_2_L6_Hy2, beta_L3_3_L6_Hy2, beta_L4_2_L6_Hy2, beta_L5_2_L6_Hy2, beta_L5_3_L6_Hy2)
  Beta_L7_Hy2 = c(beta_int_L7_Hy2, beta_L1_2_L7_Hy2, beta_L2_2_L7_Hy2, beta_L2_3_L7_Hy2, beta_L2_4_L7_Hy2, beta_L3_2_L7_Hy2, beta_L3_3_L7_Hy2, beta_L4_2_L7_Hy2, beta_L5_2_L7_Hy2, beta_L5_3_L7_Hy2)
  
  # Given Hy = 3
  # Coefficients of Y0
  beta_int_Y0_Hy3 = -1 
  beta_L1_2_Y0_Hy3 = -1
  beta_L2_2_Y0_Hy3 = -1
  beta_L2_3_Y0_Hy3 = -1
  beta_L2_4_Y0_Hy3 = -1
  beta_L3_2_Y0_Hy3 = -1
  beta_L3_3_Y0_Hy3 = -1
  beta_L4_2_Y0_Hy3 = -1
  beta_L5_2_Y0_Hy3 = -1
  beta_L5_3_Y0_Hy3 = -1
  # Coefficients of L6
  beta_int_L6_Hy3 = 1 
  beta_L1_2_L6_Hy3 = 1
  beta_L2_2_L6_Hy3 = 1
  beta_L2_3_L6_Hy3 = 1
  beta_L2_4_L6_Hy3 = 1
  beta_L3_2_L6_Hy3 = 1
  beta_L3_3_L6_Hy3 = 1
  beta_L4_2_L6_Hy3 = 1
  beta_L5_2_L6_Hy3 = 1
  beta_L5_3_L6_Hy3 = 1
  # Coefficients of L7
  beta_int_L7_Hy3 = -1 
  beta_L1_2_L7_Hy3 = -1
  beta_L2_2_L7_Hy3 = -1
  beta_L2_3_L7_Hy3 = -1
  beta_L2_4_L7_Hy3 = -1
  beta_L3_2_L7_Hy3 = -1
  beta_L3_3_L7_Hy3 = -1
  beta_L4_2_L7_Hy3 = -1
  beta_L5_2_L7_Hy3 = -1
  beta_L5_3_L7_Hy3 = -1
  
  Beta_Y0_Hy3 = c(beta_int_Y0_Hy3, beta_L1_2_Y0_Hy3, beta_L2_2_Y0_Hy3, beta_L2_3_Y0_Hy3, beta_L2_4_Y0_Hy3, beta_L3_2_Y0_Hy3, beta_L3_3_Y0_Hy3, beta_L4_2_Y0_Hy3, beta_L5_2_Y0_Hy3, beta_L5_3_Y0_Hy3)
  Beta_L6_Hy3 = c(beta_int_L6_Hy3, beta_L1_2_L6_Hy3, beta_L2_2_L6_Hy3, beta_L2_3_L6_Hy3, beta_L2_4_L6_Hy3, beta_L3_2_L6_Hy3, beta_L3_3_L6_Hy3, beta_L4_2_L6_Hy3, beta_L5_2_L6_Hy3, beta_L5_3_L6_Hy3)
  Beta_L7_Hy3 = c(beta_int_L7_Hy3, beta_L1_2_L7_Hy3, beta_L2_2_L7_Hy3, beta_L2_3_L7_Hy3, beta_L2_4_L7_Hy3, beta_L3_2_L7_Hy3, beta_L3_3_L7_Hy3, beta_L4_2_L7_Hy3, beta_L5_2_L7_Hy3, beta_L5_3_L7_Hy3)
  
  X_cat_dummies = as.matrix(cbind(1, L1_5_dummies))
  
  Y0_L6_7= array(data= NA, dim = c(n_sample, 3))
  for (i in 1: n_sample) {
    
    if (Hy[i] ==1){
      
      mu_Y0_Hy1 = X_cat_dummies %*% Beta_Y0_Hy1
      mu_L6_Hy1 = X_cat_dummies %*% Beta_L6_Hy1
      mu_L7_Hy1 = X_cat_dummies %*% Beta_L7_Hy1
      mu_Y0_L6_7_Hy1 = cbind(mu_Y0_Hy1, mu_L6_Hy1, mu_L6_Hy1)
      Sigma_Y0_L6_7_Hy1 = matrix(c(1, 0.3, 0.5, 
                                   0.3, 1, 0.5,
                                   0.5, 0.5, 1), 
                                 nrow=3, ncol = 3, byrow = T)
      Y0_L6_7[i,] = mvrnorm(n = 1, mu= mu_Y0_L6_7_Hy1[i,], Sigma = Sigma_Y0_L6_7_Hy1)
      
    } else if (Hy[i] == 2){
      
      mu_Y0_Hy2 = X_cat_dummies %*% Beta_Y0_Hy2
      mu_L6_Hy2 = X_cat_dummies %*% Beta_L6_Hy2
      mu_L7_Hy2 = X_cat_dummies %*% Beta_L7_Hy2
      mu_Y0_L6_7_Hy2 = cbind(mu_Y0_Hy2, mu_L6_Hy2, mu_L7_Hy2)
      Sigma_Y0_L6_7_Hy2 = matrix(c(1, 0.2, 0.5, 
                                   0.2, 1, 0.6,
                                   0.5, 0.6, 1), 
                                 nrow=3, ncol = 3, byrow = T)
      Y0_L6_7[i,] = mvrnorm(n = 1, mu= mu_Y0_L6_7_Hy2[i,], Sigma = Sigma_Y0_L6_7_Hy2)  
    } else if (Hy[i] == 3){
      
      mu_Y0_Hy3 = X_cat_dummies %*% Beta_Y0_Hy3
      mu_L6_Hy3 = X_cat_dummies %*% Beta_L6_Hy3
      mu_L7_Hy3 = X_cat_dummies %*% Beta_L7_Hy3
      mu_Y0_L6_7_Hy3 = cbind(mu_Y0_Hy3, mu_L6_Hy3, mu_L7_Hy3)
      Sigma_Y0_L6_7_Hy3 = matrix(c(1, 0.5, 0.1, 
                                   0.5, 1, 0,
                                   0.1, 0, 1), 
                                 nrow=3, ncol = 3, byrow = T)
      Y0_L6_7[i,] = mvrnorm(n = 1, mu= mu_Y0_L6_7_Hy3[i,], Sigma = Sigma_Y0_L6_7_Hy3)
    }
  }
  
  
  Y0 = Y0_L6_7[,1]
  L6 = Y0_L6_7[,2]
  L7 = Y0_L6_7[,3]
  L6_7 = cbind(L6, L7)
  L1_7 = cbind(L1, L2, L3, L4, L5, L6, L7)
  
  # Generate treatment indicator A
  expit_fn = function(x) { 1/(1+exp(-x)) } # inverse of the logit function
  prob_temp = expit_fn(-0.3 *(L1_7[,6]+L1_7[,7]) +3*ifelse(L5==1, 1, 0)+ 1.5*ifelse(L5==3, 1, 0)- 5*ifelse(L3==2, 1, 0) - ifelse(L3==3, 1, 0))
  TrueA = rbinom(n=n_sample, size=1, prob=prob_temp)
  
  # Generate observed potential outcome
  Y1 = Y0 + ATE_true
  Y_obs = TrueA*Y1+(1-TrueA)*Y0

	
	DATA[[i_rep]] = list(Y_obs=Y_obs, TrueA=TrueA, L1_7=L1_7)
  
}


# True mean of L6 and L7
mu_L6_rep  = mu_L7_rep = numeric(400)

for (i_rep in 1: 400) {
  
  print(paste0("i_rep=",i_rep))
  set.seed(100 + i_rep)
  
  mu_L6_rep[i_rep] = mean(DATA[[i_rep]]$L1_7[,6])
  mu_L7_rep[i_rep] = mean(DATA[[i_rep]]$L1_7[,7])

  
}

mu_L6_true = mean(mu_L6_rep)
mu_L7_true = mean(mu_L7_rep)

save(DATA, rep_seq, ATE_true, mu_L6_true, mu_L7_true, file = "1_DATA_BD.RData")


rm(list=ls()) 
load("1_DATA_BD.RData")

# Plots of continuous variables
i_rep = 1
Y = DATA[[i_rep]]$"Y_obs"
L1_7= DATA[[i_rep]]$L1_7
L6 = L1_7[, 6]
L7 = L1_7[, 7]

png(file=paste("pairwise_cont_sim2.png", sep = ""), width=1950,height=650,pointsize=20, res= 110)
par(mfrow=c(1,3), mai=c(1.5,1.5,0.4,0.2))
plot(L6, Y, xlab = "X6", ylab = "Y", pch = 16, col= "grey", xlim = c(-12, 12), ylim = c(-10, 12), cex= 1.5, cex.lab=1.5)
plot(L7, Y, xlab = "X7", ylab = "Y", pch = 16, col= "grey", xlim = c(-10, 12), ylim = c(-10, 12), cex= 1.5, cex.lab=1.5)
plot(L6, L7, pch = 16, col= "grey", xlab = "X6", ylab = "X7", xlim = c(-12, 12), ylim = c(-10, 12), cex.lab=1.5, cex= 1.5)
dev.off()

# Calculate true E(L6)
prob_z = c(0.05, 0.4, 0.3, 0.15, 0.1)

r1_z = c(0.50, 0.70, 0.40, 0.20, 0.30)
r2_z = c(0.20, 0.15, 0.20, 0.50, 0.20)
r3_z = c(0.30, 0.15, 0.40, 0.30, 0.50)

prob_r1 = r1_z %*% prob_z
prob_r2 = r2_z %*% prob_z
prob_r3 = r3_z %*% prob_z
prob_r = c(prob_r1, prob_r2, prob_r3)
sum(prob_r)

s1_z = c(0.15, 0.35, 0.25, 0.15, 0.45)
s2_z = c(0.10, 0.20, 0.15, 0.25, 0.15)
s3_z = c(0.35, 0.15, 0.15, 0.15, 0.15)
s4_z = c(0.15, 0.15, 0.35, 0.15, 0.10)
s5_z = c(0.25, 0.15, 0.10, 0.30, 0.15)

prob_s1 = s1_z %*% prob_z
prob_s2 = s2_z %*% prob_z
prob_s3 = s3_z %*% prob_z
prob_s4 = s4_z %*% prob_z
prob_s5 = s5_z %*% prob_z
prob_s = c(prob_s1, prob_s2, prob_s3, prob_s4, prob_s5)
sum(prob_s)

mean_L1_2 = c(0.40, 0.20, 0.45, 0.45, 0.75) %*% prob_s
mean_L2_2 = c(0.15, 0.35, 0.15, 0.35, 0.05) %*% prob_s
mean_L2_3 = c(0.20, 0.15, 0.05, 0.05, 0.15) %*% prob_s
mean_L2_4 = c(0.10, 0.05, 0.55, 0.35, 0.35) %*% prob_s
mean_L3_2 = c(0.10, 0.20, 0.30, 0.40, 0.20) %*% prob_s
mean_L3_3 = c(0.70, 0.60, 0.50, 0.30, 0.35) %*% prob_s
mean_L4_2 = c(0.50, 0.10, 0.50, 0.35, 0.25) %*% prob_s
mean_L5_2 = c(0.50, 0.10, 0.30, 0.20, 0.15) %*% prob_s
mean_L5_3 = c(0.20, 0.40, 0.20, 0.30, 0.10) %*% prob_s

mean_x_star = c(1, mean_L1_2, mean_L2_2, mean_L2_3, mean_L2_4, mean_L3_2, mean_L3_3, mean_L4_2, mean_L5_2, mean_L5_3)

beta_int_L6_Hy1 = 1.5 
beta_L1_2_L6_Hy1 = 0.5
beta_L2_2_L6_Hy1 = 1
beta_L2_3_L6_Hy1 = 2
beta_L2_4_L6_Hy1 = 1.5
beta_L3_2_L6_Hy1 = 2
beta_L3_3_L6_Hy1 = 0.5
beta_L4_2_L6_Hy1 = 0.5
beta_L5_2_L6_Hy1 = 0.3
beta_L5_3_L6_Hy1 = 2.5
Beta_L6_Hy1 = c(beta_int_L6_Hy1, beta_L1_2_L6_Hy1, beta_L2_2_L6_Hy1, beta_L2_3_L6_Hy1, beta_L2_4_L6_Hy1, 
                beta_L3_2_L6_Hy1, beta_L3_3_L6_Hy1, beta_L4_2_L6_Hy1, beta_L5_2_L6_Hy1, beta_L5_3_L6_Hy1)

beta_int_L6_Hy2 = -2.5 
beta_L1_2_L6_Hy2 = -0.1
beta_L2_2_L6_Hy2 = -3.1
beta_L2_3_L6_Hy2 = -2.2
beta_L2_4_L6_Hy2 = 1.5
beta_L3_2_L6_Hy2 = -2.0
beta_L3_3_L6_Hy2 = -0.5
beta_L4_2_L6_Hy2 = -2.5
beta_L5_2_L6_Hy2 = -1.3
beta_L5_3_L6_Hy2 = -0.5
Beta_L6_Hy2 = c(beta_int_L6_Hy2, beta_L1_2_L6_Hy2, beta_L2_2_L6_Hy2, beta_L2_3_L6_Hy2, beta_L2_4_L6_Hy2, 
                beta_L3_2_L6_Hy2, beta_L3_3_L6_Hy2, beta_L4_2_L6_Hy2, beta_L5_2_L6_Hy2, beta_L5_3_L6_Hy2)

beta_int_L6_Hy3 = 1 
beta_L1_2_L6_Hy3 = 1
beta_L2_2_L6_Hy3 = 1
beta_L2_3_L6_Hy3 = 1
beta_L2_4_L6_Hy3 = 1
beta_L3_2_L6_Hy3 = 1
beta_L3_3_L6_Hy3 = 1
beta_L4_2_L6_Hy3 = 1
beta_L5_2_L6_Hy3 = 1
beta_L5_3_L6_Hy3 = 1 
Beta_L6_Hy3 = c(beta_int_L6_Hy3, beta_L1_2_L6_Hy3, beta_L2_2_L6_Hy3, beta_L2_3_L6_Hy3, beta_L2_4_L6_Hy3, 
                beta_L3_2_L6_Hy3, beta_L3_3_L6_Hy3, beta_L4_2_L6_Hy3, beta_L5_2_L6_Hy3, beta_L5_3_L6_Hy3)

mu_L6_r1 = Beta_L6_Hy1 %*% mean_x_star
mu_L6_r2 = Beta_L6_Hy2 %*% mean_x_star
mu_L6_r3 = Beta_L6_Hy3 %*% mean_x_star

mean_L6 = c(mu_L6_r1, mu_L6_r2, mu_L6_r3) %*% prob_r
mean_L6 = as.numeric(mean_L6)

save(mean_L6, file = "mean_L6.Rdata")
