rm(list=ls())  

#####################################
# Dataset generated before deletion (i.e, introducing missingness)
#####################################

library(MASS) # for function mvrnorm (drawing from Multivariate Normal Distribution)

# True causal effect 
ATE_true = 1.5

rep_seq = c(1:1000) #1:40 # 1:40
n_sample = 500
DATA = NULL

expit_fn = function(x) { 1/(1+exp(-x)) } # inverse of the logit function

for (i_rep in rep_seq) {
  
  # i_rep = 1
  set.seed(100 + i_rep)
  
  ###############################################
  ########## Simulated data #####################
  ###############################################
  
  # Generating 4 continuous covariates from multivaraite normal with mean mu_L1_4 and covariance Sigma_L1_4
  mu_L1_4 = rep(0,4)
  Sigma_L1_4 = matrix(c(1, 0.3, 0.3, 0.3,                
                        0.3, 1, 0.3, 0.3, 
                        0.3, 0.3, 1, 0.3, 
                        0.3, 0.3, 0.3, 1), nrow=4, ncol = 4, byrow = T)
  L1_4 = mvrnorm(n = n_sample, mu= mu_L1_4, Sigma = Sigma_L1_4)
  dimnames(L1_4)[[2]] = c("L1","L2","L3","L4")
  
  # Genearting 2 cateogrical covariates 
  L5 = rep(0,n_sample)
  for (i in 1:n_sample){
    L5[i]= sample.int(3, size = 1, prob = c(0.5+L1_4[i,1]*0.05, 0.3-L1_4[i,1]*0.05, 0.2))
  }
  
  prob_L6 = expit_fn(0.5+L1_4[ ,3]*0.2)
  L6 = rbinom(n=n_sample, size=1, prob=prob_L6)
  
  L1_6 = cbind(L1_4, L5, L6)
  
  # Generating treatment indicator A given the covariates L1 to L6
  prob_temp = expit_fn(0.3 *(L1_4[,1]+L1_4[,2]+L1_4[,3]+L1_4[,4])+0.2*ifelse(L5==1, 1, 0)+0.1*ifelse(L5==3, 1, 0) - 0.3*ifelse(L6==1, 1, 0))
  mean(prob_temp)
  TrueA = rbinom(n=n_sample, size=1, prob=prob_temp)
  
  # Generating potential outcomes Y0 and Y1
  part1prob = exp(-2*(L1_4[,1]+1)^2)/(exp(-2*(L1_4[,1]+1)^2)+ exp(-2*(L1_4[,1]-2)^2))
  part1 = rbinom(n=n_sample,1,part1prob)
  mu1 = -4- 0.5*L1_4[,2]-L1_4[,3]+ 0.5*L1_4[,4]+0.3*ifelse(L5==1, 1, 0)-0.3*ifelse(L5==3, 1, 0)-ifelse(L6==1, 1, 0)
  mu2 = 4+ 0.5*(L1_4[,2])^2- 0.8*(L1_4[,3])*as.numeric(L1_4[,3]>0)+0.6*ifelse(L5==1, 1, 0) + 1.5*ifelse(L6==1, 1, 0)
  Y0_true = part1*rnorm(n=n_sample, mean=mu1, sd=1)+(1-part1)*rnorm(n=n_sample, mean=mu2, sd=4)
  ATE_true = 1.5
  Y1_true = Y0_true+ATE_true
  
  # Generating observed outcome
  Y_obs = TrueA*Y1_true+(1-TrueA)*Y0_true 
	
	DATA[[i_rep]] = list(Y_obs=Y_obs, TrueA=TrueA, L1_6=L1_6)
  
}

save(DATA, rep_seq, ATE_true, file = "1_DATA_BD.RData")