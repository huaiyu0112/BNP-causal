# BNP-causal
R Codes for simulations and data analysis in Zang et.al. "A Bayesian Causal Inference Approach in Observational Studies with Missingness in Covariates and Outcomes".


To run our proposed BNP causal model, (1) First, we install "HCMMcausal" package using "HCMMcausal_1.5.2.tar.gz" file, (2) Second, we make the data object using readData function to load the data structure, in which we input response variable in "Response_var" argument, treatment indicator variable in "trt_indicator" argument, continuous covariates in "Cont_pred" argument and categorical covariates in "Categ_pred" argument. (3) Third, we make a model object using createModel function to load the data object. In createModel function, we could also specify the hyperparameters and upper bounds of mixture components. (4) Finally, we run multipleImp function for proposed BNP causal model where we load the data and model object and then save the results. In multipleImp function, we could input the number of burn-in iterations in "n_burnin" argument, number of multiple imputations in "m" argument, and interval (number of iterations) between imputed data in "interval_btw_Imp" argument.

See a sample R code of Simulation 1 for the illustration of implementing the HCMMcausal package below. 

```
#################################################################################
rm(list=ls())  
library(MASS) # for function mvrnorm 
library(HCMMcausal)  # for running proposed BNP causal model
sessionInfo()
set.seed(100)
# Inverse of the logit function
expit_fn = function(x) { 1/(1+exp(-x)) } 
###############################################
########## Simulated data #####################
###############################################
# True causal effect 
ATE_true = 1.5
# Define number of observations for each dataset 
n_sample = 500
# Four continuous confounders 
mu_L1_4 = rep(0,4)
Sigma_L1_4 = matrix(c(1, 0.3, 0.3, 0.3,                
                      0.3, 1, 0.3, 0.3, 
                      0.3, 0.3, 1, 0.3, 
                      0.3, 0.3, 0.3, 1), nrow=4, ncol = 4, byrow = T)
L1_4 = mvrnorm(n = n_sample, mu= mu_L1_4, Sigma = Sigma_L1_4)
dimnames(L1_4)[[2]] = c("L1","L2","L3","L4")
# Two cateogrical confounders
L5 = rep(0,n_sample)
for (i in 1:n_sample){
  L5[i]= sample.int(3, size = 1, prob = c(0.5+L1_4[i,1]*0.05,
                                          0.3-L1_4[i,1]*0.05, 0.2))
}
prob_L6 = expit_fn(0.5+L1_4[ ,3]*0.2)
L6 = rbinom(n=n_sample, size=1, prob=prob_L6)
L1_6 = cbind(L1_4, L5, L6)
# Treatment indicator 
prob_temp = expit_fn(0.3 *(L1_4[,1]+L1_4[,2]+L1_4[,3]+L1_4[,4])+0.2*ifelse(L5==1, 1, 0)
                     +0.1*ifelse(L5==3, 1, 0) - 0.3*ifelse(L6==1, 1, 0))
A = rbinom(n=n_sample, size=1, prob=prob_temp)
# Potential outcomes 
part1prob = exp(-2*(L1_4[,1]+1)^2)/(exp(-2*(L1_4[,1]+1)^2)+ exp(-2*(L1_4[,1]-2)^2))
part1 = rbinom(n=n_sample,1,part1prob)
mu1 = (-4- 0.5*L1_4[,2]-L1_4[,3]+ 0.5*L1_4[,4]+0.3*ifelse(L5==1, 1, 0)
       -0.3*ifelse(L5==3, 1, 0)-ifelse(L6==1, 1, 0))
mu2 = (4+ 0.5*(L1_4[,2])^2- 0.8*(L1_4[,3])*as.numeric(L1_4[,3]>0)
       +0.6*ifelse(L5==1, 1, 0) + 1.5*ifelse(L6==1, 1, 0))
Y0 = (part1*rnorm(n=n_sample, mean=mu1, sd=1)
           +(1-part1)*rnorm(n=n_sample, mean=mu2, sd=4))
Y1 = Y0+ ATE_true
# Observed outcome
Y_obs = A*Y1+ (1-A)*Y0 
#####################################
########## Missingness ##############
#####################################
L1_4_missing_ind = array(0,c(n_sample,4))
L1_4_missing_ind[,2] = rbinom(n_sample, 1, expit_fn(-1.5 + 3*L1_4[,3]))
L1_4_missing_ind[,4] = rbinom(n_sample, 1, expit_fn(-1.2 - 2*L1_4[,1] - L1_4[,3]))
L5_missing_ind = rbinom(n_sample, 1, expit_fn(-1.5 + 3*L1_4[,1] - 1.5*L1_4[,3]))
L6_missing_ind = rbinom(n_sample, 1, expit_fn(-1.2 + 2.2*L1_4[,1] - 0.5*L1_4[,3]))
L1_6_missing_ind = cbind(L1_4_missing_ind, L5_missing_ind, L6_missing_ind)
L1_6_missing = L1_6 
for (i_p in (1:dim( L1_6)[2])){
  L1_6_missing[,i_p] = ifelse( L1_6_missing_ind[,i_p]==1, NA, L1_6[,i_p])
}
###############################################
###Running proposed BNP causal model ##########
###############################################
obs_response = Y_obs 
Incomplete_cont_L = L1_6_missing[,1:4]
Incomplete_cat_L = L1_6_missing[,5:6]
n_burnin=2000; m_Imp=10; interval_btw_Imp=200
data_obj = readData(Response_var=obs_response, trt_indicator= A, 
                    Cont_pred=Incomplete_cont_L, Categ_pred=Incomplete_cat_L,
                    RandomSeed=100)
# Default mixture components: max_R_S_K = c(30, 50, 20)
# Default hyperprior values: a_R = 0.5, b_R = 0.5, 
# a_S = 0.5, b_S = 0.5, a_K = 0.5, b_K = 0.5,  
# psi_0 = 1, b_theta = 100, a_tau = 0.5, b_tau = 0.5, b_delta = 100
model_obj = createModel(data_obj)
result_obj = multipleImp(data_obj= data_obj, model_obj= model_obj,
                         n_burnin = n_burnin, m= m_Imp,
                         interval_btw_Imp= interval_btw_Imp, show_iter=TRUE)
# Posterior mean and standard deviation of est. ATE from the proposed BNP causal model 
ATE_BNP_causal = mean(result_obj$est_delta)
SD_ATE_BNP_causal = sd(result_obj$est_delta)
#################################################################################                     
```


