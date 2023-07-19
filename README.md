# Introduction
This GitHub repository is to house R codes for the paper "Bayesian causal inference for Observational Studies with Missingness in Covariates and Outcomes" by Huaiyu Zang, Hang Kim, Bin Huang, and Rhonda Szczesniak. Here, you can find all of our simulation codes. Our simulation codes are structured as follows: 

1. Generating simulated dataset using the data-generating process in Section 3 of the paper (see 1_DATA_BD.R). 
2. Running the proposed BNP model into the simulated data before introducing missing values (see 2_BNPc_BD.R).
3. Generating the missing data under the missing at random (MAR) assumption (see 3_DATA_MAR.R). 
4. Running the proposed BNP model into the simulated data with missing values (see 4_BNPc_MAR.R).  
5. Running existing causal models applied to the MICE imputed data (see 5_MICE_plus_existing_method_MAR.R). 
6. Running complete-case analyses with existing causal models (see 6_Existing_methods_MAR.R). 
7. Running the bartMachine model with missing covariates (see 7_bartMachine_MAR.R).

# Implementing the proposed model
To run our proposed BNP causal model, 
1. Install the "HCMMcausal" R package using the "HCMMcausal_1.5.2.tar.gz" source file. 
2. Make a data object using the "readData" function to load the data structure, in which one can input outcome variables in the "Response_var" argument, treatment indicator variable in the "trt_indicator" argument, continuous covariates in the "Cont_pred" argument and categorical covariates in the "Categ_pred" argument. 
3. Make a model object using the "createModel" function to load the data object. In the "createModel" function, one can specify the hyperparameters and upper bounds of mixture components. 
4. Run the "multipleImp" function for the proposed BNP causal model, where one can load the data and model objects and then save the results. In the "multipleImp" function, one can input the number of burn-in iterations in the "n_burnin" argument, the number of multiple imputations in the "m" argument, and the interval (number of iterations) between imputed data in the "interval_btw_Imp" argument.

# Examples

1. Running the proposed BNP model in case of missing covariates (Simulation 1)

```
rm(list=ls())  
library(MASS) # for function mvrnorm 
library(HCMMcausal)  # for running the proposed BNP causal model
expit_fn = function(x) { 1/(1+exp(-x)) }  # Inverse of the logit function 

###############################################
########## Data generation process ############
###############################################
ATE_true = 1.5 # define true causal effect
n_sample = 500 # define sample size
# Generate continuous covariates L1 to L4
mu_L1_4 = rep(0,4)
Sigma_L1_4 = matrix(c(1, 0.3, 0.3, 0.3,                
                      0.3, 1, 0.3, 0.3, 
                      0.3, 0.3, 1, 0.3, 
                      0.3, 0.3, 0.3, 1), nrow=4, ncol = 4, byrow = T)
L1_4 = mvrnorm(n = n_sample, mu= mu_L1_4, Sigma = Sigma_L1_4)
# Generate categorical covariates L5 and L6
L5 = rep(0,n_sample)
for (i in 1:n_sample){
  L5[i]= sample.int(3, size = 1, prob = c(0.5+L1_4[i,1]*0.05, 0.3-L1_4[i,1]*0.05, 0.2))
}
prob_L6 = expit_fn(0.5+L1_4[ ,3]*0.2)
L6 = rbinom(n=n_sample, size=1, prob=prob_L6)
L1_6 = cbind(L1_4, L5, L6)
# Generate treatment indicator A given the covariates L1 to L6
prob_temp = expit_fn(0.3 *(L1_4[,1]+L1_4[,2]+L1_4[,3]+L1_4[,4])+0.2*ifelse(L5==1, 1, 0) +0.1*ifelse(L5==3, 1, 0) - 0.3*ifelse(L6==1, 1, 0))
A = rbinom(n=n_sample, size=1, prob=prob_temp)
# Generate potential outcomes Y0 and Y1
part1prob = exp(-2*(L1_4[,1]+1)^2)/(exp(-2*(L1_4[,1]+1)^2)+ exp(-2*(L1_4[,1]-2)^2))
part1 = rbinom(n=n_sample,1,part1prob)
mu1 = (-4- 0.5*L1_4[,2]-L1_4[,3]+ 0.5*L1_4[,4]+0.3*ifelse(L5==1, 1, 0) -0.3*ifelse(L5==3, 1, 0)-ifelse(L6==1, 1, 0))
mu2 = (4+ 0.5*(L1_4[,2])^2- 0.8*(L1_4[,3])*as.numeric(L1_4[,3]>0) +0.6*ifelse(L5==1, 1, 0) + 1.5*ifelse(L6==1, 1, 0))
Y0 = (part1*rnorm(n=n_sample, mean=mu1, sd=1) +(1-part1)*rnorm(n=n_sample, mean=mu2, sd=4))
Y1 = Y0+ ATE_true
# Generate observed outcome
Y_obs = A*Y1+ (1-A)*Y0 

#####################################
## Introducing missingness ##########
#####################################
L1_4_missing_ind = array(0,c(n_sample,4))
L1_4_missing_ind[,1] = rbinom(n_sample, 1, expit_fn(-1 + 3*L1_4[,2]))
L1_4_missing_ind[,2] = rbinom(n_sample, 1, expit_fn(-1.5 + 3*L1_4[,3]))
L1_4_missing_ind[,4] = rbinom(n_sample, 1, expit_fn(-1.2 - 2*L1_4[,1] - L1_4[,3]))
L5_missing_ind = rbinom(n_sample, 1, expit_fn(-1.5 + 3*L1_4[,1] - 1.5*L1_4[,3]))
L6_missing_ind = rbinom(n_sample, 1, expit_fn(-1.2 + 2.2*L1_4[,1] - 0.5*L1_4[,3]))
L1_6_missing_ind = cbind(L1_4_missing_ind, L5_missing_ind, L6_missing_ind)
L1_6[L1_6_missing_ind==1] = NA

###############################################
### Running proposed BNP causal model #########
###############################################
obs_response = Y_obs 
Incomplete_cont_L = L1_6_missing[,1:4]
Incomplete_cat_L = L1_6_missing[,5:6]
n_burnin = 2000 ; m_Imp = 10 ; interval_btw_Imp = 1000 #  The ATE estimate is obtained by averaging over 10,000 MCMC iterations after 2,000 burn-in iterations
data_obj = readData(Response_var=obs_response, trt_indicator= A, Cont_pred=Incomplete_cont_L, Categ_pred=Incomplete_cat_L, RandomSeed=100)
model_obj = createModel(data_obj)
result_obj = multipleImp(data_obj= data_obj, model_obj= model_obj, n_burnin = n_burnin, m= m_Imp, interval_btw_Imp= interval_btw_Imp, show_iter=TRUE)

# Calculate posterior mean and standard deviation for ATE estimate using the proposed model 
ATE_BNP_causal = mean(result_obj$est_delta)
SD_ATE_BNP_causal = sd(result_obj$est_delta)                    
```

2. Running the proposed BNP model in case of missing both covariates and outcome (Simulation 3). 

```
library( HCMMcausal )  # for running BNP causal 
load("Simulation3_dataset.Rdata") # Download the data from the folder of Simulation 3
# In this dataset, Y1 and Y2 are outcomes; TrueA denote treatment indicator, L1 to L5 are categorical covariates, L6 and L7 are continuous covariates
# > head(Simulation3_dataset)
#            Y1       Y2 TrueA        L6          L7 L1 L2 L3 L4 L5
# [1,] 2.794514 1.941280     0  4.499388  4.14310808  1  2  3  1  3
# [2,]       NA 4.034288     1        NA  0.08950193  1  4  1  2  1
# [3,]       NA 7.513693     0  6.707670  6.27963788  1  4  2  1  2
# [4,]       NA       NA     0  2.288125  2.95212193  2  2  1  1  1
# [5,]       NA       NA     0        NA -3.30222185  2  4  2  1  1
# [6,]       NA       NA     1 -5.409038  1.26828624  1  1  3  1  2

# BNP model specification 
obs_response = Simulation3_dataset[, c("Y1", "Y2")]
TrueA = Simulation3_dataset[, c("TrueA")]
Cont_covariate = Simulation3_dataset[, c("L6", "L7")]
Cat_covariate = Simulation3_dataset[, c("L1", "L2", "L3", "L4", "L5")]
n_burnin = 2000 ; m_Imp = 10 ; interval_btw_Imp = 1000
data_obj = readData(Response_var=obs_response, trt_indicator=TrueA, Cont_pred=Cont_covariate, Categ_pred=Cat_covariate, RandomSeed=100)
model_obj = createModel(data_obj)	
result_obj = multipleImp(data_obj, model_obj, n_burnin, m_Imp, interval_btw_Imp, show_iter=T)

# Saving posterior mean and standard deviation for ATE estimates
ATE_BNPc_vec = apply(result_obj$est_delta,2,mean) # ATE_BNPc_vec[1] is ATE estimate for Y1 and ATE_BNPc_vec[2] is ATE estimate for Y2
SD_ATE_BNPc_vec = apply(result_obj$est_delta,2,sd)
```


# Contact
Don't hesitate to contact Huaiyu Zang with any questions, complaints, requests, etc., via email: huaiyuzang [at] gmail [dot] com.
