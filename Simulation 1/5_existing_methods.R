rm(list=ls()) 

library(mice) # for mice function
library(MASS) # for "mvrnorm"
library(mvtnorm) # for dmvnorm
library(bartCause) # BART causal 
library(ipw) # for ipwpoint function
library(WeightIt) # for weightit function
library(cobalt) # for get.w function 
library(survey) # for svyglm function
library("SuperLearner") # used cross validation to weigh different prediction algorithms for TMLE
library(tmle) # for tmle function
library(norm) # for mi.inference function

load("1_DATA_BD.RData")
load("3_DATA_MAR.RData")

# True causal effect 
ATE_true = 1.5

expit_fn = function(x) { 1/(1+exp(-x)) } # inverse of the logit function

# How many times the simulation are repeated 
start_seq =0
rep_seq = c(1:400) #1:40 # 1:40
n_seq = length(rep_seq)

# Define number of observations for each dataset 
n_sample = 500

dir.create("6_Off_Shelf_MAR")

for (i_rep in rep_seq) {
  
  # i_rep = 1
  set.seed(100 + i_rep + start_seq)
  seed = 100 + i_rep + start_seq
  print(paste0("i_rep=",i_rep))
  print(paste0("seed=",seed))
  
  Y_obs = incompleteDATA[[i_rep]]$Y_obs 
  A = incompleteDATA[[i_rep]]$TrueA
  Incomplete_Y = incompleteDATA[[i_rep]]$L1_6[,1:4]
  Incomplete_X = incompleteDATA[[i_rep]]$L1_6[,5:6]
  L1_6_missing = cbind(Incomplete_Y, Incomplete_X)
  L1_6 = DATA[[i_rep]]$L1_6[,1:6]
  
  ###############################################
  ## Linear regression  ##
  ###############################################
  lm_missing = lm(Y_obs ~ as.factor(A) + L1_6_missing[,1] + L1_6_missing[,2] + L1_6_missing[,3] + L1_6_missing[,4]+
                    as.factor(L1_6_missing[,5]) + as.factor(L1_6_missing[,6]))
  ATE_lm = summary(lm_missing)$coefficients[2,1]
  sd_ATE_lm = summary(lm_missing)$coefficients[2,2]
  Check_ATE_lm = (confint(lm_missing)[2,][1] < ATE_true & ATE_true < confint(lm_missing)[2,][2])
  
  ###############################################
  ## BART  ##
  ###############################################
  
  # ATE 
  bartc_ATE_missing = bartc(response= Y_obs, treatment= A, confounders= L1_6_missing, estimand   = c("ate"))
  ATE_BART = fitted(bartc_ATE_missing)
  sd_ATE_BART = sd(extract(bartc_ATE_missing))
  Check_ATE_BART = (summary(bartc_ATE_missing)$"estimates"[3] < ATE_true & ATE_true < summary(bartc_ATE_missing)$"estimates"[4])
  
  ###############################################
  ## IPTW ##
  ###############################################
  mydata = data.frame(Y_obs, A, L1_6_missing)
  mydata$L5 = factor(mydata$L5)
  mydata$L6 = factor(mydata$L6)
  L1_6_complete.cases_ind = complete.cases(L1_6_missing)
  mydata_complete.cases = mydata[L1_6_complete.cases_ind,] 
  
  # ATE 
  weightmodel_ATE= weightit(A ~ L1 + L2+ L3 + L4+ factor(L5) + factor(L6), data = mydata_complete.cases, method = "ps", estimand = "ATE")
  mydata_complete.cases$wt_ATE = get.w(weightmodel_ATE)
  msm_ATE <- (svyglm(Y_obs ~ A, design = svydesign(~ 1, weights = ~ wt_ATE, data =mydata_complete.cases)))
  ATE_IPTW = coef(msm_ATE)[2]
  sd_ATE_IPTW = summary(msm_ATE)$"coefficients"[2,2]
  Check_ATE_IPTW = (confint(msm_ATE)[2,][1] < ATE_true & ATE_true < confint(msm_ATE)[2,][2])
  
  ###############################################
  ## TMLE  ##
  ###############################################
  Y_obs_complete.cases = mydata_complete.cases[,"Y_obs"]
  A_complete.cases = mydata_complete.cases[,"A"]
  L1_6_complete.cases = mydata_complete.cases[,c(3:8)]
  
  res.tmle = tmle(Y = Y_obs_complete.cases, A = A_complete.cases , W= L1_6_complete.cases,
                  Qform= "Y ~ A +  L1 + L2 + L3 + L4 + factor(L5) + factor(L6)",
                  gform = "A ~ L1 + L2 + L3 + L4 + factor(L5) + factor(L6)",
                  fluctuation = "logistic")
  # ATE 
  ATE_TMLE = res.tmle$estimates$ATE$psi
  sd_ATE_TMLE = sqrt(res.tmle$estimates$ATE$var.psi)
  Check_ATE_TMLE = (res.tmle$estimates$ATE$CI[1] < ATE_true & ATE_true < res.tmle$estimates$ATE$CI[2])
  
  save(ATE_BART, ATE_lm, ATE_IPTW, ATE_TMLE, 
       sd_ATE_BART, sd_ATE_lm, sd_ATE_IPTW, sd_ATE_TMLE,
       Check_ATE_BART, Check_ATE_lm, Check_ATE_IPTW, Check_ATE_TMLE, 
       file = paste("6_Off_Shelf_MAR/res_", i_rep , "_off_shelf_cc",".RData", sep = ""))
  
}


rm(list=ls())
load("1_DATA_BD.RData")

rep_seq= 1:400
ATE_BART_rep = ATE_lm_rep = ATE_IPTW_rep = ATE_TMLE_rep = numeric(length(rep_seq))
sd_ATE_BART_rep = sd_ATE_lm_rep = sd_ATE_IPTW_rep = sd_ATE_TMLE_rep = numeric(length(rep_seq))
IsCovered_ATE_BART_rep = IsCovered_ATE_lm_rep = IsCovered_ATE_IPTW_rep = IsCovered_ATE_TMLE_rep = numeric(length(rep_seq))

for (i_rep in rep_seq) {
  
  # i_rep = 1
  load( paste("6_Off_Shelf_MAR/res_", i_rep , "_off_shelf_cc",".RData", sep = ""))
  
  #################################################
  # Repeated sampling results
  #################################################  
  # Est. ATE
  ATE_BART_rep[i_rep] = ATE_BART
  ATE_lm_rep[i_rep] = ATE_lm
  ATE_IPTW_rep[i_rep] = ATE_IPTW
  ATE_TMLE_rep[i_rep] = ATE_TMLE
  
  # S.E of Est. ATE
  sd_ATE_BART_rep[i_rep] = sd_ATE_BART
  sd_ATE_lm_rep[i_rep] = sd_ATE_lm
  sd_ATE_IPTW_rep[i_rep] = sd_ATE_IPTW
  sd_ATE_TMLE_rep[i_rep] = sd_ATE_TMLE
  
  # 95% C.I of ATE
  IsCovered_ATE_BART_rep[i_rep] = Check_ATE_BART
  IsCovered_ATE_lm_rep[i_rep] = Check_ATE_lm
  IsCovered_ATE_IPTW_rep[i_rep] = Check_ATE_IPTW 
  IsCovered_ATE_TMLE_rep[i_rep] = Check_ATE_TMLE
  
}

# ###################################
# # Monte Carlo simulation results 
# ###################################
# BART
avg_ATE_BART = mean(ATE_BART_rep)
bias_ATE_BART = mean(ATE_BART_rep -ATE_true)
emp_sd_ATE_BART = sd(ATE_BART_rep)
avg_sd_ATE_BART = mean(sd_ATE_BART_rep)
coverage_ATE_BART= mean(IsCovered_ATE_BART_rep)
RMSE_ATE_BART = sqrt(mean((ATE_BART_rep - ATE_true)^2))
MAE_ATE_BART = median(abs(ATE_BART_rep - ATE_true))

# linear regression
avg_ATE_lm = mean(ATE_lm_rep)
bias_ATE_lm = mean(ATE_lm_rep -ATE_true)
emp_sd_ATE_lm = sd(ATE_lm_rep)
avg_sd_ATE_lm = mean(sd_ATE_lm_rep)
coverage_ATE_lm= mean(IsCovered_ATE_lm_rep)
RMSE_ATE_lm = sqrt(mean((ATE_lm_rep - ATE_true)^2))
MAE_ATE_lm = median(abs(ATE_lm_rep - ATE_true))

# IPTW 
avg_ATE_IPTW = mean(ATE_IPTW_rep)
bias_ATE_IPTW = mean(ATE_IPTW_rep -ATE_true)
emp_sd_ATE_IPTW = sd(ATE_IPTW_rep)
avg_sd_ATE_IPTW = mean(sd_ATE_IPTW_rep)
coverage_ATE_IPTW= mean(IsCovered_ATE_IPTW_rep)
RMSE_ATE_IPTW = sqrt(mean((ATE_IPTW_rep - ATE_true)^2))
MAE_ATE_IPTW = median(abs(ATE_IPTW_rep - ATE_true))

# TMLE
avg_ATE_TMLE = mean(ATE_TMLE_rep)
bias_ATE_TMLE = mean(ATE_TMLE_rep -ATE_true)
emp_sd_ATE_TMLE = sd(ATE_TMLE_rep)
avg_sd_ATE_TMLE = mean(sd_ATE_TMLE_rep)
coverage_ATE_TMLE= mean(IsCovered_ATE_TMLE_rep)
RMSE_ATE_TMLE = sqrt(mean((ATE_TMLE_rep - ATE_true)^2))
MAE_ATE_TMLE = median(abs(ATE_TMLE_rep - ATE_true))

method_causal = c("BART cc", "LM cc", "IPTW cc", "TMLE cc")
effect_causal = rep("ATE", times= 4)
point_est_causal = c(avg_ATE_BART, avg_ATE_lm, avg_ATE_IPTW, avg_ATE_TMLE)
bias_causal = c(bias_ATE_BART, bias_ATE_lm, bias_ATE_IPTW, bias_ATE_TMLE)
emp_sd_causal = c(emp_sd_ATE_BART, emp_sd_ATE_lm, emp_sd_ATE_IPTW, emp_sd_ATE_TMLE)
avg_sd_causal = c(avg_sd_ATE_BART, avg_sd_ATE_lm, avg_sd_ATE_IPTW, avg_sd_ATE_TMLE)
coverage_causal = c(coverage_ATE_BART, coverage_ATE_lm, coverage_ATE_IPTW, coverage_ATE_TMLE)
RMSE_causal = c(RMSE_ATE_BART, RMSE_ATE_lm, RMSE_ATE_IPTW, RMSE_ATE_TMLE)
MAE_causal = c(MAE_ATE_BART, MAE_ATE_lm, MAE_ATE_IPTW, MAE_ATE_TMLE)

off_shelf_causal_res_cc_MS = data.frame(method = method_causal, effect= effect_causal, 
                                        point_est= point_est_causal, emp_sd= emp_sd_causal, RMSE = RMSE_causal, 
                                        MAE = MAE_causal, coverage = coverage_causal, avg_sd = avg_sd_causal)

write.csv(off_shelf_causal_res_cc_MS, file = "6_Off_Shelf_MAR_causal_MS.csv")

