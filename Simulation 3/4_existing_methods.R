rm(list=ls())  

library(bartCause) # BART causal 
library(ipw) # for ipwpoint function
library(WeightIt)
library(cobalt)
library(survey) # for svyglm function
library("SuperLearner") # used cross validation to weigh different prediction algorithms for TMLE
library(tmle) # for tmle function
library(gam) # Generalized Additive Models (GAM) prediction algorithm 
library(glmnet) # Glmnet prediction algorithm
library(randomForest) # Random Forest prediction algorithm
library(norm)

load("1_DATA.RData")

n_burnin = 2000 ; m = 10 ; interval_btw = 200
rep_seq = 1: 400
n_seq = length(rep_seq)
n_sample = dim(IncompleteDATA[[1]]$L1_7)[[1]]

ATE_Y_true = 1.5
ATE_X_true = 0.5

ATE_Y_BART_rep = ATE_Y_lm_rep = ATE_Y_IPTW_rep = ATE_Y_TMLE_rep = numeric(length(rep_seq))
sd_ATE_Y_BART_rep = sd_ATE_Y_lm_rep = sd_ATE_Y_IPTW_rep = sd_ATE_Y_TMLE_rep = numeric(length(rep_seq))
IsCovered_ATE_Y_BART_rep = IsCovered_ATE_Y_lm_rep = IsCovered_ATE_Y_IPTW_rep = IsCovered_ATE_Y_TMLE_rep = numeric(length(rep_seq))

ATE_X_BART_rep = ATE_X_lm_rep = ATE_X_IPTW_rep = ATE_X_TMLE_rep = numeric(length(rep_seq))
sd_ATE_X_BART_rep = sd_ATE_X_lm_rep = sd_ATE_X_IPTW_rep = sd_ATE_X_TMLE_rep = numeric(length(rep_seq))
IsCovered_ATE_X_BART_rep = IsCovered_ATE_X_lm_rep = IsCovered_ATE_X_IPTW_rep = IsCovered_ATE_X_TMLE_rep = numeric(length(rep_seq))

dir.create("6_Off_Shelf_Est_CC")

for (i_rep in rep_seq){
  
  # i_rep = 1
  print(paste0("i_rep=",i_rep))
  set.seed(100 + i_rep)
  
  Y_missing = IncompleteDATA[[i_rep]]$Outcome[, 1] 
  X_missing = IncompleteDATA[[i_rep]]$Outcome[, 2] 
  A = IncompleteDATA[[i_rep]]$TrueA
  Incomplete_Y = IncompleteDATA[[i_rep]]$L1_7[,6:7]
  Incomplete_X = IncompleteDATA[[i_rep]]$L1_7[,1:5]
  L1_7_missing = cbind(Incomplete_X, Incomplete_Y)
  dat_sim2 = as.data.frame(cbind(Y_missing, X_missing, A, L1_7_missing))
  
  dat_sim2_cc = dat_sim2[complete.cases(dat_sim2),]
  
  ########################
  ## Linear regression  ##
  ########################
  lm_Y_missing = lm(Y_missing ~ as.factor(A) + factor(L1) + factor(L2) + factor(L3) + factor(L4) + factor(L5) + L6 + L7, data = dat_sim2 )
  ATE_Y_lm = summary(lm_Y_missing)$coefficients[2,1]
  sd_ATE_Y_lm  = summary(lm_Y_missing)$coefficients[2,2]
  Check_ATE_Y_lm = (confint(lm_Y_missing)[2,][1] < ATE_Y_true & ATE_Y_true < confint(lm_Y_missing)[2,][2])
  
  lm_X_missing = lm(X_missing ~ as.factor(A) + factor(L1) + factor(L2) + factor(L3) + factor(L4) + factor(L5) + L6 + L7, data = dat_sim2 )
  ATE_X_lm = summary(lm_X_missing)$coefficients[2,1]
  sd_ATE_X_lm  = summary(lm_X_missing)$coefficients[2,2]
  Check_ATE_X_lm = (confint(lm_X_missing)[2,][1] < ATE_X_true & ATE_X_true < confint(lm_X_missing)[2,][2])
  
  ###########
  ## BART  ##
  ###########
  bartc_ATE_Y_res = bartc(response= Y_missing, treatment= A, confounders= factor(L1) + factor(L2) + factor(L3) + factor(L4) + factor(L5) + L6 + L7, 
                          estimand  = c("ate"), data = dat_sim2_cc)
  ATE_Y_BART = fitted(bartc_ATE_Y_res)
  sd_ATE_Y_BART = sd(extract(bartc_ATE_Y_res ))
  Check_ATE_Y_BART = (summary(bartc_ATE_Y_res)$"estimates"[3] < ATE_Y_true & ATE_Y_true < summary(bartc_ATE_Y_res)$"estimates"[4])
  
  bartc_ATE_X_res = bartc(response= X_missing, treatment= A, confounders= factor(L1) + factor(L2) + factor(L3) + factor(L4) + factor(L5) + L6 + L7, 
                          estimand   = c("ate"), data = dat_sim2_cc)
  ATE_X_BART = fitted(bartc_ATE_X_res)
  sd_ATE_X_BART = sd(extract(bartc_ATE_X_res ))
  Check_ATE_X_BART = (summary(bartc_ATE_X_res)$"estimates"[3] < ATE_X_true & ATE_X_true < summary(bartc_ATE_X_res)$"estimates"[4])
  
  
  ###########
  ## IPTW ##
  ##########
  
  weightmodel_ATE= weightit(A ~ factor(L1) + factor(L2) + factor(L3) + factor(L4) + factor(L5) + L6 + L7, data = dat_sim2_cc, method = "ps", estimand = "ATE")
  dat_sim2_cc$wt_ATE = get.w(weightmodel_ATE)
  msm_ATE_Y = (svyglm(Y_missing ~ A, design = svydesign(~ 1, weights = ~ wt_ATE, data =dat_sim2_cc)))
  ATE_Y_IPTW = coef(msm_ATE_Y)[2]
  sd_ATE_Y_IPTW = summary(msm_ATE_Y)$"coefficients"[2,2]
  Check_ATE_Y_IPTW = (confint(msm_ATE_Y)[2,][1] < ATE_Y_true & ATE_Y_true < confint(msm_ATE_Y)[2,][2])
  
  msm_ATE_X = (svyglm(X_missing ~ A, design = svydesign(~ 1, weights = ~ wt_ATE, data =dat_sim2_cc)))
  ATE_X_IPTW = coef(msm_ATE_X)[2]
  sd_ATE_X_IPTW = summary(msm_ATE_X)$"coefficients"[2,2]
  Check_ATE_X_IPTW = (confint(msm_ATE_X)[2,][1] < ATE_X_true & ATE_X_true < confint(msm_ATE_X)[2,][2])
  
  ###########
  ## TMLE  ##
  ###########
  # R-package tmle (base implementation includes SL.step, SL.glm and SL.glm.interaction)
  
  res.tmle_Y = tmle(Y = dat_sim2_cc$Y_missing, A = dat_sim2_cc$A, W= dat_sim2_cc[, c("L1", "L2", "L3", "L4", "L5", "L6", "L7")],
                    Qform= "Y ~ A +  factor(L1) + factor(L2) + factor(L3) + factor(L4) + factor(L5) + L6 + L7",
                    gform = "A ~ factor(L1) + factor(L2) + factor(L3) + factor(L4) + factor(L5) + L6 + L7",
                    # g.SL.library = c("SL.mean", "SL.glm", "SL.step"),
                    # Q.SL.library = c("SL.mean", "SL.glm", "SL.step"),
                    fluctuation = "logistic")
  
  ATE_Y_TMLE = res.tmle_Y$estimates$ATE$psi
  sd_ATE_Y_TMLE = sqrt(res.tmle_Y$estimates$ATE$var.psi)
  Check_ATE_Y_TMLE = (res.tmle_Y$estimates$ATE$CI[1] < ATE_Y_true & ATE_Y_true < res.tmle_Y$estimates$ATE$CI[2])
  
  res.tmle_X = tmle(Y = dat_sim2_cc$X_missing, A = dat_sim2_cc$A, W= dat_sim2_cc[, c("L1", "L2", "L3", "L4", "L5", "L6", "L7")],
                    Qform= "Y ~ A +  factor(L1) + factor(L2) + factor(L3) + factor(L4) + factor(L5) + L6 + L7",
                    gform = "A ~ factor(L1) + factor(L2) + factor(L3) + factor(L4) + factor(L5) + L6 + L7",
                    # g.SL.library = c("SL.mean", "SL.glm", "SL.step"),
                    # Q.SL.library = c("SL.mean", "SL.glm", "SL.step"),
                    fluctuation = "logistic")
  
  ATE_X_TMLE = res.tmle_X$estimates$ATE$psi
  sd_ATE_X_TMLE = sqrt(res.tmle_X$estimates$ATE$var.psi)
  Check_ATE_X_TMLE = (res.tmle_X$estimates$ATE$CI[1] < ATE_X_true & ATE_X_true < res.tmle_X$estimates$ATE$CI[2])
  
  save(ATE_Y_BART, ATE_Y_lm, ATE_Y_IPTW, ATE_Y_TMLE, 
       sd_ATE_Y_BART, sd_ATE_Y_lm, sd_ATE_Y_IPTW, sd_ATE_Y_TMLE,
       Check_ATE_Y_BART, Check_ATE_Y_lm, Check_ATE_Y_IPTW, Check_ATE_Y_TMLE, 
       ATE_X_BART, ATE_X_lm, ATE_X_IPTW, ATE_X_TMLE, 
       sd_ATE_X_BART, sd_ATE_X_lm, sd_ATE_X_IPTW, sd_ATE_X_TMLE,
       Check_ATE_X_BART, Check_ATE_X_lm, Check_ATE_X_IPTW, Check_ATE_X_TMLE,
       file = paste("6_Off_Shelf_Est_CC/res_", i_rep , "_off_shelf_cc",".RData", sep = ""))
  
} # for (i_rep)


rm(list=ls())
load("1_DATA.RData")

rep_seq= 1:400
ATE_Y_true= 1.5 
ATE_X_true = 0.5
ATE_Y_BART_rep = ATE_Y_lm_rep = ATE_Y_IPTW_rep = ATE_Y_TMLE_rep = numeric(length(rep_seq))
sd_ATE_Y_BART_rep = sd_ATE_Y_lm_rep = sd_ATE_Y_IPTW_rep = sd_ATE_Y_TMLE_rep = numeric(length(rep_seq))
IsCovered_ATE_Y_BART_rep = IsCovered_ATE_Y_lm_rep = IsCovered_ATE_Y_IPTW_rep = IsCovered_ATE_Y_TMLE_rep = numeric(length(rep_seq))

ATE_X_BART_rep = ATE_X_lm_rep = ATE_X_IPTW_rep = ATE_X_TMLE_rep = numeric(length(rep_seq))
sd_ATE_X_BART_rep = sd_ATE_X_lm_rep = sd_ATE_X_IPTW_rep = sd_ATE_X_TMLE_rep = numeric(length(rep_seq))
IsCovered_ATE_X_BART_rep = IsCovered_ATE_X_lm_rep = IsCovered_ATE_X_IPTW_rep = IsCovered_ATE_X_TMLE_rep = numeric(length(rep_seq))

for (i_rep in rep_seq) {
  
  # i_rep = 1
  load(paste("6_Off_Shelf_Est_CC/res_", i_rep , "_off_shelf_cc",".RData", sep = ""))
  
  #################################################
  # Repeated sampling results
  #################################################  
  # Est. ATE_Y
  ATE_Y_BART_rep[i_rep] = ATE_Y_BART
  ATE_Y_lm_rep[i_rep] = ATE_Y_lm
  ATE_Y_IPTW_rep[i_rep] = ATE_Y_IPTW
  ATE_Y_TMLE_rep[i_rep] = ATE_Y_TMLE
  
  # S.E of Est. ATE_Y
  sd_ATE_Y_BART_rep[i_rep] = sd_ATE_Y_BART
  sd_ATE_Y_lm_rep[i_rep] = sd_ATE_Y_lm
  sd_ATE_Y_IPTW_rep[i_rep] = sd_ATE_Y_IPTW
  sd_ATE_Y_TMLE_rep[i_rep] = sd_ATE_Y_TMLE
  
  # 95% C.I of ATE_Y
  IsCovered_ATE_Y_BART_rep[i_rep] = Check_ATE_Y_BART
  IsCovered_ATE_Y_lm_rep[i_rep] = Check_ATE_Y_lm
  IsCovered_ATE_Y_IPTW_rep[i_rep] = Check_ATE_Y_IPTW 
  IsCovered_ATE_Y_TMLE_rep[i_rep] = Check_ATE_Y_TMLE
  
  ATE_X_BART_rep[i_rep] = ATE_X_BART
  ATE_X_lm_rep[i_rep] = ATE_X_lm
  ATE_X_IPTW_rep[i_rep] = ATE_X_IPTW
  ATE_X_TMLE_rep[i_rep] = ATE_X_TMLE
  
  # S.E of Est. ATE_X
  sd_ATE_X_BART_rep[i_rep] = sd_ATE_X_BART
  sd_ATE_X_lm_rep[i_rep] = sd_ATE_X_lm
  sd_ATE_X_IPTW_rep[i_rep] = sd_ATE_X_IPTW
  sd_ATE_X_TMLE_rep[i_rep] = sd_ATE_X_TMLE
  
  # 95% C.I of ATE_X
  IsCovered_ATE_X_BART_rep[i_rep] = Check_ATE_X_BART
  IsCovered_ATE_X_lm_rep[i_rep] = Check_ATE_X_lm
  IsCovered_ATE_X_IPTW_rep[i_rep] = Check_ATE_X_IPTW 
  IsCovered_ATE_X_TMLE_rep[i_rep] = Check_ATE_X_TMLE
  
}


# ###################################
# # Causal inference for Y
# ###################################
# BART
ATE_Y_BART = mean(ATE_Y_BART_rep)
bias_ATE_Y_BART = mean(ATE_Y_BART_rep -ATE_Y_true)
emp_sd_ATE_Y_BART = sd(ATE_Y_BART_rep)
avg_sd_ATE_Y_BART = mean(sd_ATE_Y_BART_rep)
coverage_ATE_Y_BART= mean(IsCovered_ATE_Y_BART_rep)
RMSE_ATE_Y_BART = sqrt(mean((ATE_Y_BART_rep - ATE_Y_true)^2))
MAE_ATE_Y_BART = median(abs(ATE_Y_BART_rep - ATE_Y_true))

# linear regression
ATE_Y_lm = mean(ATE_Y_lm_rep)
bias_ATE_Y_lm = mean(ATE_Y_lm_rep -ATE_Y_true)
emp_sd_ATE_Y_lm = sd(ATE_Y_lm_rep)
avg_sd_ATE_Y_lm = mean(sd_ATE_Y_lm_rep)
coverage_ATE_Y_lm= mean(IsCovered_ATE_Y_lm_rep)
RMSE_ATE_Y_lm = sqrt(mean((ATE_Y_lm_rep - ATE_Y_true)^2))
MAE_ATE_Y_lm = median(abs(ATE_Y_lm_rep - ATE_Y_true))

# IPTW 
ATE_Y_IPTW = mean(ATE_Y_IPTW_rep)
bias_ATE_Y_IPTW = mean(ATE_Y_IPTW_rep -ATE_Y_true)
emp_sd_ATE_Y_IPTW = sd(ATE_Y_IPTW_rep)
avg_sd_ATE_Y_IPTW = mean(sd_ATE_Y_IPTW_rep)
coverage_ATE_Y_IPTW= mean(IsCovered_ATE_Y_IPTW_rep)
RMSE_ATE_Y_IPTW = sqrt(mean((ATE_Y_IPTW_rep - ATE_Y_true)^2))
MAE_ATE_Y_IPTW = median(abs(ATE_Y_IPTW_rep - ATE_Y_true))

# TMLE
ATE_Y_TMLE = mean(ATE_Y_TMLE_rep)
bias_ATE_Y_TMLE = mean(ATE_Y_TMLE_rep -ATE_Y_true)
emp_sd_ATE_Y_TMLE = sd(ATE_Y_TMLE_rep)
avg_sd_ATE_Y_TMLE = mean(sd_ATE_Y_TMLE_rep)
coverage_ATE_Y_TMLE= mean(IsCovered_ATE_Y_TMLE_rep)
RMSE_ATE_Y_TMLE = sqrt(mean((ATE_Y_TMLE_rep - ATE_Y_true)^2))
MAE_ATE_Y_TMLE = median(abs(ATE_Y_TMLE_rep - ATE_Y_true))


method_causal = c("BART", "LM", "IPTW", "TMLE")
effect_causal = rep("ATE_Y", times= 4)
Est_causal = c(ATE_Y_BART, ATE_Y_lm, ATE_Y_IPTW, ATE_Y_TMLE)
bias_causal = c(bias_ATE_Y_BART, bias_ATE_Y_lm, bias_ATE_Y_IPTW, bias_ATE_Y_TMLE)
emp_sd_causal = c(emp_sd_ATE_Y_BART, emp_sd_ATE_Y_lm, emp_sd_ATE_Y_IPTW, emp_sd_ATE_Y_TMLE)
avg_sd_causal = c(avg_sd_ATE_Y_BART, avg_sd_ATE_Y_lm, avg_sd_ATE_Y_IPTW, avg_sd_ATE_Y_TMLE)
coverage_causal = c(coverage_ATE_Y_BART, coverage_ATE_Y_lm, coverage_ATE_Y_IPTW, coverage_ATE_Y_TMLE)
RMSE_causal = c(RMSE_ATE_Y_BART, RMSE_ATE_Y_lm, RMSE_ATE_Y_IPTW, RMSE_ATE_Y_TMLE)
MAE_causal = c(MAE_ATE_Y_BART, MAE_ATE_Y_lm, MAE_ATE_Y_IPTW, MAE_ATE_Y_TMLE)


method_causal = c("BART", "LM", "IPTW", "TMLE")
effect_causal = rep("ATE_Y", times= 4)
bias_causal = c(bias_ATE_Y_BART, bias_ATE_Y_lm, bias_ATE_Y_IPTW, bias_ATE_Y_TMLE)
emp_sd_causal = c(emp_sd_ATE_Y_BART, emp_sd_ATE_Y_lm, emp_sd_ATE_Y_IPTW, emp_sd_ATE_Y_TMLE)
avg_sd_causal = c(avg_sd_ATE_Y_BART, avg_sd_ATE_Y_lm, avg_sd_ATE_Y_IPTW, avg_sd_ATE_Y_TMLE)
coverage_causal = c(coverage_ATE_Y_BART, coverage_ATE_Y_lm, coverage_ATE_Y_IPTW, coverage_ATE_Y_TMLE)
RMSE_causal = c(RMSE_ATE_Y_BART, RMSE_ATE_Y_lm, RMSE_ATE_Y_IPTW, RMSE_ATE_Y_TMLE)
MAE_causal = c(MAE_ATE_Y_BART, MAE_ATE_Y_lm, MAE_ATE_Y_IPTW, MAE_ATE_Y_TMLE)

causal_res_Y = data.frame(method = method_causal, effect= effect_causal, bias= bias_causal, 
                          emp_sd= emp_sd_causal, avg_sd = avg_sd_causal, coverage = coverage_causal, 
                          RMSE = RMSE_causal, MAE = MAE_causal)


write.csv(causal_res_Y, file = "6_Off_Shelf_Est_ATE_Y.csv")


# ###################################
# # Causal inference for X
# ###################################
# BART
ATE_X_BART = mean(ATE_X_BART_rep)
bias_ATE_X_BART = mean(ATE_X_BART_rep -ATE_X_true)
emp_sd_ATE_X_BART = sd(ATE_X_BART_rep)
avg_sd_ATE_X_BART = mean(sd_ATE_X_BART_rep)
coverage_ATE_X_BART= mean(IsCovered_ATE_X_BART_rep)
RMSE_ATE_X_BART = sqrt(mean((ATE_X_BART_rep - ATE_X_true)^2))
MAE_ATE_X_BART = median(abs(ATE_X_BART_rep - ATE_X_true))

# linear regression
ATE_X_lm = mean(ATE_X_lm_rep)
bias_ATE_X_lm = mean(ATE_X_lm_rep -ATE_X_true)
emp_sd_ATE_X_lm = sd(ATE_X_lm_rep)
avg_sd_ATE_X_lm = mean(sd_ATE_X_lm_rep)
coverage_ATE_X_lm= mean(IsCovered_ATE_X_lm_rep)
RMSE_ATE_X_lm = sqrt(mean((ATE_X_lm_rep - ATE_X_true)^2))
MAE_ATE_X_lm = median(abs(ATE_X_lm_rep - ATE_X_true))

# IPTW 
ATE_X_IPTW = mean(ATE_X_IPTW_rep)
bias_ATE_X_IPTW = mean(ATE_X_IPTW_rep -ATE_X_true)
emp_sd_ATE_X_IPTW = sd(ATE_X_IPTW_rep)
avg_sd_ATE_X_IPTW = mean(sd_ATE_X_IPTW_rep)
coverage_ATE_X_IPTW= mean(IsCovered_ATE_X_IPTW_rep)
RMSE_ATE_X_IPTW = sqrt(mean((ATE_X_IPTW_rep - ATE_X_true)^2))
MAE_ATE_X_IPTW = median(abs(ATE_X_IPTW_rep - ATE_X_true))

# TMLE
ATE_X_TMLE = mean(ATE_X_TMLE_rep)
bias_ATE_X_TMLE = mean(ATE_X_TMLE_rep -ATE_X_true)
emp_sd_ATE_X_TMLE = sd(ATE_X_TMLE_rep)
avg_sd_ATE_X_TMLE = mean(sd_ATE_X_TMLE_rep)
coverage_ATE_X_TMLE= mean(IsCovered_ATE_X_TMLE_rep)
RMSE_ATE_X_TMLE = sqrt(mean((ATE_X_TMLE_rep - ATE_X_true)^2))
MAE_ATE_X_TMLE = median(abs(ATE_X_TMLE_rep - ATE_X_true))

method_causal = c("BART", "LM", "IPTW", "TMLE")
effect_causal = rep("ATE_X", times= 4)
Est_causal = c(ATE_X_BART, ATE_X_lm, ATE_X_IPTW, ATE_X_TMLE)
bias_causal = c(bias_ATE_X_BART, bias_ATE_X_lm, bias_ATE_X_IPTW, bias_ATE_X_TMLE)
emp_sd_causal = c(emp_sd_ATE_X_BART, emp_sd_ATE_X_lm, emp_sd_ATE_X_IPTW, emp_sd_ATE_X_TMLE)
avg_sd_causal = c(avg_sd_ATE_X_BART, avg_sd_ATE_X_lm, avg_sd_ATE_X_IPTW, avg_sd_ATE_X_TMLE)
coverage_causal = c(coverage_ATE_X_BART, coverage_ATE_X_lm, coverage_ATE_X_IPTW, coverage_ATE_X_TMLE)
RMSE_causal = c(RMSE_ATE_X_BART, RMSE_ATE_X_lm, RMSE_ATE_X_IPTW, RMSE_ATE_X_TMLE)
MAE_causal = c(MAE_ATE_X_BART, MAE_ATE_X_lm, MAE_ATE_X_IPTW, MAE_ATE_X_TMLE)

causal_res_X = data.frame(method = method_causal, effect= effect_causal, bias= bias_causal, 
                          emp_sd= emp_sd_causal, avg_sd = avg_sd_causal, coverage = coverage_causal, 
                          RMSE = RMSE_causal, MAE = MAE_causal)

write.csv(causal_res_X, file = "6_Off_Shelf_Est_ATE_X.csv")


