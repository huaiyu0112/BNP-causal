rm(list=ls())  

library(mice) # MICE imputation 
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

rep_seq = 1: 1000
n_seq = length(rep_seq)
n_sample = dim(IncompleteDATA[[1]]$L1_7)[[1]]
ATE_Y_true= 1.5 
ATE_X_true = 0.5

dir.create("4_MICE_existing_method_est_MAR")

for (i_rep in rep_seq){
	
  # i_rep = 1
  print(paste0("i_rep=",i_rep))
  seed_i_rep = 100 + i_rep
	set.seed(seed_i_rep)
	  
  Y_missing = IncompleteDATA[[i_rep]]$Outcome[, 1] 
  X_missing = IncompleteDATA[[i_rep]]$Outcome[, 2] 
  A = IncompleteDATA[[i_rep]]$TrueA
  Incomplete_Y = IncompleteDATA[[i_rep]]$L1_7[,6:7]
  Incomplete_X = IncompleteDATA[[i_rep]]$L1_7[,1:5]
  L1_7_missing = cbind(Incomplete_X, Incomplete_Y)
  
  # imputation model that include A
  m_Imp = 10
  # may change the order 
  dat_sim2_missing = cbind(L1_7_missing, A, X_missing, Y_missing)
  dat_sim2_imp_mice = mice(dat_sim2_missing, m=m_Imp ,maxit=50,meth='pmm',seed= seed_i_rep)
  
  # save est in each completed data
  ATE_Y_lm_imp = ATE_Y_BART_imp =  ATE_Y_IPTW_imp =  ATE_Y_TMLE_imp =  numeric(m_Imp)
  ATE_X_lm_imp = ATE_X_BART_imp =  ATE_X_IPTW_imp =  ATE_X_TMLE_imp =  numeric(m_Imp)
  sd_ATE_Y_lm_imp = sd_ATE_Y_BART_imp =  sd_ATE_Y_IPTW_imp =  sd_ATE_Y_TMLE_imp =  numeric(m_Imp)
  sd_ATE_X_lm_imp = sd_ATE_X_BART_imp =  sd_ATE_X_IPTW_imp =  sd_ATE_X_TMLE_imp =  numeric(m_Imp)
  
  for (i_imp in 1: 10){
    
    # i_imp = 1
    dat_sim2_imp = complete(dat_sim2_imp_mice, i_imp)
    L1_7_after_imp = dat_sim2_imp[, 1:7]
    
    ########################
    ## Linear regression  ##
    ########################
    lm_Y_after_imp = lm(Y_missing ~ factor(A) + factor(L1) + factor(L2) + factor(L3) + factor(L4)+ factor(L5) + L6 + L7, data = dat_sim2_imp )
    ATE_Y_lm_imp[i_imp] = summary(lm_Y_after_imp)$coefficients[2,1]
    sd_ATE_Y_lm_imp[i_imp]  = summary(lm_Y_after_imp)$coefficients[2,2]
    
    lm_X_after_imp = lm(X_missing ~ factor(A) + factor(L1) + factor(L2) + factor(L3) + factor(L4)+ factor(L5) + L6 + L7, data = dat_sim2_imp )
    ATE_X_lm_imp[i_imp] = summary(lm_X_after_imp)$coefficients[2,1]
    sd_ATE_X_lm_imp[i_imp]  = summary(lm_X_after_imp)$coefficients[2,2]
    
    ###########
    ## BART  ##
    ###########
    bartc_ATE_Y_after_imp = bartc(response= Y_missing, treatment= A, confounders= factor(L1) + factor(L2) + factor(L3) + factor(L4)+ factor(L5) + L6 + L7, 
                                  estimand   = c("ate"), data = dat_sim2_imp)
    ATE_Y_BART_imp[i_imp] = fitted(bartc_ATE_Y_after_imp)
    sd_ATE_Y_BART_imp[i_imp] = sd(extract(bartc_ATE_Y_after_imp ))
    
    bartc_ATE_X_after_imp = bartc(response= X_missing, treatment= A, confounders= factor(L1) + factor(L2) + factor(L3) + factor(L4)+ factor(L5) + L6 + L7, 
                                  estimand   = c("ate"), data = dat_sim2_imp)
    ATE_X_BART_imp[i_imp] = fitted(bartc_ATE_X_after_imp)
    sd_ATE_X_BART_imp[i_imp] = sd(extract(bartc_ATE_X_after_imp ))
    
    ###########
    ## IPTW ##
    ##########
    
    weightmodel_ATE= weightit(A ~ factor(L1) + factor(L2) + factor(L3) + factor(L4)+ factor(L5) + L6 + L7, data = dat_sim2_imp, method = "ps", estimand = "ATE")
    dat_sim2_imp$wt_ATE = get.w(weightmodel_ATE)
    msm_ATE_Y = (svyglm(Y_missing ~ A, design = svydesign(~ 1, weights = ~ wt_ATE, data =dat_sim2_imp)))
    ATE_Y_IPTW_imp[i_imp] = coef(msm_ATE_Y)[2]
    sd_ATE_Y_IPTW_imp[i_imp] = summary(msm_ATE_Y)$"coefficients"[2,2]
    
    msm_ATE_X = (svyglm(X_missing ~ A, design = svydesign(~ 1, weights = ~ wt_ATE, data =dat_sim2_imp)))
    ATE_X_IPTW_imp[i_imp] = coef(msm_ATE_X)[2]
    sd_ATE_X_IPTW_imp[i_imp] = summary(msm_ATE_X)$"coefficients"[2,2]
  
    ###########
    ## TMLE  ##
    ###########
    # R-package tmle (base implementation includes SL.step, SL.glm and SL.glm.interaction)
  
    res.tmle_Y = tmle(Y = dat_sim2_imp$Y_missing, A = dat_sim2_imp$A, W= dat_sim2_imp[, c("L1", "L2", "L3", "L4", "L5", "L6", "L7")],
                      Qform= "Y ~ A +  factor(L1) + factor(L2) + factor(L3) + factor(L4)+ factor(L5) + L6 + L7",
                      gform = "A ~ factor(L1) + factor(L2) + factor(L3) + factor(L4)+ factor(L5) + L6 + L7",
                      # g.SL.library = c("SL.mean", "SL.glm", "SL.step"),
                      # Q.SL.library = c("SL.mean", "SL.glm", "SL.step"),
                      fluctuation = "logistic")
    
    ATE_Y_TMLE_imp[i_imp] = res.tmle_Y$estimates$ATE$psi
    sd_ATE_Y_TMLE_imp[i_imp] = sqrt(res.tmle_Y$estimates$ATE$var.psi)
    
    res.tmle_X = tmle(Y = dat_sim2_imp$X_missing, A = dat_sim2_imp$A, W= dat_sim2_imp[, c("L1", "L2", "L3", "L4", "L5", "L6", "L7")],
                      Qform= "Y ~ A +  factor(L1) + factor(L2) + factor(L3) + factor(L4)+ factor(L5) + L6 + L7",
                      gform = "A ~ factor(L1) + factor(L2) + factor(L3) + factor(L4)+ factor(L5) + L6 + L7",
                      # g.SL.library = c("SL.mean", "SL.glm", "SL.step"),
                      # Q.SL.library = c("SL.mean", "SL.glm", "SL.step"),
                      fluctuation = "logistic")
    
    ATE_X_TMLE_imp[i_imp] = res.tmle_X$estimates$ATE$psi
    sd_ATE_X_TMLE_imp[i_imp] = sqrt(res.tmle_X$estimates$ATE$var.psi)
    
  }
  
  # Combining Rule
  mi_mean_output = cbind(ATE_Y_lm_imp , ATE_Y_BART_imp ,  ATE_Y_IPTW_imp ,  ATE_Y_TMLE_imp,
                         ATE_X_lm_imp , ATE_X_BART_imp ,  ATE_X_IPTW_imp ,  ATE_X_TMLE_imp)
  mi_se_output = cbind(sd_ATE_Y_lm_imp , sd_ATE_Y_BART_imp ,  sd_ATE_Y_IPTW_imp ,  sd_ATE_Y_TMLE_imp,
                       sd_ATE_X_lm_imp , sd_ATE_X_BART_imp ,  sd_ATE_X_IPTW_imp ,  sd_ATE_X_TMLE_imp)
  
  mu.list <- vector("list", m_Imp)
  ses.list <- vector("list", m_Imp)
  
  for (imp_rep in 1 : m_Imp){
    mu.list[[imp_rep]] = mi_mean_output[imp_rep, ]
    ses.list[[imp_rep]] = mi_se_output[imp_rep, ]
  }
  
  combined.imp.results= mi.inference(est= mu.list, std.err=ses.list, confidence=0.95)
  q.mi = combined.imp.results$`est`
  se.mi = combined.imp.results$std.err
  lower.mi = combined.imp.results$lower
  upper.mi = combined.imp.results$upper
  true.mi = c(rep(ATE_Y_true, times= 4), rep(ATE_X_true, times= 4))
  
  check.mi = numeric(length(q.mi))
  for (i_mi in 1: length(q.mi)) {
    check.mi[i_mi] = (lower.mi[i_mi] < true.mi[i_mi] & true.mi[i_mi] < upper.mi[i_mi])
  }
  
  names(check.mi) = names(q.mi)

  save(true.mi, q.mi, se.mi, check.mi, file = paste("4_MICE_existing_method_est_MAR/res_", i_rep , "_mice_causal_multiimp",".RData", sep = ""))
} # for (i_rep)


rm(list=ls())
load("1_DATA.RData")

rep_seq= 1:1000
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
  load(paste("4_MICE_existing_method_est_MAR/res_", i_rep, "_mice_causal_multiimp",".RData", sep = ""))
  
  ATE_Y_lm_rep[i_rep] = q.mi["ATE_Y_lm_imp"]
  ATE_Y_BART_rep[i_rep] = q.mi["ATE_Y_BART_imp"]
  ATE_Y_IPTW_rep[i_rep] = q.mi["ATE_Y_IPTW_imp"]
  ATE_Y_TMLE_rep[i_rep] = q.mi["ATE_Y_TMLE_imp"]
  
  sd_ATE_Y_lm_rep[i_rep] = se.mi["ATE_Y_lm_imp"]
  sd_ATE_Y_BART_rep[i_rep] = se.mi["ATE_Y_BART_imp"]
  sd_ATE_Y_IPTW_rep[i_rep] = se.mi["ATE_Y_IPTW_imp"]
  sd_ATE_Y_TMLE_rep[i_rep] = se.mi["ATE_Y_TMLE_imp"]
  
  names(check.mi) = names(q.mi)
  IsCovered_ATE_Y_lm_rep[i_rep] = check.mi["ATE_Y_lm_imp"]
  IsCovered_ATE_Y_BART_rep[i_rep] = check.mi["ATE_Y_BART_imp"]
  IsCovered_ATE_Y_IPTW_rep[i_rep] = check.mi["ATE_Y_IPTW_imp"]
  IsCovered_ATE_Y_TMLE_rep[i_rep] = check.mi["ATE_Y_TMLE_imp"]
  
  ATE_X_lm_rep[i_rep] = q.mi["ATE_X_lm_imp"]
  ATE_X_BART_rep[i_rep] = q.mi["ATE_X_BART_imp"]
  ATE_X_IPTW_rep[i_rep] = q.mi["ATE_X_IPTW_imp"]
  ATE_X_TMLE_rep[i_rep] = q.mi["ATE_X_TMLE_imp"]
  
  sd_ATE_X_lm_rep[i_rep] = se.mi["ATE_X_lm_imp"]
  sd_ATE_X_BART_rep[i_rep] = se.mi["ATE_X_BART_imp"]
  sd_ATE_X_IPTW_rep[i_rep] = se.mi["ATE_X_IPTW_imp"]
  sd_ATE_X_TMLE_rep[i_rep] = se.mi["ATE_X_TMLE_imp"]
  
  IsCovered_ATE_X_lm_rep[i_rep] = check.mi["ATE_X_lm_imp"]
  IsCovered_ATE_X_BART_rep[i_rep] = check.mi["ATE_X_BART_imp"]
  IsCovered_ATE_X_IPTW_rep[i_rep] = check.mi["ATE_X_IPTW_imp"]
  IsCovered_ATE_X_TMLE_rep[i_rep] = check.mi["ATE_X_TMLE_imp"]
  
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
effect_causal = rep("ATE_Y", times= 4)
Est_causal = c(ATE_Y_BART, ATE_Y_lm, ATE_Y_IPTW, ATE_Y_TMLE)
bias_causal = c(bias_ATE_Y_BART, bias_ATE_Y_lm, bias_ATE_Y_IPTW, bias_ATE_Y_TMLE)
emp_sd_causal = c(emp_sd_ATE_Y_BART, emp_sd_ATE_Y_lm, emp_sd_ATE_Y_IPTW, emp_sd_ATE_Y_TMLE)
avg_sd_causal = c(avg_sd_ATE_Y_BART, avg_sd_ATE_Y_lm, avg_sd_ATE_Y_IPTW, avg_sd_ATE_Y_TMLE)
coverage_causal = c(coverage_ATE_Y_BART, coverage_ATE_Y_lm, coverage_ATE_Y_IPTW, coverage_ATE_Y_TMLE)
RMSE_causal = c(RMSE_ATE_Y_BART, RMSE_ATE_Y_lm, RMSE_ATE_Y_IPTW, RMSE_ATE_Y_TMLE)
MAE_causal = c(MAE_ATE_Y_BART, MAE_ATE_Y_lm, MAE_ATE_Y_IPTW, MAE_ATE_Y_TMLE)

causal_res_Y = data.frame(method = method_causal, effect= effect_causal, Est= Est_causal, 
                          emp_sd= emp_sd_causal, RMSE = RMSE_causal, MAE = MAE_causal, coverage = coverage_causal, avg_sd = avg_sd_causal)
write.csv(causal_res_Y, file = "4_MICE_existing_method_est_ATE_Y.csv")

method_causal = c("BART", "LM", "IPTW", "TMLE")
effect_causal = rep("ATE_X", times= 4)
Est_causal = c(ATE_X_BART, ATE_X_lm, ATE_X_IPTW, ATE_X_TMLE)
bias_causal = c(bias_ATE_X_BART, bias_ATE_X_lm, bias_ATE_X_IPTW, bias_ATE_X_TMLE)
emp_sd_causal = c(emp_sd_ATE_X_BART, emp_sd_ATE_X_lm, emp_sd_ATE_X_IPTW, emp_sd_ATE_X_TMLE)
avg_sd_causal = c(avg_sd_ATE_X_BART, avg_sd_ATE_X_lm, avg_sd_ATE_X_IPTW, avg_sd_ATE_X_TMLE)
coverage_causal = c(coverage_ATE_X_BART, coverage_ATE_X_lm, coverage_ATE_X_IPTW, coverage_ATE_X_TMLE)
RMSE_causal = c(RMSE_ATE_X_BART, RMSE_ATE_X_lm, RMSE_ATE_X_IPTW, RMSE_ATE_X_TMLE)
MAE_causal = c(MAE_ATE_X_BART, MAE_ATE_X_lm, MAE_ATE_X_IPTW, MAE_ATE_X_TMLE)

causal_res_X = data.frame(method = method_causal, effect= effect_causal, Est= Est_causal, 
                           emp_sd= emp_sd_causal, RMSE = RMSE_causal, MAE = MAE_causal, coverage = coverage_causal, avg_sd = avg_sd_causal)
write.csv(causal_res_X, file = "4_MICE_existing_method_est_ATE_X.csv")
