rm(list=ls()) 

library(mice) # MICE imputation 
# BART causal package
library(bartCause) # BART causal 
# IPTW packages 
library(ipw) # for ipwpoint function
library(WeightIt)
library(cobalt)
# library(twang)
library(survey) # for svyglm function
# TMLE packages 
library("SuperLearner") # used cross validation to weigh different prediction algorithms for TMLE
library(tmle) # for tmle function
library(gam) # Generalized Additive Models (GAM) prediction algorithm 
library(glmnet) # Glmnet prediction algorithm
library(randomForest) # Random Forest prediction algorithm
library(norm)

# True causal effect 
ATE_true = 1.5

expit_fn = function(x) { 1/(1+exp(-x)) } # inverse of the logit function

# How many times the simulation are repeated 
start_seq = 0
rep_seq = c(1:1000) #1:40 # 1:40
n_seq = length(rep_seq)

load("1_DATA_BD.RData")
load("3_DATA_MAR.RData")

dir.create("5_MICE_plus_existing_method_MAR")

# Define number of observations for each dataset 
n_sample = 500

for (i_rep in rep_seq) {
  
  # i_rep = 1
  set.seed(100 + i_rep + start_seq)
  seed = 100 + i_rep + start_seq
  print(paste0("i_rep=",i_rep))
  print(paste0("seed=",seed))
  
  Y_obs = incompleteDATA[[i_rep]]$Y_obs 
  A = incompleteDATA[[i_rep]]$TrueA
  Incomplete_Y = incompleteDATA[[i_rep]]$L1_7[,6:7]
  Incomplete_X = incompleteDATA[[i_rep]]$L1_7[,1:5]
  L1_7_missing = cbind(Incomplete_X, Incomplete_Y)
  L1_7= DATA[[i_rep]]$L1_7
  
  # imputation model that include A
  m_Imp = 10
  dat_sim2_missing = cbind(L1_7_missing, A, Y_obs)
  dat_sim2_imp_mice = mice(dat_sim2_missing, m=m_Imp ,maxit=50,meth='pmm',seed=500)
  
  # save est in each completed data
  mu_L6_out = se_L6_out = numeric(m_Imp)
  ATE_lm_imp = ATE_BART_imp =  ATE_IPTW_imp =  ATE_TMLE_imp =  numeric(m_Imp)
  sd_ATE_lm_imp = sd_ATE_BART_imp = sd_ATE_IPTW_imp = sd_ATE_TMLE_imp = numeric(m_Imp)
  
  for (i_imp in 1: 10){
    
    # i_imp = 1
    dat_sim2_imp = complete(dat_sim2_imp_mice, i_imp)
    dat_sim2_imp$L1 = factor(dat_sim2_imp$L1)
    dat_sim2_imp$L2 = factor(dat_sim2_imp$L2)
    dat_sim2_imp$L3 = factor(dat_sim2_imp$L3)
    dat_sim2_imp$L4 = factor(dat_sim2_imp$L4)
    dat_sim2_imp$L5 = factor(dat_sim2_imp$L5)
    
    L1_7_after_imp = dat_sim2_imp[, 1:7]

    mu_L6_out[i_imp] = mean(L1_7_after_imp[, "L6"])
    se_L6_out[i_imp] =  sqrt(var(L1_7_after_imp[, "L6"])/n_sample)
    
    ########################
    ## Linear regression  ##
    ########################
    lm_after_imp = lm(Y_obs ~ as.factor(A) + L1_7_after_imp[,1] + L1_7_after_imp[,2] + L1_7_after_imp[,3] + L1_7_after_imp[,4]+
                        L1_7_after_imp[,5] + L1_7_after_imp[,6] + L1_7_after_imp[,7])
    ATE_lm_imp[i_imp] = summary(lm_after_imp)$coefficients[2,1]
    sd_ATE_lm_imp[i_imp]  = summary(lm_after_imp)$coefficients[2,2]
    
    ###########
    ## BART  ##
    ###########
    bartc_ATE_after_imp = bartc(response= Y_obs, treatment= A, confounders= L1_7_after_imp, estimand   = c("ate"))
    ATE_BART_imp[i_imp] = fitted(bartc_ATE_after_imp)
    sd_ATE_BART_imp[i_imp] = sd(extract(bartc_ATE_after_imp))
    
    ###########
    ## IPTW ##
    ##########
    weightmodel_ATE= weightit(A ~ factor(L1) + factor(L2) + factor(L3) + factor(L4) + factor(L5) 
                              + L6 + L7, data = dat_sim2_imp, method = "ps", estimand = "ATE")
    dat_sim2_imp$wt_ATE = get.w(weightmodel_ATE)
    msm_ATE <- (svyglm(Y_obs ~ A, design = svydesign(~ 1, weights = ~ wt_ATE, data =dat_sim2_imp)))
    ATE_IPTW_imp[i_imp] = coef(msm_ATE)[2]
    sd_ATE_IPTW_imp[i_imp] = summary(msm_ATE)$"coefficients"[2,2]
    
    ###########
    ## TMLE  ##
    ###########
    res.tmle = tmle(Y = Y_obs, A = A, W= L1_7_after_imp,
                    Qform= "Y ~ A + factor(L1) + factor(L2) + factor(L3) + factor(L4) + factor(L5) + L6 + L7",
                    gform = "A ~ factor(L1) + factor(L2) + factor(L3) + factor(L4)+ factor(L5) + L6 + L7",
                    # g.SL.library = c("SL.mean", "SL.glm", "SL.step"),
                    # Q.SL.library = c("SL.mean", "SL.glm", "SL.step"),
                    fluctuation = "logistic")
    ATE_TMLE_imp[i_imp] = res.tmle$estimates$ATE$psi
    sd_ATE_TMLE_imp[i_imp] = sqrt(res.tmle$estimates$ATE$var.psi)
    
  }
  
  # Combining Rule
  mi_mean_output = cbind(mu_L6_out, ATE_lm_imp, ATE_BART_imp,  ATE_IPTW_imp,  ATE_TMLE_imp)
  mi_se_output = cbind(se_L6_out, sd_ATE_lm_imp , sd_ATE_BART_imp , sd_ATE_IPTW_imp , sd_ATE_TMLE_imp)
  
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
  
  true.mi = c(mean_L6, rep(1.5, times= 4))
  
  check.mi = numeric(length(q.mi))
  for (i_mi in 1: length(q.mi)) {
    check.mi[i_mi] = (lower.mi[i_mi] < true.mi[i_mi] & true.mi[i_mi] < upper.mi[i_mi])
  }
  
  names(check.mi) = names(q.mi)
    
  save(true.mi, q.mi, se.mi, check.mi, file = paste("5_MICE_plus_existing_method_MAR/res_", i_rep + start_seq, "_mice_causal_multiimp",".RData", sep = ""))
  
}


rm(list=ls()) 
load("1_DATA_BD.RData")
load("3_DATA_MAR.RData")

rep_seq = c(1:1000) #1:40 # 1:40

ATE_BART_rep = ATE_lm_rep = ATE_IPTW_rep = ATE_TMLE_rep = numeric(length(rep_seq))
sd_ATE_BART_rep = sd_ATE_lm_rep = sd_ATE_IPTW_rep = sd_ATE_TMLE_rep = numeric(length(rep_seq))
IsCovered_ATE_BART_rep = IsCovered_ATE_lm_rep = IsCovered_ATE_IPTW_rep = IsCovered_ATE_TMLE_rep = numeric(length(rep_seq))
mu_L6_rep = se_mu_L6_rep = IsCovered_mu_L6_rep = numeric(length(rep_seq))
# mean_L6= 1.923266

for (i_rep in rep_seq) {
  
  load(paste("5_MICE_plus_existing_method_MAR/res_", i_rep, "_mice_causal_multiimp",".RData", sep = ""))
  
  mu_L6_rep[i_rep] = q.mi["mu_L6_out"]
  se_mu_L6_rep[i_rep] = se.mi["mu_L6_out"]
  IsCovered_mu_L6_rep[i_rep] = check.mi["mu_L6_out"]
  
  
  ATE_lm_rep[i_rep] = q.mi["ATE_lm_imp"]
  ATE_BART_rep[i_rep] = q.mi["ATE_BART_imp"]
  ATE_IPTW_rep[i_rep] = q.mi["ATE_IPTW_imp"]
  ATE_TMLE_rep[i_rep] = q.mi["ATE_TMLE_imp"]
  
  sd_ATE_lm_rep[i_rep] = se.mi["ATE_lm_imp"]
  sd_ATE_BART_rep[i_rep] = se.mi["ATE_BART_imp"]
  sd_ATE_IPTW_rep[i_rep] = se.mi["ATE_IPTW_imp"]
  sd_ATE_TMLE_rep[i_rep] = se.mi["ATE_TMLE_imp"]
  
  names(check.mi) = names(q.mi)
  IsCovered_ATE_lm_rep[i_rep] = check.mi["ATE_lm_imp"]
  IsCovered_ATE_BART_rep[i_rep] = check.mi["ATE_BART_imp"]
  IsCovered_ATE_IPTW_rep[i_rep] = check.mi["ATE_IPTW_imp"]
  IsCovered_ATE_TMLE_rep[i_rep] = check.mi["ATE_TMLE_imp"]
  
}

# ###################################
# # Causal inference 
# ###################################
# BART
avg_ATE_BART = mean(ATE_BART_rep)
bias_ATE_BART = mean(ATE_BART_rep -ATE_true)
emp_sd_ATE_BART = sd(ATE_BART_rep)
avg_sd_ATE_BART = mean(sd_ATE_BART_rep)
coverage_ATE_BART= mean(IsCovered_ATE_BART_rep)
RMSE_ATE_BART = sqrt(mean((ATE_BART_rep - ATE_true)^2))
MAE_ATE_BART = median(abs(ATE_BART_rep - ATE_true))

# LM
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

method_causal = c("BART", "LM", "IPTW", "TMLE")
effect_causal = rep("ATE", times= 4)
point_est_causal = c(avg_ATE_BART, avg_ATE_lm, avg_ATE_IPTW, avg_ATE_TMLE)
bias_causal = c(bias_ATE_BART, bias_ATE_lm, bias_ATE_IPTW, bias_ATE_TMLE)
emp_sd_causal = c(emp_sd_ATE_BART, emp_sd_ATE_lm, emp_sd_ATE_IPTW, emp_sd_ATE_TMLE)
avg_sd_causal = c(avg_sd_ATE_BART, avg_sd_ATE_lm, avg_sd_ATE_IPTW, avg_sd_ATE_TMLE)
coverage_causal = c(coverage_ATE_BART, coverage_ATE_lm, coverage_ATE_IPTW, coverage_ATE_TMLE)
RMSE_causal = c(RMSE_ATE_BART, RMSE_ATE_lm, RMSE_ATE_IPTW, RMSE_ATE_TMLE)
MAE_causal = c(MAE_ATE_BART, MAE_ATE_lm, MAE_ATE_IPTW, MAE_ATE_TMLE)

causal_res = data.frame(method = method_causal, effect= effect_causal, point_est= point_est_causal, 
                           emp_sd= emp_sd_causal, RMSE = RMSE_causal, MAE = MAE_causal, 
                           coverage = coverage_causal, avg_sd = avg_sd_causal)
write.csv(causal_res, file = "5_MICE_plus_existing_method_MAR_causal.csv")

###################################
# Imputation performance
###################################
bias_mu_L6 = mean(mu_L6_rep - mean_L6)
est_mu_L6 = mean(mu_L6_rep)
emp_sd_mu_L6 = sd(mu_L6_rep)
avg_sd_mu_L6 = mean(se_mu_L6_rep)
coverage_mu_L6 = mean(IsCovered_mu_L6_rep)
RMSE_mu_L6 = sqrt(mean((mu_L6_rep - mean_L6)^2))

imp_res = data.frame(effect= "E(L6)", Est= est_mu_L6, emp_sd = emp_sd_mu_L6, 
                     RMSE = RMSE_mu_L6, coverage = coverage_mu_L6, avg_sd = avg_sd_mu_L6)
write.csv(imp_res, file = "5_MICE_plus_existing_method_MAR_imputation.csv")

