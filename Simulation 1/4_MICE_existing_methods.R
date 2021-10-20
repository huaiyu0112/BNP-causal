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

# True causal effect 
ATE_true = 1.5

expit_fn = function(x) { 1/(1+exp(-x)) } # inverse of the logit function

# How many times the simulation are repeated 
start_seq =0
rep_seq = c(1:400) #1:40 # 1:40
n_seq = length(rep_seq)

load("1_DATA_BD.RData")
load("3_DATA_MAR.RData")

dir.create("5_MICE_Off_Shelf_MAR")

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
  Incomplete_Y = incompleteDATA[[i_rep]]$L1_6[,1:4]
  Incomplete_X = incompleteDATA[[i_rep]]$L1_6[,5:6]
  L1_6_missing = cbind(Incomplete_Y, Incomplete_X)
  L1_6 = DATA[[i_rep]]$L1_6[,1:6]
  
  # imputation model includes A and Y
  m_Imp = 10
  dat_sim1_missing = cbind(L1_6_missing, A, Y_obs)
  dat_sim1_imp_mice = mice(dat_sim1_missing, m=m_Imp ,maxit=50,meth='pmm',seed=500)
  
  # save est in each completed data
  mu_L2_out  = mu_L4_out = numeric(m_Imp)
  se_L2_out  = se_L4_out = numeric(m_Imp)
  ATE_lm_imp = ATE_BART_imp =  ATE_IPTW_imp =  ATE_TMLE_imp =  numeric(m_Imp)
  sd_ATE_lm_imp = sd_ATE_BART_imp = sd_ATE_IPTW_imp = sd_ATE_TMLE_imp = numeric(m_Imp)
  
  for (i_imp in 1: 10){
    
    # i_imp = 1
    dat_sim1_imp = complete(dat_sim1_imp_mice, i_imp)
    L1_6_after_imp = dat_sim1_imp[, 1:6]
    
    mu_L2_out[i_imp] = mean(L1_6_after_imp[, "L2"])
    se_L2_out[i_imp] = sqrt(var(L1_6_after_imp[, "L2"])/n_sample)
    mu_L4_out[i_imp] = mean(L1_6_after_imp[, "L4"])
    se_L4_out[i_imp] = sqrt(var(L1_6_after_imp[, "L4"])/n_sample)
    
    ########################
    ## Linear regression  ##
    ########################
    lm_after_imp = lm(Y_obs ~ as.factor(A) + L1 + L2 + L3 + L4+ as.factor(L5) + as.factor(L6), data = dat_sim1_imp)
    ATE_lm_imp[i_imp] = summary(lm_after_imp)$coefficients[2,1]
    sd_ATE_lm_imp[i_imp]  = summary(lm_after_imp)$coefficients[2,2]
    
    ###########
    ## BART  ##
    ###########
    bartc_ATE_after_imp = bartc(response= Y_obs, treatment= A, confounders= L1 + L2 + L3 + L4 + factor(L5) + factor(L6), estimand   = c("ate"), data = dat_sim1_imp)
    ATE_BART_imp[i_imp] = fitted(bartc_ATE_after_imp)
    sd_ATE_BART_imp[i_imp] = sd(extract(bartc_ATE_after_imp))
    
    ###########
    ## IPTW ##
    ##########
    weightmodel_ATE= weightit(A ~ L1 + L2+ L3 + L4+ + as.factor(L5) + as.factor(L6), data = dat_sim1_imp, method = "ps", estimand = "ATE")
    dat_sim1_imp$wt_ATE = get.w(weightmodel_ATE)
    msm_ATE = (svyglm(Y_obs ~ A, design = svydesign(~ 1, weights = ~ wt_ATE, data =dat_sim1_imp)))
    ATE_IPTW_imp[i_imp] = coef(msm_ATE)[2]
    sd_ATE_IPTW_imp[i_imp] = summary(msm_ATE)$"coefficients"[2,2]
    
    ###########
    ## TMLE  ##
    ###########
    res.tmle = tmle(Y = dat_sim1_imp$Y_obs, A = dat_sim1_imp$A, W= dat_sim1_imp[, c("L1", "L2", "L3", "L4", "L5", "L6")],
                    Qform= "Y ~ A +  L1 + L2 + L3 + L4 + factor(L5) + factor(L6)",
                    gform = "A ~ L1 + L2 + L3 + L4 + factor(L5) + factor(L6)",
                    fluctuation = "logistic")
    ATE_TMLE_imp[i_imp] = res.tmle$estimates$ATE$psi
    sd_ATE_TMLE_imp[i_imp] = sqrt(res.tmle$estimates$ATE$var.psi)
    
  }
  
  # Combining Rule
  mi_mean_output = cbind(mu_L2_out, mu_L4_out, ATE_lm_imp, ATE_BART_imp,  ATE_IPTW_imp,  ATE_TMLE_imp )
  mi_se_output = cbind(se_L2_out, se_L4_out, sd_ATE_lm_imp , sd_ATE_BART_imp , sd_ATE_IPTW_imp , sd_ATE_TMLE_imp)
  
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
  true.mi = c(0, 0, rep(1.5, times= 4))
  
  check.mi = numeric(length(q.mi))
  for (i_mi in 1: length(q.mi)) {
    check.mi[i_mi] = (lower.mi[i_mi] < true.mi[i_mi] & true.mi[i_mi] < upper.mi[i_mi])
  }
  
  names(check.mi) = names(q.mi)
    
  save(true.mi, q.mi, se.mi, check.mi, file = paste("5_MICE_Off_Shelf_MAR/res_", i_rep + start_seq, "_mice_causal_multiimp",".RData", sep = ""))
  
}


rm(list=ls()) 
load("1_DATA_BD.RData")
load("3_DATA_MAR.RData")

# True causal effect 
ATE_true = 1.5
rep_seq = c(1:400) #1:40 # 1:40

ATE_BART_rep = ATE_lm_rep = ATE_IPTW_rep = ATE_TMLE_rep = numeric(length(rep_seq))
sd_ATE_BART_rep = sd_ATE_lm_rep = sd_ATE_IPTW_rep = sd_ATE_TMLE_rep = numeric(length(rep_seq))
IsCovered_ATE_BART_rep = IsCovered_ATE_lm_rep = IsCovered_ATE_IPTW_rep = IsCovered_ATE_TMLE_rep = numeric(length(rep_seq))

mu_L2_true  = mu_L4_true = 0
mu_L2_rep  = mu_L4_rep = numeric(length(rep_seq))
se_mu_L2_rep  = se_mu_L4_rep = se_p1_L5_rep= se_p2_L5_rep= se_p3_L5_rep= se_p1_L6_rep= se_p2_L6_rep= numeric(length(rep_seq))
IsCovered_mu_L2_rep  = IsCovered_mu_L4_rep = IsCovered_p1_L5_rep= IsCovered_p2_L5_rep= IsCovered_p3_L5_rep= IsCovered_p1_L6_rep= IsCovered_p2_L6_rep= numeric(length(rep_seq))


for (i_rep in rep_seq) {
  
  # i_rep = 1
  load(paste("5_MICE_Off_Shelf_MAR/res_", i_rep, "_mice_causal_multiimp",".RData", sep = ""))
  
  mu_L2_rep[i_rep] = q.mi[1]
  mu_L4_rep[i_rep] = q.mi[2]
  se_mu_L2_rep[i_rep] = se.mi[1]
  se_mu_L4_rep[i_rep] = se.mi[2]
  IsCovered_mu_L2_rep[i_rep] = check.mi[1]
  IsCovered_mu_L4_rep[i_rep] = check.mi[2]
  
  
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

causal_res = data.frame(method = method_causal, effect= effect_causal, bias= bias_causal, 
                        emp_sd= emp_sd_causal, avg_sd = avg_sd_causal, 
                        RMSE = RMSE_causal, MAE = MAE_causal, coverage = coverage_causal)

write.csv(causal_res, file = "5_MICE_Off_Shelf_MAR_causal.csv")

causal_res_MS = data.frame(method = method_causal, effect= effect_causal, point_est= point_est_causal, 
                        emp_sd= emp_sd_causal, RMSE = RMSE_causal, MAE = MAE_causal, 
                        coverage = coverage_causal, avg_sd = avg_sd_causal)
write.csv(causal_res_MS, file = "5_MICE_Off_Shelf_MAR_causal_MS.csv")

###################################
# Imputation performance
###################################
est_mu_L2 = mean(mu_L2_rep)
bias_mu_L2 = mean(mu_L2_rep - 0)
emp_sd_mu_L2 = sd(mu_L2_rep)
avg_sd_mu_L2 = mean(se_mu_L2_rep)
coverage_mu_L2 = mean(IsCovered_mu_L2_rep)
RMSE_mu_L2 = sqrt(mean((mu_L2_rep - 0)^2))

est_mu_L4 = mean(mu_L4_rep)
bias_mu_L4 = mean(mu_L4_rep - 0)
emp_sd_mu_L4 = sd(mu_L4_rep)
avg_sd_mu_L4 = mean(se_mu_L4_rep)
coverage_mu_L4 = mean(IsCovered_mu_L4_rep)
RMSE_mu_L4 = sqrt(mean((mu_L4_rep - 0)^2))

# save the results in a csv table
effect_imp = c("mu_L2", "mu_L4")
est_imp = c(est_mu_L2, est_mu_L4)
bias_imp = c(bias_mu_L2, bias_mu_L4)
emp_sd_imp = c(emp_sd_mu_L2, emp_sd_mu_L4)
avg_sd_imp = c(avg_sd_mu_L2, avg_sd_mu_L4)
coverage_imp = c(coverage_mu_L2, coverage_mu_L4)
RMSE_imp = c(RMSE_mu_L2, RMSE_mu_L4)

imp_res = data.frame(Effect= effect_imp, point = est_imp, Emp_sd = emp_sd_imp, Avg_sd = avg_sd_imp, RMSE = RMSE_imp, Coverage = coverage_imp)
write.csv(imp_res, file = "5_MICE_Off_Shelf_MAR_imputation.csv")


