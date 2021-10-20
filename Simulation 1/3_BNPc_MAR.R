rm(list=ls())  

library( HCMMcausal )  # for running BNP causal 
library( norm ) # for mi.inference function

load("1_DATA_BD.RData")
load("3_DATA_MAR.RData")

# How many times the simulation are repeated 
rep_seq = 1:400

n_seq = length(rep_seq)
n_sample = dim(incompleteDATA[[1]]$L1_6)[[1]]

n_burnin = 2000 ; m_Imp = 10 ; interval_btw_Imp = 200

ATE_HCMM_causal_rep = SD_ATE_HCMM_causal_rep = isCovered_ATE_HCMM_causal_rep = numeric(n_seq)
dir.create("4_BNPc_MAR")

for (i_rep in rep_seq) {

  print(paste0("i_rep=",i_rep))
	set.seed(100 + i_rep)
	  
  obs_response = incompleteDATA[[i_rep]]$Y_obs 
  TrueA = incompleteDATA[[i_rep]]$TrueA
  Incomplete_Y = incompleteDATA[[i_rep]]$L1_6[,1:4]
  Incomplete_X = incompleteDATA[[i_rep]]$L1_6[,5:6]

  data_obj = readData(Response_var=obs_response, trt_indicator=TrueA, Cont_pred=Incomplete_Y, Categ_pred=Incomplete_X, RandomSeed=100+i_rep)
  model_obj = createModel(data_obj)	
	result_obj = multipleImp(data_obj, model_obj, n_burnin, m_Imp, interval_btw_Imp, show_iter=FALSE)
  
	#################################################
	# Imputed value summary 
	#################################################
	
  MI_X_cube= result_obj$MI_X_cube
  MI_X_cube[,,2] = MI_X_cube[,,2]-1
  MI_Y_cube= result_obj$MI_Y_cube
  
  # Combining rule
  mu_L2_out  = mu_L4_out = numeric(m_Imp)
  se_L2_out  = se_L4_out = numeric(m_Imp)
  
  for (imp_rep in 1:m_Imp){
    
    mu_L2_out[imp_rep] = mean(result_obj$MI_Y_cube[imp_rep,, (1+2)])
    se_L2_out[imp_rep] = sqrt(var(result_obj$MI_Y_cube[imp_rep,, (1+2)])/n_sample)
    mu_L4_out[imp_rep] = mean(result_obj$MI_Y_cube[imp_rep,, (1+4)])
    se_L4_out[imp_rep] = sqrt(var(result_obj$MI_Y_cube[imp_rep,, (1+4)])/n_sample)
    
  }
  
  true.mi = c(0, 0)
  mi_mean_output = cbind(mu_L2_out, mu_L4_out)
  mi_se_output = cbind(se_L2_out, se_L4_out)
  
  mu.list <- vector("list", m_Imp)
  ses.list <- vector("list", m_Imp)
  
  for (imp_rep in 1:m_Imp){
    mu.list[[imp_rep]] = mi_mean_output[imp_rep, ]
    ses.list[[imp_rep]] = mi_se_output[imp_rep, ]
  }
  
  combined.results= mi.inference(est= mu.list, std.err=ses.list, confidence=0.95)
  q.mi = combined.results$`est`
  se.mi = combined.results$std.err
  lower.mi = combined.results$lower
  upper.mi = combined.results$upper
  
  check.mi = numeric(length(q.mi))
  for (i_mi in 1: length(q.mi)) {
    check.mi[i_mi] = (lower.mi[i_mi] < 0 & 0 < upper.mi[i_mi])
  }
  
  ATE_HCMM_causal = mean(result_obj$est_delta)
  SD_ATE_HCMM_causal = sd(result_obj$est_delta)
  Right_ATE_HCMM_causal = ATE_HCMM_causal + qnorm(0.975)*SD_ATE_HCMM_causal
  Left_ATE_HCMM_causal = ATE_HCMM_causal - qnorm(0.975)*SD_ATE_HCMM_causal
  Check_ATE_HCMM_causal = (Left_ATE_HCMM_causal < ATE_true & ATE_true < Right_ATE_HCMM_causal)
  
  save(ATE_HCMM_causal, SD_ATE_HCMM_causal, Check_ATE_HCMM_causal, true.mi, q.mi, se.mi, check.mi, file = paste0("4_BNPc_MAR/Result_", i_rep,".RData"))
			 
} # for (i_rep)


#######################################################
# Result Summary
#######################################################
rm(list=ls())  

load("1_DATA_BD.RData")

ATE_true = 1.5 # true causal effect 
rep_seq = 1: 400
ATE_HCMM_rep = SD_ATE_HCMM_rep = IsCovered_ATE_HCMM_rep = numeric(length(rep_seq))

mu_L2_rep  = mu_L4_rep = numeric(length(rep_seq))
se_mu_L2_rep  = se_mu_L4_rep =  numeric(length(rep_seq))
IsCovered_mu_L2_rep  = IsCovered_mu_L4_rep = IsCovered_p2_L6_rep= numeric(length(rep_seq))

for (i_rep in rep_seq) {
  
  load(paste0("4_BNPc_MAR/Result_", i_rep,".RData"))
  
  ATE_HCMM_rep[i_rep] = ATE_HCMM_causal
  SD_ATE_HCMM_rep[i_rep] = SD_ATE_HCMM_causal
  IsCovered_ATE_HCMM_rep[i_rep] = Check_ATE_HCMM_causal
  
  mu_L2_rep[i_rep] = q.mi[1]
  mu_L4_rep[i_rep] = q.mi[2]
  se_mu_L2_rep[i_rep] = se.mi[1]
  se_mu_L4_rep[i_rep] = se.mi[2]
  IsCovered_mu_L2_rep[i_rep] = check.mi[1]
  IsCovered_mu_L4_rep[i_rep] = check.mi[2]
  
}

###################################
# Causal inference
###################################
method = "BNPc"
effect = "ATE"
avg_ATE_HCMM= mean(ATE_HCMM_rep)
bias_ATE_HCMM = mean(ATE_HCMM_rep -ATE_true)
emp_sd_ATE_HCMM = sd(ATE_HCMM_rep)
avg_sd_ATE_HCMM = mean(SD_ATE_HCMM_rep)
coverage_ATE_HCMM= mean(IsCovered_ATE_HCMM_rep)
RMSE_ATE_HCMM = sqrt(mean((ATE_HCMM_rep - ATE_true)^2))
MAE_ATE_HCMM = median(abs(ATE_HCMM_rep - ATE_true))

HCMM_causal_res_MS = data.frame(Method= method, effect, Point_est = avg_ATE_HCMM, Emp_sd = emp_sd_ATE_HCMM,
                                RMSE = RMSE_ATE_HCMM, MAE= MAE_ATE_HCMM, Coverage= coverage_ATE_HCMM, Avg_sd= avg_sd_ATE_HCMM)
write.csv(HCMM_causal_res_MS, file = "4_BNPc_MAR_causal_MS.csv")


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

effect_imp = c("mu_L2", "mu_L4")
est_imp = c(est_mu_L2, est_mu_L4)
bias_imp = c(bias_mu_L2, bias_mu_L4)
emp_sd_imp = c(emp_sd_mu_L2, emp_sd_mu_L4)
avg_sd_imp = c(avg_sd_mu_L2, avg_sd_mu_L4)
coverage_imp = c(coverage_mu_L2, coverage_mu_L4)
RMSE_imp = c(RMSE_mu_L2, RMSE_mu_L4)

imp_res = data.frame(Method = "BNPc", Est= effect_imp, point= est_imp, Emp_sd = emp_sd_imp, Avg_sd = avg_sd_imp, RMSE = RMSE_imp, Coverage = coverage_imp)
write.csv(imp_res, file = "4_BNPc_MAR_imputation.csv")
