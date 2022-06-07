rm(list=ls())  

library( HCMMcausal )  # for running BNP causal 
library( norm )

load("1_DATA_BD.RData")
load("3_DATA_MAR.RData")

# How many times the simulation are repeated 
rep_seq = 1:1000
n_seq = length(rep_seq)
n_sample = dim(incompleteDATA[[1]]$L1_7)[[1]]

n_burnin = 2000 ; m_Imp = 10 ; interval_btw_Imp = 1000


dir.create("4_BNPc_MAR")

for (i_rep in rep_seq) {

  print(paste0("i_rep=",i_rep))
	set.seed(100 + i_rep)
	  
  obs_response = incompleteDATA[[i_rep]]$Y_obs 
  TrueA = incompleteDATA[[i_rep]]$TrueA
  Incomplete_Y = incompleteDATA[[i_rep]]$L1_7[,6:7]
  Incomplete_X = incompleteDATA[[i_rep]]$L1_7[,1:5]

  data_obj = readData(Response_var=obs_response, trt_indicator=TrueA, Cont_pred=Incomplete_Y, Categ_pred=Incomplete_X, RandomSeed=100+i_rep)
  model_obj = createModel(data_obj)	
	result_obj = multipleImp(data_obj, model_obj, n_burnin, m_Imp, interval_btw_Imp, show_iter=FALSE)
  
	#################################################
	# Imputed value summary 
	#################################################
	
  MI_X_cube= result_obj$MI_X_cube
  MI_X_cube[,,2] = MI_X_cube[,,2]-1
  MI_Y_cube= result_obj$MI_Y_cube
  
  mu_L6_out = se_L6_out = numeric(m_Imp)
  
  for (imp_rep in 1 : m_Imp){
    
    mu_L6_out[imp_rep] = mean(result_obj$MI_Y_cube[imp_rep,, (1+1)])
    se_L6_out[imp_rep] = sqrt(var(result_obj$MI_Y_cube[imp_rep,, (1+1)])/n_sample)
  }
  
  mi_mean_output = cbind(mu_L6_out)
  mi_se_output = cbind(se_L6_out)
  
  mu.list <- vector("list", m_Imp)
  ses.list <- vector("list", m_Imp)
  
  for (imp_rep in 1 : m_Imp){
    mu.list[[imp_rep]] = mi_mean_output[imp_rep, ]
    ses.list[[imp_rep]] = mi_se_output[imp_rep, ]
  }
  
  combined.results= mi.inference(est= mu.list, std.err=ses.list, confidence=0.95)
  q.mi = combined.results$`est`
  se.mi = combined.results$std.err
  lower.mi = combined.results$lower
  upper.mi = combined.results$upper
  true.mi = mean_L6
  
  check.mi = numeric(length(q.mi))
  for (i_mi in 1: length(q.mi)) {
    check.mi[i_mi] = (lower.mi[i_mi] < true.mi[i_mi] & true.mi[i_mi] < upper.mi[i_mi])
  }
  
  names(check.mi) = names(q.mi)
  
  ATE_BNPc = mean(result_obj$est_delta)
  SD_ATE_BNPc = sd(result_obj$est_delta)
  Right_ATE_BNPc = ATE_BNPc + qnorm(0.975)*SD_ATE_BNPc
  Left_ATE_BNPc = ATE_BNPc - qnorm(0.975)*SD_ATE_BNPc
  Check_ATE_BNPc = (Left_ATE_BNPc < ATE_true & ATE_true < Right_ATE_BNPc)
  
  save(ATE_BNPc, SD_ATE_BNPc, Check_ATE_BNPc, q.mi, se.mi, check.mi, file = paste0("4_BNPc_MAR/Result_", i_rep,".RData"))
			 
} # for (i_rep)



#######################################################
# Result Summary
#######################################################
rm(list=ls())  
load("1_DATA_BD.RData")

rep_seq = 1:1000
ATE_BNPc_rep = SD_ATE_BNPc_rep = IsCovered_ATE_BNPc_rep = numeric(length(rep_seq))
mu_L6_rep = se_mu_L6_rep = IsCovered_mu_L6_rep = numeric(length(rep_seq))

for (i_rep in rep_seq) {
  
  load(paste0("4_BNPc_MAR/Result_", i_rep,".RData"))
  
  ATE_BNPc_rep[i_rep] = ATE_BNPc
  SD_ATE_BNPc_rep[i_rep] = SD_ATE_BNPc
  IsCovered_ATE_BNPc_rep[i_rep] = Check_ATE_BNPc
  
  mu_L6_rep[i_rep] = q.mi["mu_L6_out"]
  se_mu_L6_rep[i_rep] = se.mi["mu_L6_out"]
  IsCovered_mu_L6_rep[i_rep] = check.mi["mu_L6_out"]
}

###################################
# Causal effect result
###################################
method = "BNPc"
effect = "ATE"
avg_ATE_BNPc= mean(ATE_BNPc_rep)
bias_ATE_BNPc = mean(ATE_BNPc_rep -ATE_true)
emp_sd_ATE_BNPc = sd(ATE_BNPc_rep)
avg_sd_ATE_BNPc = mean(SD_ATE_BNPc_rep)
coverage_ATE_BNPc= mean(IsCovered_ATE_BNPc_rep)
RMSE_ATE_BNPc = sqrt(mean((ATE_BNPc_rep - ATE_true)^2))
MAE_ATE_BNPc = median(abs(ATE_BNPc_rep - ATE_true))

BNPc_res = data.frame(Method= method, effect, Point_est = avg_ATE_BNPc, Emp_sd = emp_sd_ATE_BNPc,
                                RMSE = RMSE_ATE_BNPc, MAE= MAE_ATE_BNPc, Coverage= coverage_ATE_BNPc, Avg_sd= avg_sd_ATE_BNPc)

write.csv(BNPc_res, file = "4_BNPc_MAR_res.csv")

###################################
# Imputation performance
###################################
effect_imp = "E(L6)"
est_mu_L6 = mean(mu_L6_rep)
bias_mu_L6 = mean(mu_L6_rep - mean_L6)
emp_sd_mu_L6 = sd(mu_L6_rep)
avg_sd_mu_L6 = mean(se_mu_L6_rep)
coverage_mu_L6 = mean(IsCovered_mu_L6_rep)
RMSE_mu_L6 = sqrt(mean((mu_L6_rep - mean_L6)^2))

imp_res = data.frame(Effect= effect_imp, Est= est_mu_L6, Emp_sd = emp_sd_mu_L6, 
                     RMSE = RMSE_mu_L6, Coverage = coverage_mu_L6, Avg_sd = avg_sd_mu_L6)
write.csv(imp_res, file = "4_BNPc_MAR_imp.csv")


