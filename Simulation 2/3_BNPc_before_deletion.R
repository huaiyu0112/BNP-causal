rm(list=ls())  

library( HCMMcausal )  # for running HCMM-LD causal 
# library(mice) # for function pool in mice (Rubin's combining rule 1987)

load("1_DATA_BD.RData")
dir.create("2_BNPc_BD")

# How many times the simulation are repeated 
rep_seq = 1:400
n_seq = length(rep_seq)
n_burnin = 2000 ; m_Imp = 10 ; interval_btw_Imp = 200

for (i_rep in rep_seq) {

  print(paste0("i_rep=",i_rep))
	set.seed(100 + i_rep)
	  
  obs_response = DATA[[i_rep]]$Y_obs 
  TrueA = DATA[[i_rep]]$TrueA
  Incomplete_Y = DATA[[i_rep]]$L1_7[,6:7]
  Incomplete_X = DATA[[i_rep]]$L1_7[,1:5]

  data_obj = readData(Response_var=obs_response, trt_indicator=TrueA, Cont_pred=Incomplete_Y, Categ_pred=Incomplete_X, RandomSeed=100+i_rep)
  model_obj = createModel(data_obj)
	result_obj = multipleImp(data_obj, model_obj, n_burnin, m_Imp, interval_btw_Imp, show_iter=FALSE)
  	
	ATE_HCMM_causal = mean(result_obj$est_delta)
	SD_ATE_HCMM_causal = sd(result_obj$est_delta)
	Right_ATE_HCMM_causal = ATE_HCMM_causal + qnorm(0.975)*SD_ATE_HCMM_causal
	Left_ATE_HCMM_causal = ATE_HCMM_causal - qnorm(0.975)*SD_ATE_HCMM_causal
	Check_ATE_HCMM_causal = (Left_ATE_HCMM_causal < ATE_true & ATE_true < Right_ATE_HCMM_causal)

	save(ATE_HCMM_causal, SD_ATE_HCMM_causal, Check_ATE_HCMM_causal, file = paste0("2_BNPc_BD/Result_", i_rep,".RData"))
			 
} # for (i_rep)


###################################
# Causal inference result
###################################

rm(list=ls())  
load("1_DATA_BD.RData")

# True causal effect 
rep_seq = 1: 400
ATE_HCMM_rep = SD_ATE_HCMM_rep = IsCovered_ATE_HCMM_rep = numeric(length(rep_seq))

for (i_rep in rep_seq) {
  
  load(paste0("2_BNPc_BD/Result_", i_rep,".RData"))
  ATE_HCMM_rep[i_rep] = ATE_HCMM_causal
  SD_ATE_HCMM_rep[i_rep] = SD_ATE_HCMM_causal
  IsCovered_ATE_HCMM_rep[i_rep] = Check_ATE_HCMM_causal
  
}


###################################
# Monte Carlo simulation results
###################################
method = "BNPc"
effect = "ATE"
avg_ATE_HCMM= mean(ATE_HCMM_rep)
bias_ATE_HCMM = mean(ATE_HCMM_rep -ATE_true)
emp_sd_ATE_HCMM = sd(ATE_HCMM_rep)
ave_sd_ATE_HCMM = mean(SD_ATE_HCMM_rep)
coverage_ATE_HCMM= mean(IsCovered_ATE_HCMM_rep)
RMSE_ATE_HCMM = sqrt(mean((ATE_HCMM_rep - ATE_true)^2))
MAE_ATE_HCMM = median(abs(ATE_HCMM_rep - ATE_true))
###################################

HCMM_causal_res_MS = data.frame(Method = method, Effect = effect, Point_est = avg_ATE_HCMM, Emp_sd = emp_sd_ATE_HCMM, 
                             RMSE = RMSE_ATE_HCMM, MAE= MAE_ATE_HCMM, Coverage= coverage_ATE_HCMM, Avg_sd= ave_sd_ATE_HCMM)

write.csv(HCMM_causal_res_MS, file = "2_BNPc_BD_MS.csv")
