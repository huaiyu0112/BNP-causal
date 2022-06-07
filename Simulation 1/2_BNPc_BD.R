rm(list=ls())  

library( HCMMcausal )  # for running BNP causal 

load("1_DATA_BD.RData")
dir.create("2_BNPc_BD")

# How many times the simulation are repeated 
rep_seq = 1: 1000
n_seq = length(rep_seq)
n_burnin = 2000 ; m_Imp = 10 ; interval_btw_Imp = 1000

for (i_rep in rep_seq) {

  print(paste0("i_rep=",i_rep))
	set.seed(100 + i_rep)
	  
  obs_response = DATA[[i_rep]]$Y_obs 
  TrueA = DATA[[i_rep]]$TrueA
  Incomplete_Y = DATA[[i_rep]]$L1_6[,1:4]
  Incomplete_X = DATA[[i_rep]]$L1_6[,5:6]

  data_obj = readData(Response_var=obs_response, trt_indicator=TrueA, Cont_pred=Incomplete_Y, Categ_pred=Incomplete_X, RandomSeed=100+i_rep)
  model_obj = createModel(data_obj)
	result_obj = multipleImp(data_obj, model_obj, n_burnin, m_Imp, interval_btw_Imp, show_iter=FALSE)
  	
	ATE_BNPc = mean(result_obj$est_delta)
	SD_ATE_BNPc = sd(result_obj$est_delta)
	Right_ATE_BNPc = ATE_BNPc + qnorm(0.975)*SD_ATE_BNPc
	Left_ATE_BNPc = ATE_BNPc - qnorm(0.975)*SD_ATE_BNPc
	Check_ATE_BNPc = (Left_ATE_BNPc < ATE_true & ATE_true < Right_ATE_BNPc)

	save(ATE_BNPc, SD_ATE_BNPc, Check_ATE_BNPc, file = paste0("2_BNPc_BD/Result_", i_rep,".RData"))
	
} # for (i_rep)

###################################
# Causal inference result
###################################

rm(list=ls())  
load("1_DATA_BD.RData")

rep_seq = 1: 1000
ATE_BNPc_rep = SD_ATE_BNPc_rep = IsCovered_ATE_BNPc_rep = numeric(length(rep_seq))

for (i_rep in rep_seq) {
  
  load(paste0("2_BNPc_BD/Result_", i_rep,".RData"))
  ATE_BNPc_rep[i_rep] = ATE_BNPc
  SD_ATE_BNPc_rep[i_rep] = SD_ATE_BNPc
  IsCovered_ATE_BNPc_rep[i_rep] = Check_ATE_BNPc
  
}

###################################
# Monte Carlo simulation results
###################################
method = "BNPc"
effect = "ATE"
avg_ATE_BNPc= mean(ATE_BNPc_rep)
bias_ATE_BNPc = mean(ATE_BNPc_rep -ATE_true)
emp_sd_ATE_BNPc = sd(ATE_BNPc_rep)
ave_sd_ATE_BNPc = mean(SD_ATE_BNPc_rep)
coverage_ATE_BNPc= mean(IsCovered_ATE_BNPc_rep)
RMSE_ATE_BNPc = sqrt(mean((ATE_BNPc_rep - ATE_true)^2))
MAE_ATE_BNPc = median(abs(ATE_BNPc_rep - ATE_true))
###################################

BNPc_res = data.frame(Method = method, Effect = effect, Point_est = avg_ATE_BNPc, Emp_sd = emp_sd_ATE_BNPc, 
                                RMSE = RMSE_ATE_BNPc, MAE= MAE_ATE_BNPc, Coverage= coverage_ATE_BNPc, Avg_sd= ave_sd_ATE_BNPc)

write.csv(BNPc_res, file = "2_BNPc_BD_res.csv")
