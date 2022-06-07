rm(list=ls())  

library( HCMMcausal )  # for running BNP causal 
library( norm ) # for mi.inference function

load("1_DATA.RData")

n_burnin = 2000 ; m_Imp = 10 ; interval_btw_Imp = 1000
rep_seq = 1:1000
n_seq = length(rep_seq)
n_sample = dim(CompleteDATA[[1]]$L1_7)[[1]]

ATE_BNPc_rep = SD_ATE_BNPc_rep = isCovered_ATE_BNPc_rep = numeric(n_seq)
dir.create("2_BNP_Est_BD")

for (i_rep in rep_seq){
	
  print(paste0("i_rep=",i_rep))
	set.seed(100 + i_rep)
	  
  obs_response = cbind(CompleteDATA[[i_rep]]$"Y_obs", CompleteDATA[[i_rep]]$"X_obs")
  TrueA = CompleteDATA[[i_rep]]$TrueA
  Incomplete_Y = CompleteDATA[[i_rep]]$L1_7[,6:7]
  Incomplete_X = CompleteDATA[[i_rep]]$L1_7[,1:5]

  data_obj = readData(Response_var=obs_response, trt_indicator=TrueA, Cont_pred=Incomplete_Y, Categ_pred=Incomplete_X, RandomSeed=100+i_rep)
  model_obj = createModel(data_obj)	
	result_obj = multipleImp(data_obj, model_obj, n_burnin, m_Imp, interval_btw_Imp, show_iter=FALSE)
	
  ATE_BNPc = apply(result_obj$est_delta,2,mean)
  SD_ATE_BNPc = apply(result_obj$est_delta,2,sd)
  Right_ATE_BNPc = ATE_BNPc + qnorm(0.975)*SD_ATE_BNPc
  Left_ATE_BNPc = ATE_BNPc - qnorm(0.975)*SD_ATE_BNPc
  Check_ATE_BNPc = rep(FALSE,2)
	for (i_resp in 1:2){
		Check_ATE_BNPc[i_resp] = (Left_ATE_BNPc[i_resp] < ATE_true[i_resp] & ATE_true[i_resp] < Right_ATE_BNPc[i_resp])
	}


  save(ATE_BNPc, SD_ATE_BNPc, Check_ATE_BNPc,  file = paste0("2_BNP_Est_BD/Result_", i_rep,".RData"))

} # for (i_rep)

##################################
# BNPc results with data before deletion 
##################################

rm(list=ls())
load("1_DATA.RData")

rep_seq = 1:1000
n_rep = length(rep_seq)
est_ATE = SE_ATE = isCov_ATE = array(0,c(n_rep,2))
est_predictor = SE_predictor = isCov_predictor = array(0,c(n_rep,2))

for (i_rep in rep_seq){
  
  load(paste0("2_BNP_Est_BD/Result_",i_rep,".RData"))
  est_ATE[i_rep,] = ATE_BNPc
  SE_ATE[i_rep,] = SD_ATE_BNPc
  isCov_ATE[i_rep,] = Check_ATE_BNPc
  
} # for

RESULT_ATE = array(0,c(5,2)) ; dimnames(RESULT_ATE)[[1]] = c("True","Est","SD","Coverage", "RMSE")
RESULT_ATE[1,] = c(1.5,0.5)
RESULT_ATE[2,] = apply(est_ATE,2,mean)
RESULT_ATE[3,] = apply(SE_ATE,2,mean)
RESULT_ATE[4,] = apply(isCov_ATE,2,mean)
# RMSE
RMSE_ATE1 = sqrt(mean((est_ATE[,1] - ATE_true[1])^2))
RMSE_ATE2 = sqrt(mean((est_ATE[,2] - ATE_true[2])^2))
RESULT_ATE[5,] = c(RMSE_ATE1, RMSE_ATE2)
round(RESULT_ATE,3)
write.csv(RESULT_ATE, file= "RESULT_ATE_BNP_BD.csv")
