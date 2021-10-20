rm(list=ls())  

library( HCMMcausal )  # for running HCMM-LD causal 
library( norm ) # for mi.inference function

load("1_DATA.RData")

n_burnin = 2000 ; m_Imp = 10 ; interval_btw_Imp = 200
rep_seq = 1: 400
n_seq = length(rep_seq)
n_sample = dim(CompleteDATA[[1]]$L1_7)[[1]]

ATE_HCMM_causal_rep = SD_ATE_HCMM_causal_rep = isCovered_ATE_HCMM_causal_rep = numeric(n_seq)
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
	
  ATE_HCMM_causal = apply(result_obj$est_delta,2,mean)
  SD_ATE_HCMM_causal = apply(result_obj$est_delta,2,sd)
  Right_ATE_HCMM_causal = ATE_HCMM_causal + qnorm(0.975)*SD_ATE_HCMM_causal
  Left_ATE_HCMM_causal = ATE_HCMM_causal - qnorm(0.975)*SD_ATE_HCMM_causal
  Check_ATE_HCMM_causal = rep(FALSE,2)
	for (i_resp in 1:2){
		Check_ATE_HCMM_causal[i_resp] = (Left_ATE_HCMM_causal[i_resp] < ATE_true[i_resp] & ATE_true[i_resp] < Right_ATE_HCMM_causal[i_resp])
	}


  save(ATE_HCMM_causal, SD_ATE_HCMM_causal, Check_ATE_HCMM_causal,  file = paste0("2_BNP_Est_BD/Result_", i_rep,".RData"))

} # for (i_rep)

