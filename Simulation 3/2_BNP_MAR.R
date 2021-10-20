rm(list=ls())  

library( HCMMcausal )  # for running HCMM-LD causal 
library( norm ) # for mi.inference function

load("1_DATA.RData")

n_burnin = 2000 ; m_Imp = 10 ; interval_btw_Imp = 200
rep_seq = 1:400
n_seq = length(rep_seq)
n_sample = dim(IncompleteDATA[[1]]$L1_7)[[1]]

ATE_HCMM_causal_rep = SD_ATE_HCMM_causal_rep = isCovered_ATE_HCMM_causal_rep = numeric(n_seq)
dir.create("2_BNP_Est_MAR")

for (i_rep in rep_seq){
	
  print(paste0("i_rep=",i_rep))
	set.seed(100 + i_rep)
	  
  obs_response = IncompleteDATA[[i_rep]]$Outcome 
  TrueA = IncompleteDATA[[i_rep]]$TrueA
  Incomplete_Y = IncompleteDATA[[i_rep]]$L1_7[,6:7]
  Incomplete_X = IncompleteDATA[[i_rep]]$L1_7[,1:5]

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

  #################################################
  # Imputation performance summary
  #################################################
  
  # Combining rule
  mu_Y_out  = mu_X_out = mu_L6_out = numeric(m_Imp)
  se_Y_out  = se_X_out = se_L6_out = numeric(m_Imp)

  for (imp_rep in 1:m_Imp){

    mu_Y_out[imp_rep] = mean(result_obj$MI_Y_cube[imp_rep,, 1])
    se_Y_out[imp_rep] = sqrt(var(result_obj$MI_Y_cube[imp_rep,, 1])/n_sample)
    
    mu_X_out[imp_rep] = mean(result_obj$MI_Y_cube[imp_rep,, 2])
    se_X_out[imp_rep] = sqrt(var(result_obj$MI_Y_cube[imp_rep,, 2])/n_sample)
    
    mu_L6_out[imp_rep] = mean(result_obj$MI_Y_cube[imp_rep,, (2+1)])
    se_L6_out[imp_rep] = sqrt(var(result_obj$MI_Y_cube[imp_rep,, (2+1)])/n_sample)

  }

  true.mi = c(mean_Y_true, mean_X_true, mean_L6_true)
  mi_mean_output = cbind(mu_Y_out, mu_X_out, mu_L6_out)
  mi_se_output = cbind(se_Y_out, se_X_out, se_L6_out)

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
    check.mi[i_mi] = (lower.mi[i_mi] < true.mi[i_mi] & true.mi[i_mi] < upper.mi[i_mi])
  }
  
  save(ATE_HCMM_causal, SD_ATE_HCMM_causal, Check_ATE_HCMM_causal, true.mi, q.mi, se.mi, check.mi, file = paste0("2_BNP_Est_MAR/Result_", i_rep,".RData"))

} # for (i_rep)

