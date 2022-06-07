rm(list=ls())  

library( HCMMcausal )  # for running BNP causal 
library( norm ) # for mi.inference function

load("1_DATA.RData")

n_burnin = 2000 ; m_Imp = 10 ; interval_btw_Imp = 1000
rep_seq = 1: 1000
n_seq = length(rep_seq)
n_sample = dim(IncompleteDATA[[1]]$L1_7)[[1]]

ATE_BNPc_rep = SD_ATE_BNPc_rep = isCovered_ATE_BNPc_rep = numeric(n_seq)
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
	
  ATE_BNPc = apply(result_obj$est_delta,2,mean)
  SD_ATE_BNPc = apply(result_obj$est_delta,2,sd)
  Right_ATE_BNPc = ATE_BNPc + qnorm(0.975)*SD_ATE_BNPc
  Left_ATE_BNPc = ATE_BNPc - qnorm(0.975)*SD_ATE_BNPc
  Check_ATE_BNPc = rep(FALSE,2)
	for (i_resp in 1:2){
		Check_ATE_BNPc[i_resp] = (Left_ATE_BNPc[i_resp] < ATE_true[i_resp] & ATE_true[i_resp] < Right_ATE_BNPc[i_resp])
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
  
  save(ATE_BNPc, SD_ATE_BNPc, Check_ATE_BNPc, true.mi, q.mi, se.mi, check.mi, file = paste0("2_BNP_Est_MAR/Result_", i_rep,".RData"))

} # for (i_rep)

##################################
# BNPc results with MAR data
##################################

rm(list=ls())

load("1_DATA.RData")

rep_seq = 1:1000
n_rep = length(rep_seq)

ATE_Y_BNPc_rep = sd_ATE_Y_BNPc_rep = IsCovered_ATE_Y_BNPc_rep = numeric(length(rep_seq))
ATE_X_BNPc_rep = sd_ATE_X_BNPc_rep = IsCovered_ATE_X_BNPc_rep = numeric(length(rep_seq))

for (i_rep in rep_seq){
  
  # i_rep = 1
  load(paste0("2_BNP_Est_MAR/Result_",i_rep,".RData"))
  
  # ATE result
  ATE_Y_BNPc_rep[i_rep] = ATE_BNPc[1]
  ATE_X_BNPc_rep[i_rep] = ATE_BNPc[2]
  sd_ATE_Y_BNPc_rep[i_rep] = SD_ATE_BNPc[1]
  sd_ATE_X_BNPc_rep[i_rep] = SD_ATE_BNPc[2]
  IsCovered_ATE_Y_BNPc_rep[i_rep] = Check_ATE_BNPc[1]
  IsCovered_ATE_X_BNPc_rep[i_rep] = Check_ATE_BNPc[2]
  
} # for 

ATE_Y_true = 1.5
ATE_X_true = 0.5

# Causal inference 
ATE_Y_BNPc = mean(ATE_Y_BNPc_rep)
bias_ATE_Y_BNPc = mean(ATE_Y_BNPc_rep -ATE_Y_true)
emp_sd_ATE_Y_BNPc = sd(ATE_Y_BNPc_rep)
avg_sd_ATE_Y_BNPc = mean(sd_ATE_Y_BNPc_rep)
coverage_ATE_Y_BNPc= mean(IsCovered_ATE_Y_BNPc_rep)
RMSE_ATE_Y_BNPc = sqrt(mean((ATE_Y_BNPc_rep - ATE_Y_true)^2))
MAE_ATE_Y_BNPc = median(abs(ATE_Y_BNPc_rep - ATE_Y_true))

ATE_X_BNPc = mean(ATE_X_BNPc_rep)
bias_ATE_X_BNPc = mean(ATE_X_BNPc_rep -ATE_X_true)
emp_sd_ATE_X_BNPc = sd(ATE_X_BNPc_rep)
avg_sd_ATE_X_BNPc = mean(sd_ATE_X_BNPc_rep)
coverage_ATE_X_BNPc= mean(IsCovered_ATE_X_BNPc_rep)
RMSE_ATE_X_BNPc = sqrt(mean((ATE_X_BNPc_rep - ATE_X_true)^2))
MAE_ATE_X_BNPc = median(abs(ATE_X_BNPc_rep - ATE_X_true))

Est_BNPc = c(ATE_Y_BNPc, ATE_X_BNPc)
bias_BNPc = c(bias_ATE_Y_BNPc, bias_ATE_X_BNPc)
emp_sd_BNPc = c(emp_sd_ATE_Y_BNPc, emp_sd_ATE_X_BNPc)
avg_sd_BNPc = c(avg_sd_ATE_Y_BNPc, avg_sd_ATE_X_BNPc)
coverage_BNPc = c(coverage_ATE_Y_BNPc, coverage_ATE_X_BNPc)
RMSE_BNPc = c(RMSE_ATE_Y_BNPc, RMSE_ATE_X_BNPc)
MAE_BNPc = c(MAE_ATE_Y_BNPc, MAE_ATE_X_BNPc)

causal_res = data.frame(method = c("BNP", "BNP"), effect= c("ATE_Y", "ATE_X"), Est= Est_BNPc, 
                        emp_sd= emp_sd_BNPc, RMSE = RMSE_BNPc, MAE = MAE_BNPc , coverage = coverage_BNPc, 
                        avg_sd = avg_sd_BNPc)

write.csv(causal_res, file = "3_BNP_causal_res.csv")

###########################################################

rm(list=ls())
load("1_DATA.RData")

rep_seq = 1:1000
n_rep = length(rep_seq)

est_ATE = SE_ATE = isCov_ATE = array(0,c(n_rep,2))
est_predictor = SE_predictor = isCov_predictor = array(0,c(n_rep,3))

for (i_rep in rep_seq){
  
  load(paste0("2_BNP_Est_MAR/Result_",i_rep,".RData"))
  
  # ATE result
  est_ATE[i_rep,] = ATE_BNPc
  SE_ATE[i_rep,] = SD_ATE_BNPc
  isCov_ATE[i_rep,] = Check_ATE_BNPc
  # Imputation result
  est_predictor[i_rep,] = q.mi
  SE_predictor[i_rep,] = se.mi
  isCov_predictor[i_rep,] = check.mi
  
} # for

RESULT_ATE = array(0,c(5,2)) ; dimnames(RESULT_ATE)[[1]] = c("True","Est","SD","Coverage", "RMSE")
RESULT_ATE[1,] = c(1.5,1.0)
RESULT_ATE[2,] = apply(est_ATE,2,mean)
RESULT_ATE[3,] = apply(SE_ATE,2,mean)
RESULT_ATE[4,] = apply(isCov_ATE,2,mean)
round(RESULT_ATE,3)
# RMSE
RMSE_ATE1 = sqrt(mean((est_ATE[,1] - ATE_true[1])^2))
RMSE_ATE2 = sqrt(mean((est_ATE[,2] - ATE_true[2])^2))
RESULT_ATE[5,] = c(RMSE_ATE1, RMSE_ATE2)
round(RESULT_ATE,3)
write.csv(RESULT_ATE, file= "RESULT_ATE_BNP_MAR.csv")

RESULT_predictor = array(0,c(4,3))
dimnames(RESULT_predictor)[[1]] = c("True","Est","SD","Coverage")
dimnames(RESULT_predictor)[[2]] = names(se.mi)
RESULT_predictor[1,] = true.mi
RESULT_predictor[2,] = apply(est_predictor,2,mean)
RESULT_predictor[3,] = apply(SE_predictor,2,mean)
RESULT_predictor[4,] = apply(isCov_predictor,2,mean)
round(RESULT_predictor,3)
write.csv(RESULT_predictor, file= "3_BNP_imp_res.csv")