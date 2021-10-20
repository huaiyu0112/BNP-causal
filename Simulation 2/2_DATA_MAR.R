rm(list=ls())  

load("1_DATA_BD.RData")

ls(DATA[[1]])

n_sample = dim(DATA[[1]]$L1_7)[[1]]
incompleteDATA = NULL

expit_fn = function(x) { 1/(1+exp(-x)) } # inverse of the logit function

for (i_rep in rep_seq) {
  
  set.seed(100 + i_rep )
  	
	L1_7 = DATA[[i_rep]]$L1_7 

	###############################################
	########## Missingness ##############
	###############################################
	L1_7_missing_ind = array(0,c(n_sample, dim(L1_7)[2]))
	
	L1_7_missing_ind[,1] = rbinom(n_sample, 1, expit_fn(-2.5 + 2.9*ifelse(L1_7[,2]== 2, 1, 0) + 2*ifelse(L1_7[,2]== 3, 1, 0) + 1.5*ifelse(L1_7[,2]== 4, 1, 0)))
	L1_7_missing_ind[,2] = rbinom(n_sample, 1, expit_fn(-2.8 + 1.7* ifelse(L1_7[,3]== 2, 1, 0) + 1.9*ifelse(L1_7[,3]== 3, 1, 0)))
	L1_7_missing_ind[,3] = rbinom(n_sample, 1, expit_fn(-1.8 + 1.2*ifelse(L1_7[,4]== 2, 1, 0) ))
	L1_7_missing_ind[,4] = rbinom(n_sample, 1, expit_fn(-2 + 1*ifelse(L1_7[,5]== 2, 1, 0) + 1.5*ifelse(L1_7[,5]== 3, 1, 0) ))
	L1_7_missing_ind[,5] = rbinom(n_sample, 1, expit_fn(-3.7 + 0.7*L1_7[,6] ))
	L1_7_missing_ind[,6] = rbinom(n_sample, 1, expit_fn(-2.5 + 0.5*L1_7[,7] ))
	L1_7[L1_7_missing_ind==1] = NA
	
	Y_obs = DATA[[i_rep]]$Y_obs ; TrueA = DATA[[i_rep]]$TrueA
	
	incompleteDATA[[i_rep]] = list(Y_obs=Y_obs, TrueA=TrueA, L1_7=L1_7)
  
} # for 

save(incompleteDATA, rep_seq, ATE_true, file = "3_DATA_MAR.RData")
