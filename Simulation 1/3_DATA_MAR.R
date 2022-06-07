rm(list=ls())  

load("1_DATA_BD.RData")

ls(DATA[[1]])

n_sample = dim(DATA[[1]]$L1_6)[[1]]
incompleteDATA = NULL

complete.cases_rate = numeric(length(rep_seq))
  
expit_fn = function(x) { 1/(1+exp(-x)) } # inverse of the logit function

for (i_rep in rep_seq) {
  
  set.seed(100 + i_rep )
  	
	L1_6 = DATA[[i_rep]]$L1_6 
	L1_4 = L1_6[,1:4]
	
	L1_4_missing_ind = array(0,c(n_sample,4))
	L1_4_missing_ind[,1] = rbinom(n_sample, 1, expit_fn(-1 + 3*L1_4[,2]))
  L1_4_missing_ind[,2] = rbinom(n_sample, 1, expit_fn(-1.5 + 3*L1_4[,3]))
  L1_4_missing_ind[,4] = rbinom(n_sample, 1, expit_fn(-1.2 - 2*L1_4[,1] - L1_4[,3]))
  L5_missing_ind = rbinom(n_sample, 1, expit_fn(-1.5 + 3*L1_4[,1] - 1.5*L1_4[,3]))
  L6_missing_ind = rbinom(n_sample, 1, expit_fn(-1.2 + 2.2*L1_4[,1] - 0.5*L1_4[,3]))
  L1_6_missing_ind = cbind(L1_4_missing_ind, L5_missing_ind, L6_missing_ind)
  L1_6[L1_6_missing_ind==1] = NA
	
	Y_obs = DATA[[i_rep]]$Y_obs ; TrueA = DATA[[i_rep]]$TrueA
	incompleteDATA[[i_rep]] = list(Y_obs=Y_obs, TrueA=TrueA, L1_6=L1_6)
  
	complete.cases_rate[i_rep] = sum(complete.cases(L1_6))/n_sample
}
# for 

mean(complete.cases_rate)
# [1] 0.091066

save(incompleteDATA, rep_seq, ATE_true, file = "3_DATA_MAR.RData")
