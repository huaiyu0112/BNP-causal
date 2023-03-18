rm(list=ls())  

load("1_DATA_BD.RData")

n_sample = dim(DATA[[1]])[[1]]
incompleteDATA = NULL

rep_seq = 1: 1000
complete.cases_rate = numeric(length(rep_seq))

expit_fn = function(x) { 1/(1+exp(-x)) } # inverse of the logit function

for (i_rep in rep_seq) {
  
  # i_rep = 1
  set.seed(100 + i_rep )
  	
	L1_4 = DATA[[i_rep]]

	###############################################
	########## Missingness ##############
	###############################################
	L1_4_missing_ind = array(0,c(n_sample, dim(L1_4)[2]))
	
	L1_4_missing_ind[,1] = rbinom(n_sample, 1, expit_fn(-1.3 + 0.5 * L1_4[,4] + 2.5*ifelse(L1_4[,2]== 2, 1, 0)))
	L1_4_missing_ind[,3] = rbinom(n_sample, 1, expit_fn(-1.2 - 0.8 * L1_4[,4]- 1.5*ifelse(L1_4[,2]== 3, 1, 0)))
	L1_4_missing_ind[,4] = rbinom(n_sample, 1, 0.15)
	
	L1_4[L1_4_missing_ind==1] = NA
	incompleteDATA[[i_rep]] = L1_4
	
	complete.cases_rate[i_rep] = sum(complete.cases(L1_4))/500
  
} # for 

mean(complete.cases_rate)


head(L1_4, n =50)
# [1] 0.43938 

save(incompleteDATA, file = "3_DATA_MAR.RData")
