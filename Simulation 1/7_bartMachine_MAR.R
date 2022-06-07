rm(list=ls()) 

# https://github.com/kapelner/bartMachine

# install.packages("rJava")
# install.packages("bartMachine")
# install.packages("~/Downloads/rJava_0.9-12.tar.gz", repos = NULL, type = "source")

Sys.setenv(JAVA_HOME='C:\\Program Files (x86)\\Java\\jdk-12') # for 64-bit version


library(rJava); packageVersion("rJava")
library(randomForest)
library(bartMachine); packageVersion("bartMachine")
library(tidytreatment)
library(dplyr)
options(java.parameters = "-Xmx2500m")

load("1_DATA_BD.RData")
load("3_DATA_MAR.RData")

# True causal effect 
ATE_true = 1.5

expit_fn = function(x) { 1/(1+exp(-x)) } # inverse of the logit function

# How many times the simulation are repeated 
start_seq =0
rep_seq = c(1:1000) #1:40 # 1:40
n_seq = length(rep_seq)

ATE_bartMachine_rep = sd_ATE_bartMachine_rep = IsCovered_ATE_bartMachine_rep = numeric(length(rep_seq))



dir.create("7_bartMachine_MAR")

# Define number of observations for each dataset 
n_sample = 500

for (i_rep in rep_seq) {
  
  # i_rep = 1
  set.seed(100 + i_rep + start_seq)
  seed = 100 + i_rep + start_seq
  print(paste0("i_rep=",i_rep))
  print(paste0("seed=",seed))
  
  Y_obs = incompleteDATA[[i_rep]]$Y_obs 
  A = incompleteDATA[[i_rep]]$TrueA
  Incomplete_Y = incompleteDATA[[i_rep]]$L1_6[,1:4]
  Incomplete_X = incompleteDATA[[i_rep]]$L1_6[,5:6]
  L1_6_missing = data.frame(cbind(Incomplete_Y, Incomplete_X))
  L1_6 = DATA[[i_rep]]$L1_6[,1:6]
  
  dat = data.frame(Y_obs, A, L1_6_missing)
  
  bartM <- bartMachine(
    X = select(dat,-Y_obs),
    y = select(dat, Y_obs)[[1]], use_missing_data = TRUE) 
  
  posterior_treat_eff <- treatment_effects(bartM, treatment = "A")
  ATE_bartMachine_res = aggregate(posterior_treat_eff$cte, list(posterior_treat_eff$.draw), FUN=mean)
  ATE_bartMachine_draw = ATE_bartMachine_res$x
  
  ATE_bartMachine = mean(ATE_bartMachine_draw)
  SD_ATE_bartMachine = sd(ATE_bartMachine_draw)
  Right_ATE_bartMachine = ATE_bartMachine + qnorm(0.975)*SD_ATE_bartMachine
  Left_ATE_bartMachine = ATE_bartMachine - qnorm(0.975)*SD_ATE_bartMachine
  Check_ATE_bartMachine = (Left_ATE_bartMachine < ATE_true & ATE_true < Right_ATE_bartMachine)
  
  save(ATE_bartMachine, SD_ATE_bartMachine, Check_ATE_bartMachine,  file = paste0("7_bartMachine_MAR/Result_", i_rep,".RData"))
  
}

#######################################################
# Result Summary
#######################################################
rm(list=ls())  

load("1_DATA_BD.RData")

ATE_true = 1.5 # true causal effect 
rep_seq = 1: 1000
ATE_bartMachine_rep = SD_ATE_bartMachine_rep = IsCovered_ATE_bartMachine_rep = numeric(length(rep_seq))

for (i_rep in rep_seq) {
  
  load(paste0("7_bartMachine_MAR/Result_", i_rep,".RData"))
  
  ATE_bartMachine_rep[i_rep] = ATE_bartMachine
  SD_ATE_bartMachine_rep[i_rep] = SD_ATE_bartMachine
  IsCovered_ATE_bartMachine_rep[i_rep] = Check_ATE_bartMachine
  
}

###################################
# Causal inference
###################################
method = "bartMachine"
effect = "ATE"
avg_ATE_bartMachine= mean(ATE_bartMachine_rep)
bias_ATE_bartMachine = mean(ATE_bartMachine_rep -ATE_true)
emp_sd_ATE_bartMachine = sd(ATE_bartMachine_rep)
avg_sd_ATE_bartMachine = mean(SD_ATE_bartMachine_rep)
coverage_ATE_bartMachine= mean(IsCovered_ATE_bartMachine_rep)
RMSE_ATE_bartMachine = sqrt(mean((ATE_bartMachine_rep - ATE_true)^2))
MAE_ATE_bartMachine = median(abs(ATE_bartMachine_rep - ATE_true))

bartMachine_res = data.frame(Method= method, effect, Point_est = avg_ATE_bartMachine, Emp_sd = emp_sd_ATE_bartMachine,
                             RMSE = RMSE_ATE_bartMachine, MAE= MAE_ATE_bartMachine, Coverage= coverage_ATE_bartMachine, Avg_sd= avg_sd_ATE_bartMachine)
write.csv(bartMachine_res, file = "7_bartMachine_MAR_res.csv")
