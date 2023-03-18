rm(list=ls()) 

library(MASS) # for "mvrnorm"
library(mvtnorm) # for dmvnorm
# HCMM-LD package
library(MixedDataImpute) # HCMM-LD causal function
library( norm ) # for mi.inference function

load("1_DATA_BD.RData")
load("3_DATA_MAR.RData")

k_z = 2 ; k_y = 30 ; k_x = 50
num.impute = 20
num.burnin = 5000
num.skip = 50
thin.trace = 10
n_sample = 500
rep_seq = 1: 1000
n_rep = length(rep_seq)

mu_L3_rep = mu_L4_rep = se_mu_L3_rep = se_mu_L4_rep = IsCovered_mu_L3_rep = IsCovered_mu_L4_rep = rep(NA, n_rep)

for (i_rep in rep_seq) {
  
  # i_rep = 1
  set.seed(100 + i_rep )
  seed = 100 + i_rep
  print(paste0("i_rep=",i_rep))
  print(paste0("seed=",seed))
  
  L1_4 = incompleteDATA[[i_rep]]
  L4= L1_4[, 4]
  L3 = L1_4[,3]
  
  ###############################################
  ########## HCMM impute function ############
  ###############################################
  
  # Imputing potential outcome Y0 and Y1
  Y_cont = L1_4[, 3:4]
  L1 = L1_4[, 1]
  L2 = L1_4[, 2]
  L1 = factor(L1)
  L2 = factor(L2)
  X_cat = data.frame(L1 = L1, L2 = L2)
  
  HCMM_imp_res = hcmm_impute(X= X_cat, Y= Y_cont, kz= k_z, kx=k_x, ky= k_y, num.impute=num.impute, num.burnin=num.burnin, num.skip=num.skip, thin.trace=thin.trace)
  
  L3_mean_imp = L4_mean_imp = numeric(num.impute)
  L3_sd_imp = L4_sd_imp = numeric(num.impute)
  
  for (i in 1: num.impute){
    # Mean of L3 and L4
    L3_mean_imp[i] = mean(HCMM_imp_res$imputations[[i]]$L3)
    L3_sd_imp[i] = sqrt(var((HCMM_imp_res$imputations[[i]]$L3))/500)
    L4_mean_imp[i] = mean(HCMM_imp_res$imputations[[i]]$L4)
    L4_sd_imp[i] = sqrt(var((HCMM_imp_res$imputations[[i]]$L4))/500)
  }
  
  #################################################
  # Multiple imputation combining rule (Rubin 1987)
  #################################################
  
  true.mi = c(mean_L3, mean_L4)
  mi_mean_output = cbind(L3_mean_imp, L4_mean_imp)
  mi_se_output = cbind(L3_sd_imp, L4_sd_imp)
  
  mu.list <- vector("list", num.impute)
  ses.list <- vector("list", num.impute)
  
  for (imp_rep in 1:num.impute){
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
  
  #################################################
  # Repeated sampling results
  ################################################# 
  mu_L3_rep[i_rep] = q.mi[1]
  mu_L4_rep[i_rep] = q.mi[2]
  se_mu_L3_rep[i_rep] = se.mi[1]
  se_mu_L4_rep[i_rep] = se.mi[2]
  IsCovered_mu_L3_rep[i_rep] = check.mi[1]
  IsCovered_mu_L4_rep[i_rep] = check.mi[2]
  
} 

save(mu_L3_rep, mu_L4_rep, se_mu_L3_rep, se_mu_L4_rep, IsCovered_mu_L3_rep, IsCovered_mu_L4_rep, file = "HCCMimpute_r30s50k2.Rdata")

###################################
# Imputation performance
###################################
est_mu_L3 = mean(mu_L3_rep)
bias_mu_L3 = mean(mu_L3_rep - mean_L3)
emp_sd_mu_L3 = sd(mu_L3_rep)
avg_sd_mu_L3 = mean(se_mu_L3_rep)
coverage_mu_L3 = mean(IsCovered_mu_L3_rep)
RMSE_mu_L3 = sqrt(mean((mu_L3_rep - mean_L3)^2))

est_mu_L4 = mean(mu_L4_rep)
bias_mu_L4 = mean(mu_L4_rep - mean_L4)
emp_sd_mu_L4 = sd(mu_L4_rep)
avg_sd_mu_L4 = mean(se_mu_L4_rep)
coverage_mu_L4 = mean(IsCovered_mu_L4_rep)
RMSE_mu_L4 = sqrt(mean((mu_L4_rep - mean_L4)^2))

effect_imp = c("mu_L3", "mu_L4")
est_imp = c(est_mu_L3, est_mu_L4)
bias_imp = c(bias_mu_L3, bias_mu_L4)
emp_sd_imp = c(emp_sd_mu_L3, emp_sd_mu_L4)
avg_sd_imp = c(avg_sd_mu_L3, avg_sd_mu_L4)
coverage_imp = c(coverage_mu_L3, coverage_mu_L4)
RMSE_imp = c(RMSE_mu_L3, RMSE_mu_L4)

imp_res = data.frame(Method = "HCMM_impute", Est= effect_imp, point= est_imp, Emp_sd = emp_sd_imp, RMSE = RMSE_imp, Coverage = coverage_imp, Avg_sd = avg_sd_imp)
write.csv(imp_res, file = "HCMM_imputation_r30s50k2.csv")

save.image(file = "HCMM_impute_all_r30s50k2.RData")

#######################################################3
i_rep = length(rep_seq)
L1_4_orginal = DATA[[i_rep]]
L1_orginal = L1_4_orginal[, 1]
L3_orginal = L1_4_orginal[, 3]
L4_orginal = L1_4_orginal[, 4]
L1_4_missing = incompleteDATA[[i_rep]]
miss_indicator = which(rowSums(is.na(L1_4_missing[, c(3, 4)]))>0)
seq = seq(from= 1, to= num.impute, by= 1)

dir.create("check_r30s50k2") 

for (which_imp in seq ){
  
  # which_imp = 1  
  L3_imp = HCMM_imp_res$imputations[[which_imp]]$L3
  L4_imp = HCMM_imp_res$imputations[[which_imp]]$L4
  
  png(file=paste0("check_r30s50k2/Bivariate_plot",which_imp,".png"), width=1200,height=800,pointsize=20)
  plot(L3_orginal, L4_orginal, xlab = "L3", ylab = "L4", col= "grey", pch=16, 
       main = "",  xlim = c(-9, 9), ylim = c(-9, 9), cex= 1, cex.lab=1)
  points(L3_orginal[miss_indicator], L4_orginal[miss_indicator], col= "blue", pch=1, cex= 1)
  points(L3_imp[miss_indicator], L4_imp[miss_indicator], col= "red", pch=16, cex= 1)
  dev.off()
  
  # which_imp = 1
  par(mfrow=c(1,1))
  png(file=paste0("check_r30s50k2/Barplot",which_imp,".png"), width=1200,height=800,pointsize=20)
  L1_true_dat = data.frame(L1 = L1_orginal, group= rep("True", n= length(L1)))
  L1_imp_dat = data.frame(L1 = HCMM_imp_res$imputations[[which_imp]]$L1, group= rep("Completed", n= length(L1)))
  L1_combined_dat = rbind(L1_true_dat, L1_imp_dat)
  counts = with(L1_combined_dat, table(group, L1))
  barplot(counts, main="L1", xlab="L1", col=c("blue","grey"), legend = rownames(counts), beside=TRUE,
          args.legend = list(x = "topright", bty = "n", inset=c(0, 0)), ylim = c(0, 400))
  title(sub= paste("Imputed data", which_imp))
  dev.off()
  
}



n_trace= num.impute*num.skip/thin.trace
X_alloc_trace = Y_alloc_trace = Z_alloc_trace = numeric(n_trace)

for (i in 1: n_trace) {
  X_alloc_trace[i]= length(which(HCMM_imp_res$trace$X_alloc[i,]>0))
  Y_alloc_trace[i]= length(which(HCMM_imp_res$trace$Y_alloc[i,]>0))
  Z_alloc_trace[i]= length(which(HCMM_imp_res$trace$Z_alloc[i,]>0))
}

png(file=paste0("check_r30s50k2/traceplot_k2_s_r.png"), width=1200,height=800,pointsize=20)
par(mfrow=c(1,1))
plot(Z_alloc_trace, type = "l", ylim = c(0, 25), main = "n_comp: k-black, r-red, s-blue ", ylab = "# of occupied comp", xlab = "Iteration")
lines(Y_alloc_trace, col="red")
lines(X_alloc_trace, col="blue")
# abline( h=c(k_z,k_y,k_x), lty="dashed", col=c("black","red","blue"), lwd=2 )
dev.off()
