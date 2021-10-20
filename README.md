# BNP-causal
R Codes used to simulate and analyze the data for the manuscript "A Bayesian Causal Inference Approach in Observational Studies with Missingness in Covariates and Outcomes".


To run our proposed BNP causal model, (1) we install "HCMMcausal" package using "HCMMcausal_1.5.2.tar.gz" file, (2) we then make the data object using readData function to load the data structure, in which we input response variable in "Response_var" argument, treatment indicator variable in "trt_indicator" argument, continuous covariates in "Cont_pred" argument and categorical covariates in "Categ_pred" argument. (3) we make a model object using createModel function to load the data object. In createModel function, we could also specify the hyperparameters and upper bounds of mixture components. (4) we run multipleImp function for proposed BNP causal model where we load the data and model object and then save the results. In multipleImp function, we could input the number of burn-in iterations in "n_burnin" argument, number of multiple imputations in "m_Imp" argument, and interval (number of iterations) between imputed data in "interval_btw_Imp" argument.

