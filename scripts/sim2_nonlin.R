dgp_nonlin <- function(n, sigma = 1){
  W1 <- rnorm(n)
  W2 <- rnorm(n, -2, sd = 2)
  p  <- plogis(1+0.5*W1 - 0.75*W2 + 0.3*W1*W2)      #non-linearity in the interaction term
  A  <- rbinom(n, 1, p)
  
  mu0<- sin(W1) + 0.5*W2^2                          #control outcome
  mu1<- sin(W1) + 0.5*W2^2 + 0.3*W1^2               #treatment outcome
  
  Y  <- rnorm(n, mean = ifelse(A ==1, mu1, mu0), sd = sigma)
  data.table(W1, W2, A, Y)
}


ATE_true <- 0.3 # E(mu1 - mu0) = 0.3 E(w1^2) = 0.3

######################################################################### Helper

ATE_helper <- function(Y, A, g_hat, Q0_hat, Q1_hat){##################### 1-step
  ## one-step estimator
  psi_vec <- (A / g_hat) * (Y - Q1_hat) -
    ((1 - A) / (1 - g_hat)) * (Y - Q0_hat) +
    (Q1_hat - Q0_hat)
  psi <- mean(psi_vec)
  se  <- sd(psi_vec)/sqrt(length(psi_vec))
  ci  <- psi + qnorm(0.975) * c(-1,1) * se
  
  list(ate = psi, se = se, ci = ci)
}

clip01 <- function(x, eps = 0.0001) pmin(pmax(x, eps), 1 - eps)######## ensures positivity


#################################################### Cross-fit 1-step w/ XGBoost


ATE_crossfit_xgb <- function(df, K = 5L, nrounds_g = 250, nrounds_Q = 250){
  
  n <- nrow(df)
  
  folds <- make_folds(df, fold_fun =  folds_vfold, V = K) #### using origami
  
  for (fold in seq_len(K)){
    idx_tr <- folds[[fold]]$training_set
    idx_te <- folds[[fold]]$validation_set
    
    #####################################################propensity score estimation
    
    dtrain_g <- xgb.Dmatrix(as.matrix(df[idx_tr, .(W1, W2)]), label = df$A[idx_tr])
    
    mod_g <- xgb.train(params = list(objective = "binary:logistic", 
                                     eval_metric = "logloss", 
                                     eta = 0.05, max_depth = 3),
                       data = dtrain_g,
                       nrounds = nrounds_g,
                       verbose = 0
    )
    
    g_hat[idx_te] <- predict(mod_g, newdata = as.matrix(df[idx_te, .(W1, W2)]))
    
    g_hat[idx_te] <- clips01(g_hat[idx_te]) #################ensure positivity
    
    ##################################################outcome regression estimation
    
    dtrain_Q <- xgb.Dmatrix(as.matrix(df[idx_tr, .(A, W1, W2)]), label = df$Y[idx_tr])
    
    mod_Q <- xgb.train(params = list(objective = "reg:squarederror",
                                     eval_metric = "rmse",
                                     eta = 0.05, max_depth = 4),
                       data = dtrain_Q,
                       nrounds = nrounds_Q,
                       verbose = 0
    )
    
    ##predict for both treatment levels
    
    df_te1 <- copy(df[idx_te]); df_te1[, A:= 1]
    df_te0 <- copy(df[idx_te]); df_te0[, A:= 0]
    
    Q1_hat[idx_te] <- predict(mod_Q, newdata = as.matrix(dt_te1[, .(A, W1, W2)]))
    Q0_hat[idx_te] <- predict(mod_Q, newdata = as.matrix[dt_te0[, .(A, W1, W2)]])
    
  }  
  
  ATE_helper(df$Y, df$A, g_hat, Q0_hat, Q1_hat) 
}





###############################################################TMLE with XGboost 

lrnr_xgb_bin <- Lrnr_xgboost$new(objective = "binary:logistic",
                                 eval_metric = "logloss",
                                 nrounds = 250, eta = 0.05, max_depth = 3)

lrnr_xgb_reg <- Lrnr_xgboost$new(objective = "reg:squarederror", 
                                 metric = "rmse",
                                 nrounds = 250, eta = 0.05, max_depth = 4)



learners <- list(Y = lrnr_xgb_reg, A = lrnr_xgb_bin) 

ATE_spec <- tmle_ATE(treatment_level = 1, control_level = 0)  




#####################################################################example run


n <- 4000
sim <- dgp_nonlin(n, sigma = 1)

##one-step
results_os <- ATE_crossfit_xgb(sim, K = 5)
print(results_os)

##TMLE
nodes <- list(W = c("W1", "W2"), A = "A", Y = "Y")
tmle_fit <- tmle3(ATE_spec, sim, nodes, learners)

print(tmle_fit$summary[, c("type", "psi_transformed", "se")])


cat(sprintf("\nTrue ATE = %.2f | One-step = %.3f | TMLE = %.3f\n", ATE_true, res_os$ate, tmle_fit$summary$psi_transformed))




##not sure why but one-step is really bad here
