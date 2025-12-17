library(glmnet)
library(data.table)
library(ggplot2)
library(gridExtra)
library(here)
library(patchwork)


########################################################################### DGP

# Utility functions to generate A and Y
expit <- function(x) 1/(1+exp(-x))
yerr <- function(n){4*(rbeta(n,2,2) - 0.5)}
lgrid <- c(0.001, 0.01, 0.1, 1)

dgp_vdl <- function(T = 400, m = 10, d = 10){  

  # storage
  A  <- Y <- integer(T)
  W1 <- W3 <- integer(T)
  W2 <- integer(T)              # categorical 1,2,3
  lpY <- numeric(T)
  lpA <- numeric(T)
  
  # coefficients (Simulation 1a, vdL 2021)
  coef_A <- list(W1 =  0.25, W2 = -0.20, Y = 0.30,  A = -0.20, W3_2 = 0.20)
  coef_Y <- list(W1 = -0.80, W2 =  0.10, W3 = 0.20, A =  1.00,
                 W1_2 = -0.50, W3_2 = 0.20, intercept = 0.30)

  # Initialize first values that are all correlated with subsequent
  init <- m 

  A[1:init]  <- rbinom(init, 1, 0.5)
  lpA[1:init] <- 0
  Y[1:init]  <- yerr(init)
  lpY[1:init] <- 0

  W1[1:init] <- rbinom(init, 1, 0.5)
  W2[1:init] <- sample(1:3, init, TRUE)
  W3[1:init] <- rbinom(init, 1, 0.5) 
  
  for (t in (init + 1):T) {
    # treatment
    lpA[t] <- 0
    lpA[t] <- lpA[t] + coef_A$W1*sum(W1[(t-m):(t-1)]) + coef_A$W2*sum(W2[(t-m):(t-1)]) +
        coef_A$Y*sum(Y[(t-m):(t-1)])  + coef_A$A*sum(A[(t-m):(t-1)]) + coef_A$W3_2 * sum(W3[(t-m):(t-1)])
    A[t] <- rbinom(1, 1, expit(lpA[t] / (m)))
    
    # outcomes
    lpY[t] <- coef_Y$intercept + coef_Y$A * A[t]
    lpY[t] <- lpY[t] + coef_Y$W1*sum(W1[(t-m):(t-1)]) + coef_Y$W2*sum(W2[(t-m):(t-1)]) + coef_Y$W3*sum(W3[(t-m):(t-1)]) + coef_Y$W1_2*sum(W1[(t-m):(t-1)]) + coef_Y$W3_2*sum(W3[(t-m):(t-1)])
    Y[t] <- lpY[t] + yerr(1)
    
    # exogenous contemporaneous W
    W1[t] <- rbinom(1,1,0.5)
    W2[t] <- sample(1:3,1)
    W3[t] <- rbinom(1,1,0.5)
  }
  
  # finite lags
  DT <- data.table(t = 1:T, A, Y, W1, W2, W3)
  for (h in 1:d) {
    DT[, paste0("A_lag",h)  := shift(A,  h)]
    DT[, paste0("Y_lag",h)  := shift(Y,  h)]
    DT[, paste0("W1_lag",h) := shift(W1, h)]
    DT[, paste0("W2_lag",h) := shift(W2, h)]
    DT[, paste0("W3_lag",h) := shift(W3, h)]
  }
  setnafill(DT, type = "const", fill = 0)
  
  list(Y = DT$Y,
       A = DT$A,
       X = as.matrix(DT[, !c("t","A","Y")]),
       m = m,
       g = expit(lpA / (m)),
       Q1 = lpY)
}

########################################################## Chernozhukov Crossfit
ate_score <- function(Y, A, g, Q0, Q1, tau){
  g <- pmin(pmax(g, 1e-3), 1 - 1e-3)
  (A - g)/(g*(1-g)) * (Y - ifelse(A == 1, Q1, Q0) + Q1 - Q0 - tau) + (Q1 - Q0 - tau)
}

crossfit_iid_ts <- function(data, K = 2, seed = 258){
  set.seed(seed)
  n <- nrow(data$X)
  folds <- sample(rep(seq_len(K), length.out = n))
  g_hat <- Q0_hat <- Q1_hat <- numeric(n)
  
  for (k in seq_len(K)){
    tr <- folds != k; te <- !tr
    fit_g <- cv.glmnet(data$X[tr,,drop = FALSE], data$A[tr], family = "binomial", lambda = lgrid)
    fit_Q <- cv.glmnet(cbind(A = data$A[tr], data$X[tr,,drop = FALSE]), data$Y[tr], family = "gaussian", lambda = lgrid)
    g_hat[te] <- drop(predict(fit_g, data$X[te,, drop = FALSE], type = "response", s = "lambda.min"))
    Q1_hat[te] <- drop(predict(fit_Q, cbind(A=1, data$X[te,, drop = FALSE]), type = "response", s= "lambda.min"))
    Q0_hat[te] <- drop(predict(fit_Q, cbind(A = 0, data$X[te,,drop = FALSE]), s= "lambda.min"))
  }
  
  psi <- ate_score(data$Y, data$A, g_hat, Q0_hat, Q1_hat, 0)
  c(theta = mean(psi), se = sd(psi)/sqrt(n))
}

############################################################# Semenova Crossfit
make_nlo_folds <- function(T, K, gap) {
  block  <- floor(T / K)
  starts <- seq(1, by = block, length.out = K)
  
  lapply(seq_len(K), function(k) {
    te <- starts[k] : min(T, starts[k] + block - 1)
    
    gap_idx <- seq(
      max(1,  min(te) - gap),
      min(T, max(te) + gap)    
    )
    
    tr <- setdiff(seq_len(T), union(te, gap_idx))    
    list(tr = tr, te = te)
  })
}

crossfit_nlo_ts <- function(data, K = 5, seed = 258) {
  set.seed(seed)
  T  <- nrow(data$X)
  fd <- make_nlo_folds(T, K, gap = data$m)
  
  psi_list <- vector("list", K)
  for (k in seq_len(K)) {
    tr <- fd[[k]]$tr          # not $train
    te <- fd[[k]]$te          # not $test
    
    fit_g <- cv.glmnet(data$X[tr,, drop = FALSE], data$A[tr], family = "binomial", lambda = lgrid)
    fit_Q <- cv.glmnet(cbind(A = data$A[tr], data$X[tr,, drop = FALSE]),
                       data$Y[tr], family = "gaussian", lambda = lgrid)

    g_te  <- drop(predict(fit_g, data$X[te,, drop = FALSE],
                          type = "response", s = "lambda.min"))
    Q1_te <- drop(predict(fit_Q, cbind(A = 1, data$X[te,, drop = FALSE]),
                          type = "response",s = "lambda.min"))
    Q0_te <- drop(predict(fit_Q, cbind(A = 0, data$X[te,, drop = FALSE]),
                          type = "response",s = "lambda.min"))
    
    psi_list[[k]] <- ate_score(data$Y[te], data$A[te],
                               g_te, Q0_te, Q1_te, 0)
  }
  
  psi_all <- unlist(psi_list)
  c(theta = mean(psi_all), se = sd(psi_all) / sqrt(T))
}

crossfit_none_ts <- function(data, seed = 258) {
  set.seed(seed)
  T  <- nrow(data$X)

  fit_g <- cv.glmnet(data$X, data$A, family = "binomial", lambda = lgrid)
  fit_Q <- cv.glmnet(cbind(A = data$A, data$X), data$Y, family = "gaussian", lambda = lgrid)

    
  g_te  <- drop(predict(fit_g, data$X,
                        type = "response", s = "lambda.min"))
  Q1_te <- drop(predict(fit_Q, cbind(A = 1, data$X),
                        s = "lambda.min"))
  Q0_te <- drop(predict(fit_Q, cbind(A = 0, data$X),
                        s = "lambda.min"))

  psi_all <- ate_score(data$Y, data$A, g_te, Q0_te, Q1_te, 0)
    
  c(theta = mean(psi_all), se = sd(psi_all) / sqrt(T))
}


#################################################################### simulation!

perf_ts <- function(T_vec = c(400, 900, 1600, 2500), reps = 200, m = 10, K, seed = 42, true_tau = 0){
  set.seed(seed)
  res <- list(); idx <- 1
  for (T in T_vec){
    print(paste("Running T =", T))
    theta_iid <- theta_nlo <- theta_none <- time_iid <- time_nlo <- time_none <- numeric(reps)
    for (r in seq_len(reps)){
      data <- dgp_vdl(T, m, 0.5 * T)
      print(paste0("IID Replicate ", r))
      t1 <- system.time(est_iid <- crossfit_iid_ts(data, K, seed+T+r))["elapsed"]
      print(paste0("NLO Replicate ", r))
      t2 <- system.time(est_nlo <- crossfit_nlo_ts(data, K, seed+T+r))["elapsed"]
      print(paste0("No Cross-fit Replicate ", r))
      t3 <- system.time(est_none <- crossfit_none_ts(data, seed+T+r))["elapsed"]
      theta_iid[r] <- est_iid["theta"]; time_iid[r] <- t1
      theta_nlo[r] <- est_nlo["theta"]; time_nlo[r] <- t2
      theta_none[r] <- est_none["theta"]; time_none[r] <- t3
    }
    res[[idx]] <- data.table(method = "iid", T = T,
                             bias = mean(theta_iid - true_tau),
                             variance = var(theta_iid),
                             mse = mean((theta_iid-true_tau)^2),
                             mean_time = mean(time_iid)); idx <- idx+1
    res[[idx]] <- data.table(method = "nlo", T = T, 
                             bias = mean(theta_nlo - true_tau),
                             variance = var(theta_nlo),
                             mse = mean((theta_nlo - true_tau)^2),
                             mean_time = mean(time_nlo)); idx <- idx +1
    res[[idx]] <- data.table(method = "none", T = T, 
                             bias = mean(theta_none - true_tau),
                             variance = var(theta_none),
                             mse = mean((theta_none - true_tau)^2),
                             mean_time = mean(time_none)); idx <- idx +1
  }
  rbindlist(res)
}

set.seed(2025)


res_dt <- perf_ts(T_vec = c(400, 900, 1600, 2500), reps = 10, m = 10, K = 5, true_tau = 1)
write.csv(res_dt, here("data", "results_ts_hard.csv"), row.names = FALSE)

p_bias <- ggplot(res_dt, aes(T, bias, color = method)) +
  geom_line() + geom_point() +
  labs(title = "Bias", x = "T", y = "Bias") +
  theme_minimal()

p_var <- ggplot(res_dt, aes(T, variance, color = method)) +
  geom_line() + geom_point() +
  labs(title = "Variance", x = "T", y = "Variance") +
  theme_minimal()

p_mse <- ggplot(res_dt, aes(T, mse, color = method)) +
  geom_line() + geom_point() +
  labs(title = "MSE", x = "T", y = "MSE") +
  theme_minimal()

p_time <- ggplot(res_dt, aes(T, mean_time, color = method)) +
  geom_line() + geom_point() +
  labs(title = "Mean Time", x = "T", y = "Seconds") +
  theme_minimal()

p_final <- (p_bias | p_var) / (p_mse | p_time) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(here("figures", "plots_ts_hard.png"), p_final, width = 8, height = 6)

