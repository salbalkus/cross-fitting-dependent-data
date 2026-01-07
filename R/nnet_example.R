library(nnet)
library(data.table)
library(here)

expit <- function(x) 1/(1+exp(-x))
maxit <- 300

# Data generating process
dgp <- function(n, d){
    X <- matrix(rbeta(n * d, 2, 2), nrow = n, ncol = d)
    linear_model <- X[, 1] + X[, 2] + X[, 3] + X[, 4]
    g <- expit(linear_model)
    A <- rbinom(n, 1, g)
    Q <- (A + 2)*linear_model
    Y <- rnorm(n, mean = Q, sd = 1)
    list(X = X, A = A, Y = Y, g = g, Q = Q)
}

# ATE influence function
ate_score <- function(Y, A, g, Q0, Q1, tau){
  g <- pmin(pmax(g, 1e-2), 1 - 1e-2)
  (A - g)/(g*(1-g)) * (Y - ifelse(A == 1, Q1, Q0) + Q1 - Q0 - tau) + (Q1 - Q0 - tau)
}

# Neural network cross-fitting for IID data
crossfit_iid <- function(data, K = 5, seed = 258){
  set.seed(seed)
  n <- nrow(data$X)
  folds <- sample(rep(seq_len(K), length.out = n))
  g_hat <- Q0_hat <- Q1_hat <- Q_hat <- numeric(n)
  
  for (k in seq_len(K)){
    tr <- folds != k; te <- !tr

    fit_g <- nnet(data$X[tr,,drop = FALSE], data$A[tr], size = 5, maxit = maxit, decay = 0.1, trace = F)
    fit_Q <- nnet(cbind(A = data$A[tr], data$X[tr,,drop = FALSE]), data$Y[tr], linout = T, size = 5, decay = 0.1, maxit = maxit, trace = F)
    g_hat[te] <- predict(fit_g, data$X[te,, drop = FALSE])
    Q1_hat[te] <- predict(fit_Q, cbind(A=1, data$X[te,, drop = FALSE]), type = "raw")
    Q0_hat[te] <- predict(fit_Q, cbind(A=0, data$X[te,,drop = FALSE]), type = "raw")

  }
  
  psi <- ate_score(data$Y, data$A, g_hat, Q0_hat, Q1_hat, 0)
  c(theta = mean(psi), se = sd(psi)/sqrt(n))
}

# Neural network without cross-fitting for IID data
crossfit_none <- function(data, seed = 258) {
  set.seed(seed)
  T  <- nrow(data$X)

  fit_g <- nnet(data$X, data$A, size = 5, maxit = maxit, decay = 0.05, trace = F)
  fit_Q <- nnet(cbind(A = data$A, data$X), data$Y, linout = T, size = 5, decay = 0.05, maxit = maxit, trace = F)
  g_te <- predict(fit_g, data$X, type = "raw")
  Q1_te <- predict(fit_Q, cbind(A=1, data$X), type = "raw")
  Q0_te <- predict(fit_Q, cbind(A = 0, data$X), type = "raw")

  psi_all <- ate_score(data$Y, data$A, g_te, Q0_te, Q1_te, 0)
    
  c(theta = mean(psi_all), se = sd(psi_all) / sqrt(T))
}

perf_ts <- function(n_vec = c(400, 900, 1600, 2500), reps = 10, K = 5, d = 5, seed = 42, true_tau = 0){
  set.seed(seed)
  res <- list(); idx <- 1
  for (n in n_vec){
    print(paste("Running n =", n))
    theta_iid <- theta_none <- time_iid <- time_none <- numeric(reps)
    for (r in seq_len(reps)){
      data <- dgp(n, d)
      print(paste0("IID Replicate ", r))
      t1 <- system.time(est_iid <- crossfit_iid(data, K, seed+n+r))["elapsed"]
      print(paste0("No Cross-fit Replicate ", r))
      t3 <- system.time(est_none <- crossfit_none(data, seed+n+r))["elapsed"]
      theta_iid[r] <- est_iid["theta"]; time_iid[r] <- t1
      theta_none[r] <- est_none["theta"]; time_none[r] <- t3
    }
    res[[idx]] <- data.table(method = "iid", n = n,
                             bias = mean(theta_iid - true_tau),
                             variance = var(theta_iid),
                             mse = mean((theta_iid-true_tau)^2),
                             mean_time = mean(time_iid)); idx <- idx+1
    res[[idx]] <- data.table(method = "none", n = n,
                             bias = mean(theta_none - true_tau),
                             variance = var(theta_none),
                             mse = mean((theta_none - true_tau)^2),
                             mean_time = mean(time_none)); idx <- idx +1
  }
  rbindlist(res)
}

set.seed(2025)

d = 20
true_tau = 2
res_dt <- perf_ts(n_vec = c(400, 900, 1600, 2500), reps = 20, K = 5, d = d, true_tau = true_tau)

write.csv(res_dt, here("data", "nnet_example.csv"), row.names = FALSE)
