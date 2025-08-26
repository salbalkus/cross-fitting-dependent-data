library(ggplot2)
library(igraph)
library(gridExtra)
library(xgboost)
library(data.table)


### Data Generating Process
dgp_network <- function(n, scale=1, degree_mean = 3) {
  # Confounders
  W1 <- rbeta(n, 2, 2)
  W2 <- rpois(n, 10)
  W3 <- rbinom(n, 1, 0.3)
  mu <- 5*(W1 > 0.4) - 2*(W1 > 0.6) + 3*(W1 > 0.7) + (2*W3 - 1)*W2 + 20
  p_treat  <- plogis(mu/20 - 1)
  A  <- rbinom(n, 1, p_treat)
  mu_treat <- (2*A + 1) * mu 
  G <- as_adjacency_matrix(sample_gnp(n, degree_mean / n))
  cor_err <- (G %*% (scale * 2*(rbeta(n, 6, 6) - 0.5)))[,1] # Subtract 0.5 to make error mean-0
  Y  <- mu_treat + cor_err
  list(W1 = W1, W2 = W2, W3 = W3, A = A, Y = Y, p_treat = p_treat, mu = mu, mu_treat = mu_treat, cor_err = cor_err, G = G)
}

### Score Function for estimand
ate_score <- function(Y, A, g, Q0, Q1, tau){
  g <- pmin(pmax(g, 1e-6), 1 - 1e-6) # truncate to avoid numerical issues
  (Q1 - Q0) + (2*A - 1) / g * (Y - A * Q1 - (1 - A) * Q0)
}


### Compute an estimate using Cross-fitting from Emmenegger et al. (2024)
crossfit_emm_net <- function(data, g_params, Q_params, K = 5, seed = 258){
  n <- length(data[[1]])
  folds <- sample(rep(seq_len(K), length.out = n))
  g_hat <- Q0_hat <- Q1_hat <- numeric(n)
  X <- data.frame(W1 = data$W1, W2 = data$W2, W3 = data$W3)
  A <- data$A
  Y <- data$Y
  
  D = (diag(rep(1, n)) + data$G + data$G %*% data$G) > 0
  for (k in seq_len(K)){
    # Split into training and testing
    # Throw out the training observations correlated with test observations

    te <- folds == k
    tr <- !te
    tr <- which(!mapply(any, asplit(D[tr, te], 1)))
    te <- which(te)

    # Fit models
    g_train <- xgb.DMatrix(data = as.matrix(X[tr,,]), label = A[tr])
    g_test <- as.matrix(X[te,,])

    Q_train <- xgb.DMatrix(data = as.matrix(cbind(X, A)[tr,,]), label = Y[tr])
    Q1_test <- as.matrix(cbind(X, A = rep(1, n))[te,,])
    Q0_test <- as.matrix(cbind(X, A = rep(0, n))[te,,])

    fit_g <- xgboost(data = g_train, params = g_params, nrounds = 100, verbose = 0)
    fit_Q <- xgboost(data = Q_train, params = Q_params, nrounds = 100, verbose = 0)

    g_hat[te] <- predict(fit_g, g_test)
    Q1_hat[te] <- predict(fit_Q, Q1_test)
    Q0_hat[te] <- predict(fit_Q, Q0_test)
  }
  
  psi <- ate_score(data$Y, data$A, g_hat, Q0_hat, Q1_hat, 0)
  c(theta = mean(psi), se = sd(psi)/sqrt(n))
}

# Implement as-IID cross-fitting
crossfit_iid_net <- function(data, g_params, Q_params, K = 5, seed = 258){
  n <- length(data[[1]])
  folds <- sample(rep(seq_len(K), length.out = n))
  g_hat <- Q0_hat <- Q1_hat <- numeric(n)
  X <- data.frame(W1 = data$W1, W2 = data$W2, W3 = data$W3)
  A <- data$A
  Y <- data$Y
  
  for (k in seq_len(K)){
    # Split into training and testing
    tr <- folds != k; te <- !tr

    # Fit models
    g_train <- xgb.DMatrix(data = as.matrix(X[tr,,]), label = A[tr])
    g_test <- as.matrix(X[te,,])

    Q_train <- xgb.DMatrix(data = as.matrix(cbind(X, A)[tr,,]), label = Y[tr])
    Q1_test <- as.matrix(cbind(X, A = rep(1, n))[te,,])
    Q0_test <- as.matrix(cbind(X, A = rep(0, n))[te,,])

    fit_g <- xgboost(data = g_train, params = g_params, nrounds = 100, verbose = 0)
    fit_Q <- xgboost(data = Q_train, params = Q_params, nrounds = 100, verbose = 0)

    g_hat[te] <- predict(fit_g, g_test)
    Q1_hat[te] <- predict(fit_Q, Q1_test)
    Q0_hat[te] <- predict(fit_Q, Q0_test)
  }
  
  psi <- ate_score(data$Y, data$A, g_hat, Q0_hat, Q1_hat, 0)
  c(theta = mean(psi), se = sd(psi)/sqrt(n))
}

g_params <- list(
  booster = "gbtree",
  eta = 0.1,
  objective = "reg:logistic"
)

Q_params <- list(
  booster = "gbtree",
  eta = 0.1,
  objective = "reg:squarederror"
)

# Compute true ATE
data_for_truth = dgp_network(1000000)
truth = mean(2*data_for_truth$mu) # When A = 1, mu_treat = 3*mu; when A = 0, mu_treat = mu. The difference is 2*mu.

### Run simulations

perf_network <- function(truth, g_params, Q_params, n_vec = c(400, 900, 1600, 2500), reps = 1, K = 5, seed = 42){
  set.seed(seed)
  res <- list(); idx <- 1
  # Iterate through each sample size
  for (n in n_vec){
    se_iid <- theta_iid <- time_iid <- se_emm <- theta_emm <- time_emm <- numeric(reps)
    # Iterate through each replicate
    for (r in seq_len(reps)){
      data <- dgp_network(n)
      t1 <- system.time(est_iid <- crossfit_iid_net(data, g_params, Q_params, K = K, seed = seed+n+r))["elapsed"]
      t2 <- system.time(est_emm <- crossfit_emm_net(data, g_params, Q_params, K = K, seed = seed+n+r))["elapsed"]
      theta_iid[r] <- est_iid["theta"]; se_iid[r] <- est_iid["se"]; time_iid[r] <- t1
      theta_emm[r] <- est_emm["theta"]; se_emm[r] <- est_emm["se"]; time_emm[r] <- t2
    }
    res[[idx]] <- data.table(method = "iid", n = n,
                             bias = mean(theta_iid - truth),
                             variance = var(theta_iid),
                             mse = mean((theta_iid-truth)^2),
                             mean_est_sd = mean(se_iid),
                             mean_time = mean(time_iid)); idx <- idx+1
    res[[idx]] <- data.table(method = "emm", n = n, 
                             bias = mean(theta_emm - truth),
                             variance = var(theta_emm),
                             mse = mean((theta_emm - truth)^2),
                             mean_est_sd = mean(se_emm),
                             mean_time = mean(time_emm)); idx <- idx +1
  }
  rbindlist(res)
}

out <- perf_network(truth, g_params, Q_params, reps = 10)

res_dt <- out
library(patchwork)

p_bias <- ggplot(res_dt, aes(n, bias, color = method)) +
  geom_line() + geom_point() +
  labs(title = "Bias", x = "T", y = "Bias") +
  theme_minimal()

p_var <- ggplot(res_dt, aes(n, variance, color = method)) +
  geom_line() + geom_point() +
  labs(title = "Variance", x = "T", y = "Variance") +
  theme_minimal()

p_mse <- ggplot(res_dt, aes(n, mse, color = method)) +
  geom_line() + geom_point() +
  labs(title = "MSE", x = "T", y = "MSE") +
  theme_minimal()

p_time <- ggplot(res_dt, aes(n, mean_time, color = method)) +
  geom_line() + geom_point() +
  labs(title = "Mean Time", x = "T", y = "Seconds") +
  theme_minimal()

(p_bias | p_var) / (p_mse | p_time) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
