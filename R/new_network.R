library(ggplot2)
library(igraph)
library(gridExtra)
library(xgboost)


### Data Generating Process
dgp_network <- function(n, scale=5) {
  # Confounders
  W1 <- rbeta(n, 2, 2)
  W2 <- rpois(n, 10)
  W3 <- rbinom(n, 1, 0.3)
  mu <- 5*(W1 > 0.4) - 2*(W1 > 0.6) + 3*(W1 > 0.7) + (2*W3 - 1)*W2 + 20
  p_treat  <- plogis(mu/20 - 1)
  A  <- rbinom(n, 1, p_treat)
  mu_treat <- (2*A + 1) * mu 
  G <- as_adjacency_matrix(sample_gnp(n, 3 / n))
  cor_err <- (G %*% (scale * rbeta(n, 6, 6)))[,1]
  Y  <- mu_treat + cor_err
  list(W1 = W1, W2 = W2, W3 = W3, A = A, Y = Y, p_treat = p_treat, mu = mu, mu_treat = mu_treat, G = G)
}

### Score Function for estimand
ate_score <- function(Y, A, g, Q0, Q1, tau){
  g <- pmin(pmax(g, 1e-6), 1 - 1e-6) # truncate to avoid numerical issues
  (Q1 - Q0) + (2*A - 1) / g * (Y - A * Q1 - (1 - A) * Q0)
}

xgb_train <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
xgb_test <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)


### Compute an estimate using Cross-fitting from Emmenegger et al. (2024)
crossfit_emm_net <- function(data, g_params, Q_params, K = 2, seed = 258){
  n <- nrow(data$X)
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
crossfit_iid_net <- function(data, g_params, Q_params, K = 2, seed = 258){
  n <- nrow(data$X)
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

data_for_truth = dgp_network(1000000)
truth = mean(data_for_truth$mu * (2*data_for_truth$A + 1))
