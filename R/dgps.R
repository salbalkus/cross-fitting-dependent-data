library(data.table)
library(stats)
library(MASS)
library(devtools)
library(sl3)   #super learner
library(tmle3)
library(origami) #for cross-validation/fitting 
library(xgboost)
library(ggplot2) 
library(igraph)


################################# linear, iid

dgp_iid <- function(n, a0 = 0, a1 =1, b0 = 0, b1 = 1, b2 = 1, sigma = 1){
  #confounder:
  W <- rnorm(n)
  
  #treatment: 
  p <- plogis(a0 + a1 * W)
  A <- rbinom(n, 1, prob = p)
  
  #outcome: 
  Y <- rnorm(n, mean = b0 + b1*A + b2*W, sd = sigma)
  
  
  data.table(W,A,Y)
}


################################# Non-linear 

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



################################# two-way cluster

dgp_cluster <- function(N = 25, M = 25, p = 100,
                        omega = c(.25, .25),
                        s_X = .25, s_ev = .25) {
  # following matrices and vectors are all straight from section 5.1 of Chiang 2019
  Sigma.p <- toeplitz(s_X^(0:(p-1)))
  
  # shocks for X
  aX  <- mvrnorm(N*M, rep(0, p), Sigma.p)
  aXi <- mvrnorm(N,   rep(0, p), Sigma.p)
  aXj <- mvrnorm(M,   rep(0, p), Sigma.p)
  
  # shocks for errors & instruments
  aEv <- mvrnorm(N*M, c(0,0), matrix(c(1,s_ev,s_ev,1), 2))
  aE  <- mvrnorm(N,   c(0,0), matrix(c(1,s_ev,s_ev,1), 2))
  aV  <- mvrnorm(M,   c(0,0), matrix(c(1,s_ev,s_ev,1), 2))
  
  idx <- expand.grid(i = 1:N, j = 1:M)
  
  make_2way <- function(cell, row_eff, col_eff) {
    sub_rows <- function(mat_or_vec, rows) {
      if (is.matrix(mat_or_vec))
        mat_or_vec[rows, , drop = FALSE]
      else
        mat_or_vec[rows]
    }
    row_part <- sub_rows(row_eff, idx$i)
    col_part <- sub_rows(col_eff, idx$j)
    
    (1 - sum(omega)) * cell +
      omega[1] * row_part +
      omega[2] * col_part
  }
  
  X   <- make_2way(aX,       aXi,     aXj)
  eps <- make_2way(aEv[,1],  aE[,1],  aV[,1])
  nu  <- make_2way(aEv[,2],  aE[,2],  aV[,2])
  V   <- make_2way(rnorm(N*M), rnorm(N), rnorm(M))
  
  theta0 <- 1
  beta0  <- rep(0.5, p)
  
  D <- X %*% beta0 + nu
  Z <- X %*% beta0 + V
  Y <- D*theta0    + X %*% beta0 + eps
  
  list(Y = Y, D = D, Z = Z, X = X, idx = idx, theta = theta0)
}



########################################## m-dependent time series

dgp_ts <- function(T = 600,     #series length
                   p = 30,      #number of covariates
                   m = 30,      #dependence range (lag)
                   s_X = 0.5,   #autodecay for feature covariance
                   theta0 = 1,
                   beta_0 = rep(0.5, p),
                   seed = 258
){
  Sigma.p <- as.matrix(toeplitz(s_X^(0:(p - 1))))
  stopifnot(length(rep(0, p)) == ncol(Sigma.p))   # mean length vs Σ
  stopifnot(all(dim(Sigma.p) == c(p, p)))
  
  X <- mvrnorm(T, rep(0,p), Sigma.p)   # correlated features
  
  MA_m <- function(n,m){
    u <- rnorm(n+m)
    as.numeric(stats::filter(u, rep(1, m+1), sides = 1)[(m+1):(n+m)])
  }
  
  nu <- MA_m(T,m)
  epsilon <- MA_m(T,m)
  V <- MA_m(T,m)
  
  D <- as.vector(X %*% beta_0 + nu)
  Z <- as.vector(X %*% beta_0 + V)
  Y <- theta0 * D + as.vector(X %*% beta_0) + epsilon
  
  idx <- data.frame(t = seq_len(T))                         #keep a time index
  
  list(Y = Y, D = D, Z = Z, X = X, idx = idx, theta = theta0, m = m)
}

################################# network setting

dgp_network <- function(n, scale=5) {
  # Confounders
  W1 <- rbeta(n, 2, 2)
  W2 <- rpois(n, 10)
  W3 <- rbinom(n, 1, 0.3)
  mu <- 5*(W1 > 0.4) - 2*(W1 > 0.6) + 3*(W1 > 0.7) + (2*W3 - 1)*W2 + 20
  p_treat  <- plogis(mu/20 - 1)
  A  <- rbinom(n, 1, p_treat)
  mu_treat <- (2*A + 1) * mu 
  G <- as_adjacency_matrix(sample_gnp(n, 10 / n))
  cor_err <- (G %*% (scale * rbeta(n, 6, 6)))[,1]
  Y  <- mu_treat + cor_err
  list(W1 = W1, W2 = W2, W3 = W3, A = A, Y = Y, p_treat = p_treat, mu_treat = mu_treat, G = G)
}


