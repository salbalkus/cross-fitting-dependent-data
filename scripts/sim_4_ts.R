
library(MASS)
library(glmnet)
library(Matrix)

################################### m-dependent time series generation


# Y_t = \theta_0 D_t +X_t \beta_0 + \epsilon_t   ;   
# D_t (endogenous regressor)  = X_t \beta_0 + \nu_t
# Z_t (instrument) = X_t \beta_0 + V_t
# \epsilon, \nu, V are independent Moving Average ( MA(m) ) so the cov(u_t, u_{t+h}) = 0 for h > m

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

######################################################chernozhukov cross-fitting

crossfit_iid_ts <- function(df, K = 2){
  n <- nrow(df$X)
  f <- sample(rep(1:K, length.out = n))        #random fold
  
  psi_list <- vector("list", K)
  for (k in 1:K){
    test_idx <- which(f == k)
    train_idx <- which(f != k)
    
    gY <- cv.glmnet(df$X[train_idx, ], df$Y[train_idx], alpha = 1)
    gD <- cv.glmnet(df$X[train_idx, ], df$D[train_idx], alpha = 1)
    gZ <- cv.glmnet(df$X[train_idx, ], df$Z[train_idx], alpha = 1)
    
    gY_hat <- drop(predict(gY, df$X[test_idx, ]))
    gD_hat <- drop(predict(gD, df$X[test_idx, ]))
    gZ_hat <- drop(predict(gZ, df$X[test_idx, ]))
    
    psi_list[[k]] <- data.frame(psi = (df$Y[test_idx] - gY_hat) * (df$Z[test_idx] - gZ_hat),
                                u = df$D[test_idx] - gD_hat,
                                v = df$Z[test_idx] - gZ_hat
                                
    )
  }
  
  psi_all <- do.call(rbind, psi_list)
  
  J0 <- mean(psi_all$u * psi_all$v)
  theta_hat <- mean(psi_all$ps) / J0
  se <- sqrt(var(psi_all$psi - theta_hat * psi_all$u * psi_all$v) / (n * J0^2))
  
  c(theta = theta_hat, se = se)
}

############################################## neighbours left out cross-fitting

make_nlo_folds <- function(T, K, gap){
  block_len <- floor(T/K)
  starts <- seq(1, by = block_len, length.out = K)
  
  lapply(seq_len(K), function(k) {
    test_idx <- starts[k]:min(T, starts[k] + block_len - 1)
    # figuring out the gap:
    gap_idx <- seq(max(1, min(test_idx) - gap), min (T, max(test_idx) + gap))
    train_idx<- setdiff(seq_len(T) ,union(test_idx, gap_idx))
    
    list(train_idx = train_idx, test_idx = test_idx)
  })
}

crossfit_nlo_ts <- function(df, K = 5){
  T <- nrow(df$X)
  m <- df$m           # assuming you know the dependence span
  fd <- make_nlo_folds(T, K, gap = m)
  
  psi_list <- vector("list", K)
  
  for (k in seq_len(K)){
    train_idx <- fd[[k]]$train_idx
    test_idx <- fd[[k]]$test_idx
    
    gY <- cv.glmnet(df$X[train_idx, ], df$Y[train_idx], alpha = 1)
    gD <- cv.glmnet(df$X[train_idx, ], df$D[train_idx], alpha = 1)
    gZ <- cv.glmnet(df$X[train_idx, ], df$Z[train_idx], alpha = 1)
    
    gY_hat <- drop(predict(gY, df$X[test_idx, ]))
    gD_hat <- drop(predict(gD, df$X[test_idx, ]))
    gZ_hat <- drop(predict(gZ, df$X[test_idx, ]))
    
    psi_list[[k]] <- data.frame(psi = (df$Y[test_idx] - gY_hat) * (df$Z[test_idx] - gZ_hat),
                                u = df$D[test_idx] - gD_hat,
                                v = df$Z[test_idx] - gZ_hat
                                
    )
  }
  
  psi_all <- do.call(rbind, psi_list)
  
  J0 <- mean(psi_all$u * psi_all$v)
  theta_hat <- mean(psi_all$ps) / J0
  
  fold_means <- sapply(psi_list, function(df) mean(df$psi - theta_hat * df$u * df$v))
  
  se <- sqrt(mean(fold_means^2)/K * J0^2)
  
  c(theta = theta_hat, se = se)
  
}


################################################################# simulation run
sim_ts <- function(B = 25, T = 500, p = 100, m = 2, K_iid = 2, K_nlo = 5){
  res_iid <- matrix(NA, B, 2)
  res_nlo <- matrix(NA, B, 2)
  
  for (b in 1:B){
    df <- dgp_ts(T, p, m)
    res_iid[b, ] <- crossfit_iid_ts(df, K = K_iid)
    res_nlo[b, ] <- crossfit_nlo_ts(df, K = K_nlo)
  }
  
  list(iid = res_iid, nlo = res_nlo, theta0 = 1)
}


set.seed(42)

out <- sim_ts(B = 20, T = 500, p = 50, m = 2)

theta0 <- out$theta0

bias_iid <- mean(out$iid[,1] - theta0)
var_iid <- var(out$iid[,1])

bias_nlo <- mean(out$nlo[,1]) - theta0
var_nlo <- var(out$nlo[,1])


cat(sprintf("IID   bias = % .4f | var = %.4f\n", bias_iid, var_iid))
cat(sprintf("NLO(m) bias = % .4f | var = %.4f\n", bias_nlo, var_nlo))
