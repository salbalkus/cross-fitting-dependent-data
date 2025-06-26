```{r}
# install.packages(c("MASS", "glmnet", "data.table", "progressr"))
# 
# library(MASS)
# library(glmnet)
# library(Matrix)                                                                
# library(progressr)
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



#########################################################chernozhukov cross-fit

crossfit_iid <- function(df, K = 2) {
  n <- nrow(df$X)
  f <- sample(rep(1:K, length.out = n))
  
  psi_list <- vector("list", K)
  for (k in 1:K) {
    test  <- which(f == k)
    train <- which(f != k)
    
    gY <- cv.glmnet(df$X[train, ], df$Y[train], alpha = 1)
    gD <- cv.glmnet(df$X[train, ], df$D[train], alpha = 1)
    gZ <- cv.glmnet(df$X[train, ], df$Z[train], alpha = 1)
    
    gYhat <- drop(predict(gY, df$X[test, ]))
    gDhat <- drop(predict(gD, df$X[test, ]))
    gZhat <- drop(predict(gZ, df$X[test, ]))
    
    psi_list[[k]] <- data.frame(
      psi = (df$Y[test] - gYhat) * (df$Z[test] - gZhat),
      u   =  df$D[test] - gDhat,
      v   =  df$Z[test] - gZhat
    )
  }
  
  psi_all <- do.call(rbind, psi_list)
  
  J0        <- mean(psi_all$u * psi_all$v)
  theta_hat <- mean(psi_all$psi) / J0
  se        <- sqrt( var(psi_all$psi - theta_hat * psi_all$u * psi_all$v) /
                       (n * J0^2) )
  
  c(theta = theta_hat, se = se)
}


################################################################chiang cross-fit

crossfit_cluster <- function(df, K = 2) {
  
  N <- max(df$idx$i);  M <- max(df$idx$j)
  
  # split rows/cols into K roughly equal folds
  I <- split(sample(1:N), rep(1:K, length.out = N))
  J <- split(sample(1:M), rep(1:K, length.out = M))
  
  psi_list <- list()
  
  for (k1 in 1:K) for (k2 in 1:K) {
    
    in_cells  <- which(df$idx$i %in% I[[k1]] &
                         df$idx$j %in% J[[k2]])
    out_cells <- setdiff(seq_len(nrow(df$X)), in_cells)
    
    # nuisance fits on out of fold cells
    gY <- cv.glmnet(df$X[out_cells, ], df$Y[out_cells], alpha = 1)
    gD <- cv.glmnet(df$X[out_cells, ], df$D[out_cells], alpha = 1)
    gZ <- cv.glmnet(df$X[out_cells, ], df$Z[out_cells], alpha = 1)
    
    # drop() to get vectors rather than 1 column matrices
    gYhat <- drop(predict(gY, df$X[in_cells, ]))
    gDhat <- drop(predict(gD, df$X[in_cells, ]))
    gZhat <- drop(predict(gZ, df$X[in_cells, ]))
    
    y_tilde <- df$Y[in_cells] - gYhat
    d_tilde <- df$D[in_cells] - gDhat
    z_tilde <- df$Z[in_cells] - gZhat
    
    psi_list[[length(psi_list) + 1]] <- data.frame(
      psi = y_tilde * z_tilde,       
      u   = d_tilde,
      v   = z_tilde,
      i   = df$idx$i[in_cells],
      j   = df$idx$j[in_cells]
    )
  }
  
  psi_all <- do.call(rbind, psi_list)   # clean data-frame, 5 columns
  
  
  J0        <- mean(psi_all$u * psi_all$v)
  theta_hat <- mean(psi_all$psi) / J0
  
  # two-way-cluster variance (Chiang et al. eqn 3.5)
  # helper to get cluster-level means of psi times v and psi times u within each cluster
  gammaN <- tapply(seq_along(psi_all$psi), psi_all$i,
                   function(idx) mean(psi_all$psi[idx] * psi_all$v[idx]))
  gammaM <- tapply(seq_along(psi_all$psi), psi_all$j,
                   function(idx) mean(psi_all$psi[idx] * psi_all$u[idx]))
  
  muN <- min(N, M) / N
  muM <- min(N, M) / M
  Gamma <- muN * var(gammaN) + muM * var(gammaM)
  
  se <- sqrt(Gamma / (nrow(psi_all) * J0^2))
  
  c(theta = theta_hat, se = se)
}


######################################################################sample run

sim <- function(B = 500, N = 25, M = 25, p = 100){
  res_iid <- matrix(NA, B, 2); res_cluster <- res_iid
  for (b in 1:B){
    df <- dgp_cluster(N, M, p)
    res_iid[b,] <- crossfit_iid(df)
    res_cluster[b,] <- crossfit_cluster(df)
  }
  list(iid = res_iid, cluster = res_cluster)
}


out <- sim(B = 20)

theta0 <- 1                                # true value in the DGP

bias_iid      <- mean(out$iid[,1])      - theta0
var_iid       <- var( out$iid[,1] )

bias_cluster  <- mean(out$cluster[,1])  - theta0
var_cluster   <- var( out$cluster[,1] )

cat("IID bias:", bias_iid, "   /   IID variance:", var_iid, "\n")
cat("K^2 bias:", bias_cluster, "   /   K^2 variance:", var_cluster)


```
