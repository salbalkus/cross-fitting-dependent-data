  
devtools::install_github("tlverse/tlverse")
suppressPackageStartupMessages({
  library(data.table)
  library(origami)      # folding utilities for cross-fitting
  library(glmnet)       # Lasso / Elastic-Net
  library(xgboost)      # gradient boosted trees
  library(tmle3)        # TMLE framework
  library(sl3)          # SuperLearner wrappers
})

ATE_helper <- function(Y, A, g_hat, Q0_hat, Q1_hat) {
  # One-step influence function & Wald CI
  IC  <- (A / g_hat) * (Y - Q1_hat) -
    ((1 - A) / (1 - g_hat)) * (Y - Q0_hat) +
    (Q1_hat - Q0_hat)
  psi <- mean(IC)
  se  <- sd(IC) / sqrt(length(IC))
  list(est = psi, se = se)
}

clip01 <- function(x, eps = 1e-4) pmin(pmax(x, eps), 1 - eps)
feat_cols <-setdiff(names(data), c("A", "Y"))


######################################### GLM one-step with K-fold cross-fitting
# est_ate_glm <- function(data, K = 5L, seed = 258) {
#   set.seed(seed)
#   n      <- nrow(data)
#   g_hat  <- Q0_hat <- Q1_hat <- numeric(n)
#   folds  <- sample(rep(seq_len(K), length.out = n))
#   
#   for (k in seq_len(K)) {
#     tr <- folds != k
#     te <- !tr
#     
#     # Propensity
#     fit_g <- glm(A ~ ., data = data[tr, c("A", feat_cols), with = FALSE], family = binomial)
#     g_hat[te] <- predict(fit_g, data[te, .(W)], type = "response")
#     
#     # Outcome
#     fit_Q <- glm(c("A", feat_cols), response = "Y", data = data[tr, c("A", "Y", feat_cols), with = FALSE])
#     data_te1 <- copy(data[te]); data_te1[, A := 1]
#     data_te0 <- copy(data[te]); data_te0[, A := 0]
#     Q1_hat[te] <- predict(fit_Q, data_te1)
#     Q0_hat[te] <- predict(fit_Q, data_te0)
#   }
#   
#   ATE_helper(data$Y, data$A, g_hat, Q0_hat, Q1_hat)
# }


# ---------------------------------------------------------------- 1. GLM one-step
est_ate_glm <- function(data, K = 5L, seed = 1) {
  set.seed(seed)
  
  ## -------- helpers
  covars <- setdiff(names(data), c("A", "Y"))          # all W-columns
  if (length(covars) == 0)
    stop("est_ate_glm: no covariate columns found (expect names other than A,Y)")
  
  # build formulas programmatically
  form_g <- as.formula(
    paste("A ~", paste(covars, collapse = " + "))
  )
  form_Q <- as.formula(
    paste("Y ~ A +", paste(covars, collapse = " + "))
  )
  
  ## -------- K-fold cross-fitting
  n      <- nrow(data)
  g_hat  <- Q0_hat <- Q1_hat <- numeric(n)
  fold_id <- sample(rep(seq_len(K), length.out = n))
  
  for (k in seq_len(K)) {
    idx_tr <- fold_id != k
    idx_te <- !idx_tr
    
    # Propensity
    fit_g <- glm(form_g, data = data[idx_tr], family = binomial)
    g_hat[idx_te] <- predict(fit_g, data[idx_te], type = "response")
    
    # Outcome
    fit_Q <- glm(form_Q, data = data[idx_tr])
    
    d1 <- copy(data[idx_te]); d1[, A := 1]
    d0 <- copy(data[idx_te]); d0[, A := 0]
    Q1_hat[idx_te] <- predict(fit_Q, d1)
    Q0_hat[idx_te] <- predict(fit_Q, d0)
  }
  
  ATE_helper(data$Y, data$A, g_hat, Q0_hat, Q1_hat)
}


########################################## XGBoost one-step (for non-linear DGP)
est_ate_xgb <- function(data, K = 5L, nrounds_g = 250, nrounds_Q = 250, seed = 258){
  set.seed(seed)
  
  n <- nrow(data)
  
  g_hat <- Q0_hat < Q1_hat <- numeric(n)
  
  folds <- make_folds(data, fold_fun = folds_vfold, V = K) # using origami
  
  for (k in seq_len(K)){
    idx_tr <- folds[[k]]$training_set
    idx_te <- folds[[k]]$validation_set
    
    
    # Propensity score estimation
    
    dg <- xgb.Dmatrix(as.matrix(data[idx_tr, ..feat_cols]), label = data$A[idx_tr])
    mod_g <- xgb.train(list(objective = "binary:logistic",
                            eval_metric = "logloss",
                            eta = 0.05, max_depth = 3),             #should tune at some point
                          dg, nrounds = nrounds_g, verbose = 0)
    g_hat[idx_te] <- clip01(predict(mod_g,
                             as.matrix(data[te, .(W1, W2)])))
    
    
    
    # Outcome regression estimation
    
    dQ <- xgb.Dmatrix(as.matrix(data[idx_tr, c(A, feat_cols), with = FALSE]), label = data$Y[tr])
    mod_Q <- xgb.train(list(objective = "reg:squarederror",
                            eval_metric = "rmse",
                            eta = 0.05, max_depth = 4),
                       dQ, nrounds = nrounds_Q, verbose = 0)
    
    
    ## predict for both treatment levels
    df_te1 <- copy(data[idx_te]); df_te1[, A:= 1]
    df_te0 <- copy(data[idx_te]); df_te0[, A:= 0 ]
    
    Q1_hat[te] <- predict(mod_Q, newdata = as.matrix(df_te1[, .(A, W1, W2)]))
    Q0_hat[te] <- predict(mod_Q, newdata = as.matrix(df_te0[, .(A, W1, W2)]))
  }
  
  ATE_helper(data$Y, data$A, g_hat, Q0_hat, Q1_hat)
}



######################################### TMLE with GLM learners (as a baseline)
est_tmle_glm <- function(data){
  nodes <-list(W = grep("^W", names(data), value = TRUE), A = "A", Y = "Y")
  
  lrnr <- make_learner(Lrnr_glm)
  tmle3(tmle_ATE(1,0), data, nodes, list(Y = lrnr, A = lrnr)) |> 
    (\(fit) list(est = fit$summary$psi_transofrmed, se = fit$summary$se))()
}



##################################################### TMLE with XGboost learners
est_tmle_xgb <- function(data){
  nodes <- list(W = c("W1", "W2"), A = "A", Y = "Y")
  lrnr_bin <- Lrnr_xgboost$new(objective = "binary:logistic",
                               eval_metric = "logloss",
                               nrounds = 250, eta = 0.05, max_depth = 3)
  
  lrnr_reg <- Lrnr_xgboost$new(objective = "reg:squarederror",
                               metric = "rmse",
                               nrounds = 250, eta = 0.05, max_depth = 4)
  tmle3(tmle_ATE(1,0), data, nodes, list(Y = lrnr_reg, A = lrnr_bin)) |>
    (\(fit) list(est = fit$summary$psi_transformed, se = fit$summary$se))()
}




#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#! DML


################################################## IID (Cherozhukov et al. 2018)
dml_iid <- function(data, K = 2L, seed = 258){
  set.seed(seed)
  n <- nrow(data$X)
  folds <- sample(rep(seq_len(K), length.out = n))
  
  collect <- vector("list", K)
  
  for (k in seq_len(K)){
    te <- folds == k; tr <- !idx_te
    
    gY <- cv.glmnet(data$X[tr, ], data$Y[tr], alpha = 1)
    gD <- cv.glmnet(data$X[tr, ], data$G[tr], alpha = 1)
    gZ <- cv.glmnet(data$X[tr, ], data$Z[tr], alpha = 1)
    
    collect[[k]] <- data.frame(
      psi = (data$Y[te] - drop(predict(gY, data$X[te, ]))) *
            (data$Z[te] - drop(predict(gZ, data$X[te, ]))),
      u = data$D[te] - drop(predict(gD, data$X[te, ])),
      v = data$Z[te] - drop(predict(gZ, data$X[te, ]))
    )
    
    psi <- do.call(rbind, collect)
    J0 <- mean(psi$u * psi$v)
    th <- mean(psi$psi) / J0
    se <- sqrt(var(psi$psi - th* psi$u * psi$v) / (n * J0^2))
    
    list(est = th, se = se)
    
  }
  
}



############################################# 2-way cluster (Chiang et al. 2023)
dml_2way <- function(data, K = 2, seed = 258){
  set.seed(seed)
  
  N <- max(data$idx$i); M <- max(data$idx$j)
  I <- split(sample(N), rep(seq_len(K), length.out = N))
  J <- split(sample(M), rep(seq_len(K), length.out = M))
  
  bag <- list()
  
  for (k1 in seq_len(K)) for (k2 in seq_len(K)){
    in_cells <- which(data$idx$i %in% I[[k1]] &
                        data$idx$j %in% I[[k2]])
    
    out_cells <- setdiff(seq_len(nrow(data$X)), in_cells)
    
    
    gY <- cv.glmnet(data$X[out_cells, ], data$Y[out_cells], alpha = 1)
    gD <- cv.glmnet(data$X[out_cells, ], data$D[out_cells], alpha = 1)
    gZ <- cv.glmnet(data$X[out_cells, ], data$Z[out_cells], alpha = 1)
    
    y_tilde <- data$Y[in_cells] - drop(predict(gY, data$X[in_cells, ]))
    d_tilde <- data$D[in_cells] - drop(predict(gD, data$X[in_cells, ]))
    z_tilde <- data$Z[in_cells] - drop(predict(gZ, data$X[in_cells, ]))
    
    bag[[length(bag) + 1]] <- data.frame(psi = y_tilde * z_tilde,
                                          u = d_tilde, v = z_tilde,
                                         i = data$idx$i[in_cells],
                                         j = data$idx$j[in_cells]
                                         )
  }
  psi <- do.call(rbind, bag)
  J0 <- mean(psi$u * psi$v)
  th <- mean(psi$psi) / J0
  
  gammaI <- tapply(psi$psi * psi$v, psi$i, mean)
  gammaJ <- tapply(psi$psi, psi$v, psi$j, mean)
  muI <- min(N< M) / N; muJ <- min(N,M)/M
  Gamma <- muI*var(gammaI) + muJ * var(gammaJ)
  se <- sqrt(Gamma/nrow(psi)* J0^2)
  list(est = th, se = se)
}



################################################ IID for m-dependent time series
dml_iid_ts <- function(data, K = 2L, seed = 258){
  set.seed(seed)
  n <- nrow(data$X)
  folds <- sample(rep(seq_len(K), length.out = n))
  
  bag <- vector("list", K)
  for (k in seq_len(K)){
    te <- folds == k; tr <- !te
    gY <- cv.glmnet(data$X[tr, ], data$Y[tr], alpha = 1)
    gD <- cv.glmnet(data$X[tr, ], data$D[tr], alpha = 1)
    gZ <- cv.glmnet(data$X[tr, ], data$Z[tr], alpha = 1)
    
    bag[[k]] <- data.frame(
      psi  = (data$Y[te] - drop(predict(gY, data$X[te, ])))*
              (data$Z[te] - drop(predict(gZ, data$X[te, ]))),
      u = data$D[te] - drop(predict(gZ, data$X[te, ])),
      v = data$Z[te] - drop(predict(gZ, data$X[te, ]))
    )
  }
  psi <- do.call(rbind, bag)
  J0 <- mean(psi$u * psi$v)
  th <- mean(psi$psi) /J0
  se <- sqrt(var(psi$psi - th * psi$u * psi$v) / (n * J0^2))
  list(est = th, se = se)
}




########################################### Neighbors left out (NLO)
make_nlo_folds <- function(T, K, gap){
  block <- floor(T/K)
  stats <- seq(1, by = block, length.out = K)
  
  lapply(seq_len(K), function(k){
    test <- starts[k] : min(T, starts[k] + block - 1)
    gap <- seq(max(1, min(test) - gap), min(T, max(test) + gap))
    list(train = setdiff(seq_len(T), union(test, gap)), test = test)
  }
    )
}
dml_nlo <- function(data, K = 5L, seed = 258){
  set.seed(seed)
  
  T <- nrow(data$X)
  f <- make_nlo_folds(T, K, gap = data$m)
  
  bag <- vector("list", K)
  for (k in seq_len(K)){
    tr <- fd[[k]]$train; te <- fd[[k]]$test
    gY <- cv.glmnet(data$X[tr, ], data$Y[tr], alpha = 1)
    gD <- cv.glmnet(data$X[tr, ], data$D[tr], alpha = 1)
    gZ <- cv.glmnet(data$X[tr, ], data$Z[tr], alpha = 1)
    
    bag[[k]] <- data.frame(
      psi = (data$Y[te] - drop(predict(gY, data$X[te, ]))) *
        data$Z[te] - drop(predict(gZ, data$X[te, ])),
      u = data$D[te] - drop(predict(gD, data$X[te, ])),
      v = data$Z[te] - drop(predict(gZ, data$X[te, ]))
    )
  }
  psi <- do.call(rbind, bag)
  Kbar <- mean(sapply(bag, nrow))
  J0 <- mean(psi$u * psi$v)
  th <- mean(psi$psi) / J0
  fold_means <- sapply(bag, \(df) mean (df$psi - th * df$u * df$v))
  se <- sqrt(mean(fold_means^2)/Kbar*J0^2)
  list(est = th, se = se)
}

