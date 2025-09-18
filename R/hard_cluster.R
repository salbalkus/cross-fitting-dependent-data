library(earth)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(here)

dgp_two_way <- function(N, M, sigma_Y, seed = NULL){
  ## row and column effects (i.i.d. standard Ns)
  a.X <- rnorm(N); b.X <- rnorm(M)     # covariates
  a.m <- rnorm(N); b.m <- rnorm(M)     # propensiy
  a.t <- rnorm(N); b.t <- rnorm(M)     # treatment effect (it'll be heterogenous)
  a.0 <- rnorm(N); b.0 <- rnorm(M)     # baseline
  
  
  grid <- CJ(i = 1:N, j = 1:M)
  
  # Covariates
  X1 <- a.X[grid$i] + b.X[grid$j] + rnorm(N*M)
  X2 <- rnorm(N*M)
  X3 <- sin(X1) + 0.3*X2 + rnorm(N*M, 0, 0.5)
  X4 <- (X1 > 0)*X2 + rnorm(N*M, 0, 0.5)
  X5 <- rnorm(N*M, 1, 2)
  
  lin_ps <- -0.4 + 0.8*X1 -0.7*X2^2 +0.5*sin(X3) + 0.4*X1*X2 - 0.5*X4*X5 + 0.6*a.m[grid$i] - 0.6*b.m[grid$j]
  
  q_true <- plogis(lin_ps/5)
  A <- rbinom(N*M, 1, q_true)
  
  
  # treatment effect
  tau_true <- 1 + 0.5 * sin(X1) - 0.5 * sin(X2) + 0.3*X3*X4 + 0.4*a.t[grid$i]*b.t[grid$j]
  
  
  # baseline outcome
  Q0 <- 2 + 0.5 * X1 - 0.4 * X2 + 0.3 * X3^2 - 0.5*sin(X4) + 0.4* X1 * X2 +
    0.6*a.0[grid$i] + 0.6*b.0[grid$j]
  
  Y <- Q0 + A * tau_true + rnorm(N*M, 0, sigma_Y)
  
  dt <- data.table(i = grid$i, j = grid$j, Y = Y, A = A,
                   X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5)
  
  attr(dt, "ATE_true") <- mean(tau_true)
  dt
}

##################################################################### helpers
featurize <- function(dt){
  X <- as.matrix(dt[, .(X1, X2, X3, X4, X5)])
  cbind(X,
        X[,1]^2, X[,2]^2, X[,3]^2,
        sin(X[,1]), sin(X[,2]), sin(X[,3]),
        X[,1] * X[,2], X[,2] * X[,3], X[,3]*X[,4], X[,4]*X[,5])
  
}

make_multi_folds <- function(N, M, K){
  list(rf = sample(rep(1:K, length.out = N)),
       cf = sample(rep(1:K, length.out = M)))
}

################################################################## iid crossfit

iid_cf <- function(dt, K = 2){
  n <- nrow(dt)
  cv <- sample(rep(1:K, length.out = n))
  dt[, ':='(m_hat = NA_real_, Q1_hat = NA_real_, Q0_hat = NA_real_)]
  
  for (k in 1:K){
    tr <- cv != k; te <- cv == k
    #Xtr <- featurize(dt[tr]); Xte <- featurize(dt[te])
    Xtr <- dt[tr]; Xte <- dt[te]
    
    g_fit <- earth(A ~ X1 + X2 + X3 + X4 + X5, data = Xtr, degree = 2, glm = list(family="binomial"))
    dt$g_hat[te] <- predict(g_fit, Xte, type = "response")
    
    y_fit <- earth(Y ~ A + X1 + X2 + X3 + X4 + X5, data = Xtr, degree = 2, glm = list(family="gaussian"))

    Xte1 <- cbind(subset(Xte, select=-c(A)), A = 1)
    dt$Q1_hat[te] <- predict(y_fit, Xte1)
    Xte0 <- cbind(subset(Xte, select=-c(A)), A = 0)
    dt$Q0_hat[te] <- predict(y_fit, Xte0)
  }
  
  
  g <- dt$g_hat; Q1 <- dt$Q1_hat; Q0 <- dt$Q0_hat
  g <- pmin(pmax(g, 1e-3), 1 - 1e-3) # Trim propensity scores to avoid extreme weights

  A <- dt$A; Y <- dt$Y
  dr <- A*(Y-Q1)/g - (1-A)*(Y-Q0)/(1-g) + (Q1 - Q0)
  theta <- mean(dr); psi <- dr - theta
  list(theta = theta, se = sqrt(mean(psi^2))/sqrt(n), psi = psi)
  
}

############################################################### 2-way crossfit

twoWay_cf <- function(dt, K = 2){
  N <- max(dt$i); M <- max(dt$j)
  f <- make_multi_folds(N, M, K)
  dt[, ':='(rf = f$rf[i], cf = f$cf[j], g_hat = NA_real_, Q0_hat = NA_real_, Q1_hat = NA_real_ )]
  
  for (r in 1:K) for (c in 1:K){
    test <- dt$rf == r & dt$cf == c
    train <- !(dt$rf == r | dt$cf == c)  # strict complement!! no rows and columns!
    
    #Xtr <- featurize(dt[train]); Xte <- featurize(dt[test])
    Xtr <- dt[train]; Xte <- dt[test]
    
    g_fit <- earth(A ~ X1 + X2 + X3 + X4 + X5, data = Xtr, degree = 2, glm = list(family="binomial"))
    dt$g_hat[test] <- predict(g_fit, Xte, type = "response")
    
    y_fit <- earth(Y ~ A + X1 + X2 + X3 + X4 + X5, data = Xtr, degree = 2, glm = list(family="gaussian"))

    Xte1 <- cbind(subset(Xte, select=-c(A)), A = 1)
    dt$Q1_hat[test] <- predict(y_fit, Xte1)
    Xte0 <- cbind(subset(Xte, select=-c(A)), A = 0)
    dt$Q0_hat[test] <- predict(y_fit, Xte0)
  }
  
  g <- dt$g_hat; Q1 <- dt$Q1_hat; Q0 <- dt$Q0_hat
  g <- pmin(pmax(g, 1e-3), 1 - 1e-3) # Trim propensity scores to avoid extreme weights
  A <- dt$A; Y <- dt$Y
  dr <- A*(Y-Q1)/g - (1-A)*(Y-Q0)/(1-g) + (Q1 - Q0)
  theta <- mean(dr); psi <- dr - theta
  
  #row_sum <- dt[,.(s = sum(psi)), by = i]$s ; col_sum <- dt[,.(s = sum(psi)), by = j]$s
  #var_hat <- (sum(row_sum^2) + sum(col_sum^2) - sum(psi^2))/ (N*M)^2
  var_hat <- mean(psi^2)
  list(theta = theta, se = sqrt(var_hat), psi  = psi)
}

sim <- function(N_vec, M_vec, R, K_iid, K_2w){
  
  # Generate the large-sample truth
  dt_truth = dgp_two_way(N = 10000, M = 10000, sigma_Y = 0.5)
  truth = attr(dt_truth, "ATE_true")
  remove(dt_truth)

  out <- list(); i <- 1
  for (N in N_vec) for (M in M_vec) if (N == M){
    print(paste0("Running ", N, " x ", M))
    res <- replicate(R, {
      dt <- dgp_two_way(N,M, 0.5)
      #truth <- attr(dt, "ATE_true") #within-sample truth
      
      print(paste0("Simulating IID..."))
      t1 <- system.time(est_iid <- iid_cf(copy(dt), K_iid))[3]
      print(paste0("Simulating Cluster..."))
      t2 <- system.time(est_2w <- twoWay_cf(copy(dt), K_iid))[3]
      
      c(theta_iid = est_iid$theta,
        theta_2w = est_2w$theta,
        time_iid = t1,
        time_2w = t2,
        truth = truth)
    })
    
    res <- as.data.table(t(res))
    out[[i]] <- data.table(
      N = N, M = M,
      bias_iid = mean((res$theta_iid - res$truth) / res$truth),
      var_iid = var(res$theta_iid),
      mse_iid = mean((res$theta_iid - res$truth)^2),
      mean_t_iid = mean(res$time_iid),
      
      bias_2w = mean((res$theta_2w - res$truth) / res$truth),
      var_2w = var(res$theta_2w),
      mse_2w = mean((res$theta_2w - res$truth)^2),
      mean_t_2w = mean(res$time_2w)
    )
    i <- i+1
  }
  rbindlist(out)
}

### Run Simulation ###
set.seed(1)

tab <- sim(N_vec = c(30,40,50,60), M_vec = c(30,40,50,60), R = 10, K_iid = 2, K_2w = 2)

write.csv(tab, here("data", "results_cluster_hard.csv"), row.names = FALSE)

nm_levels <- tab %>% arrange(N, M) %>% 
  transmute(NM = paste0(N, "×", M)) %>% distinct() %>% pull()

plot_df <- tab %>% 
  transmute(NM = factor(paste0(N, "×", M), levels = nm_levels),
            bias_iid,  bias_2w,
            var_iid,   var_2w,
            mse_iid,   mse_2w,
            time_iid = mean_t_iid,
            time_2w  = mean_t_2w) %>%
  pivot_longer(-NM,
               names_to = c("metric", "estimator"),
               names_sep = "_",
               values_to = "value") %>%
  mutate(metric    = factor(metric,
                            levels = c("bias","var","mse","time"),
                            labels = c("Bias","Variance","MSE","CPU time (s)")),
         estimator = recode(estimator,
                            iid  = "IID DML",
                            tw   = "2‑Way DML"))

write.csv(plot_df, here("data", "plotdf_cluster_hard.csv"), row.names = FALSE)

ggplot(plot_df,
       aes(x = NM, y = value,
           colour = estimator, group = estimator)) +
  geom_line(position = position_dodge(width = .5)) +
  geom_point(position = position_dodge(width = .5), size = 3) +
  facet_wrap(~metric, ncol = 2, scales = "free_y") +
  labs(x = "Grid size  (N × M)",
       y = NULL,
       colour = NULL,
       title = "IID vs TwoWay DML across grid sizes") +
  theme_bw() +
  theme(strip.text  = element_text(face = "bold"),
        plot.title  = element_text(size = 14, face = "bold", hjust = .5),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here("figures", "fig_cluster_hard.png"), width = 8, height = 6)
