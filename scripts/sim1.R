## Packages

#             first time installs only:
# install.packages(c("devtools", "tmle3", "origami", "SuperLearner", "data.table", "ggplot2"))

library(devtools)
#install_github("tlverse/sl3@devel")
#install_github("tlverse/tmle3")

#install.packages("ranger")
#install.packages("xgboost")


library(sl3)   #super learner
library(tmle3)
library(origami) #for cross-validation/fitting 
library(ranger)
library(xgboost)

library(data.table)
library(stats)
library(ggplot2)




#setup/data-generation process

set.seed(258)


dgp_WAY <- function(n, a0 = 0, a1 =1, b0 = 0, b1 = 1, b2 = 1, sigma = 1){
  #confounder:
  W <- rnorm(n)
  
  #treatment: 
  p <- plogis(a0 + a1 * W)
  A <- rbinom(n, 1, prob = p)
  
  #outcome: 
  Y <- rnorm(n, mean = b0 + b1*A + b2*W, sd = sigma)
  
  
  data.table(W,A,Y)
}

true_ATE = 1



#estimation (one-step)

ATE_helper <- function(Y, A, g_hat, Q0_hat, Q1_hat){   #prints out the estimate with a 95% confidence interval
  IC <- (A/ g_hat) * (Y - Q1_hat) - ((1-A)/(1-g_hat)) * (Y - Q0_hat) + (Q1_hat - Q0_hat)
  psi <- mean(IC)
  se <- sd(IC)/sqrt(length(IC))
  ci <- psi + c(-1,1) * qnorm(0.975) * se
  list(ate = psi, se = se, ci = ci)
}



ATE_crossfit <- function(df, K = 5L){
  n <- nrow(df)
  g_hat <- Q0_hat <- Q1_hat <- numeric(n)
  
  #create folds of (roughly) equal sizes 
  fold_id <- sample(rep(1:K, length.out = n))
  
  
  for (k in 1:K){
    idx_tr <- which(fold_id != k)
    idx_te <- which(fold_id == k)
    
    #calculate propensity score with a logistic GLM
    fit_g <- glm (A ~ W, data = df[idx_tr], family = binomial)
    g_hat[idx_te] <- predict(fit_g, newdata = df[idx_te], type = "response")
    
    
    #calculate outgome with a regular LM
    fit_Q <- lm(Y~A + W, data = df[idx_tr])
    
    
    #predict for both treatment assignments
    te1 <- copy(df[idx_te]); te1[, A:=1]
    te0 <- copy(df[idx_te]); te0[, A:=0]
    
    Q1_hat[idx_te] <- predict(fit_Q, newdata = te1)
    Q0_hat[idx_te] <- predict(fit_Q, newdata = te0)
  }
  
  ATE_helper(df$Y, df$A, g_hat, Q0_hat, Q1_hat)
}



###example sim

n <- 2000
sim <- dgp_WAY(n)

res <- ATE_crossfit(sim, K = 5)
print(res)




#estimation (TMLE)  **yields wider confidence interval on same data**

nodes <- list(W = "W", A = "A", Y = "Y")

ATE_spec = tmle_ATE(treatment_level = 1, control_level = 0)

lrnr_glm <- make_learner(Lrnr_glm)

lrnr_list <- list(Y = lrnr_glm,
                  A = lrnr_glm)

tmle_fit <- tmle3(ATE_spec, sim, nodes, lrnr_list)
print(tmle_fit$summary)
