library(simChef)
source("R/dgps.R")
source("R/estimators.R")

dgp_list <- list(
  Linear = create_dgp(dgp_iid, .name = "iid", n = 2000),
  Nonlinear = create_dgp(dgp_nonlin, .name = "Nonlin", n = 4000, sigma = 1),
  TwoWayCluster = create_dgp(dgp_cluster, .name = "Clustered"),
  mDepTS = create_dgp(dgp_ts, .name = "TS_mdep")
)



method_list <- list(
  ATE_GLM = create_method(est_ate_glm, .name = "ATE_GLM"),
  ATE_XGB = create_method(est_ate_xgb, .name = "ATE_XGB"),
  TMLE_GLM = create_method(est_tmle_glm, .name = "TMLE_GLM"),
  TMLE_XGB = create_method(est_tmle_xgb, .name = "TMLE_XGB"),
  DML_IID = create_method(dml_iid, .name = "DML_chernozhukov"),
  DML_CLUSTER = create_method(dml_2way, .name = "DML_chiang"),
  DML_TS_IID = create_method(dml_iid_ts, .name = "DML_iid_TS"),
  DML_TS_NLO = create_method(dml_nlo, .name = "DML_semenova")
)


eval_bias <- create_evaluator(function(fit_results, theta0 =1){
  data.table(bias = mean(vapply(fit_results, '[[', numeric(1), "est") - theta0))
},
.name = "Bias"
)


exp_cf <- create_experiment(name = "crossfit_dep", dgp_list = dgp_list, 
                            method_list = method_list, evaluator_list = list(eval_bias))








# Experiment 1: ATE DGPs × ATE methods ---------------------
exp_ate <- create_experiment(
  name         = "ate_only",
  dgp_list     = dgp_list[c("LAY", "Nonlin")],
  method_list  = method_list[c("ATE_GLM", "TMLE_GLM", "ATE_XGB", "TMLE_XGB")]
)

# Experiment 2: IV   DGPs × IV  methods ---------------------
exp_iv  <- create_experiment(
  name         = "iv_only",
  dgp_list     = dgp_list[c("Clustered", "TS_mdep")],
  method_list  = method_list[c("IV_IID", "IV_Cluster", "IV_TS_IID", "IV_TS_NLO")]
)


