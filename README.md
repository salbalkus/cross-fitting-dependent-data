[ ] use sl3 or another selector to compare xgboost and other models for the ATEs

[ ] fit estimators on one of the datasets suggested by nima

  [ ] look at what papers try to simulate, maybe simplify, create dgp and estimation process that mimic paper. if not just find uses in other papers
  
[ ] take iid set up and generate clusters; draw a set of clusters and within each cluster draw units -- each cluster has a different set of parameters for this model. correlated covariates and correlated treaments, but not correlated outcomes. 

generate latent variablet hat affects treatment but not the outcome

[ ] ,,,for time series, at each time point, as long as you condition the previous treatment, the errors are independent -- uncontrollable factors lead to dependence in the error. same as clusters but L is an m-dep time series. 
