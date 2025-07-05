# Cross-Fitting on Dependent Data

The goal of this project is to provide a guide to practitioners on how to perform cross-fitting when conducting causal inference in settings where units are correlated with each other (i.e. clustered, spatial, network, longitudinal data). The final product will include:

1. A literature review surveying how cross-fitting has been used in causal machine learning, and what its purpose is
2. The main theorem showing that as-iid cross-fitting is still valid in the network setting
3. Specializations of the main theorem to several practical data settings, with simulations demonstrating the results.

## To do
- [X] Set up recurring meeting room in the Smith Campus center on Tuesdays or Thursdays.
- [X] Literature review. Gather sources on cross-fitting and/or causal inference in dependent data (clustered, longitudinal, network, spatial, etc.), and add them to refs.bib on Overleaf with one sentence describing what they do and one sentence on why they're relevant to our project.
- [X] Set up simulations using the `here` package in R. Generate data from $L$ - $A$ - $Y$, and use machine learning model to learn the ATE. Use "Causal Inference: What if?" as a reference, and put your code in the "scripts" folder.
- [X] Implement nonlinear simulation and compute one-step and TMLE with a machine learning model like XGBoost or earth/MARS
- [X] Look at the dependent data papers from the literature and create simulation settings with similar dependency structures
  - [X] clusters
  - [X] time-series
  - [ ] networks
- [X] Implement different variations of cross-fitting strategies from the papers in the literature review and Lahiri (2003) "Resampling Methods for Dependent Data".
  - [X] clusters
  - [X] time-series
  - [ ] networks
- [ ] Run simulations that record the bias and variance of the estimator across different sample sizes using different sample splitting strategies.
- [X] Derive corollaries based on Theorem 1 for specific dependency structures. For example, a clustered data structure should admit a CLT with rate $\sqrt{n}$.
