# Publication code for "On the use of cross-fitting in causal machine learning with correlated data"

This code is meant to reproduce the simulation results in Balkus, Laith, and Hejazi (2026). 

* ```R``` contains functions and other R code used for running experiments. 
* ```scripts``` contains code to generate publication-ready plots and to run code on the cluster. 
* ```data``` contains .csv files storing simulation results. 
* ```figures``` visualize the simulation results from ```data```, created by ```create_plots.R``` in the ```scripts``` folder. 

To reproduce results from the paper, rerun the each file in ```scripts```. Each ```run_``` script corresponds to one of the example DGPs in the paper (clustered, network, and time series data). The ```create_plots.R``` script produces the plots in the manuscript using the data generated from each ```run_``` script. 
