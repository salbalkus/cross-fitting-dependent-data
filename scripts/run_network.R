
# Read in args from the command line
args = commandArgs(trailingOnly=TRUE)
n = as.numeric(args[1])     # First argument is the sample size
reps = as.numeric(args[2])  # Second argument is the number of replicates

library(here)
source(here("R", "new_network.R"))
print(paste("Running network simulation with n =", n, "and reps =", reps))
perf_network(truth, g_params, Q_params, n, reps = reps, K = 5, seed = 42)

