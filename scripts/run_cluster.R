renv::restore(prompt = FALSE)
renv::load(here::here())

library(here)
print("Running cluster simulation")
source(here("R", "hard_cluster.R"))
