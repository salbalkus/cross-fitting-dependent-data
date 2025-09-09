renv::restore(prompt = FALSE)
renv::load(here::here())

library(here)
print("Running cluster simulation")
source(here("R", "new_cluster.R"))
