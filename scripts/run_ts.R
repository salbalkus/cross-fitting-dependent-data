renv::restore(prompt = FALSE)
renv::load(here::here())

library(here)
print("Running time series simulation")
source(here("R", "new_ts.R"))
