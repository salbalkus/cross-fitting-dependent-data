library(ggplot2)
library(dplyr)
library(gridExtra)
library(here)
library(patchwork)

read_data <- function(filename, version = 1, ns){
    dfs = lapply(ns, \(n) read.csv(here("data", paste0("sim_net", n, "_", version, ".csv"))))
    return(do.call(rbind, dfs))
}

create_plot <- function(filename, version, ns){

    full_data <- read_data(filename, version = version, ns = ns)

    res_dt <- full_data %>% mutate(bias = theta - truth) %>%
    group_by(method, n) %>%
    summarise(bias = mean(bias),
                variance = var(theta),
                mse = mean((theta - truth)^2),
                mean_time = mean(time))

    p_bias <- ggplot(res_dt, aes(n, bias, color = method)) +
    geom_line() + geom_point() +
    labs(title = "Bias", x = "T", y = "Bias") +
    theme_minimal()

    p_var <- ggplot(res_dt, aes(n, variance, color = method)) +
    geom_line() + geom_point() +
    labs(title = "Variance", x = "T", y = "Variance") +
    theme_minimal()

    p_mse <- ggplot(res_dt, aes(n, mse, color = method)) +
    geom_line() + geom_point() +
    labs(title = "MSE", x = "T", y = "MSE") +
    theme_minimal()

    p_time <- ggplot(res_dt, aes(n, mean_time, color = method)) +
    geom_line() + geom_point() +
    labs(title = "Mean Time", x = "T", y = "Seconds") +
    theme_minimal()

    (p_bias | p_var) / (p_mse | p_time) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}
