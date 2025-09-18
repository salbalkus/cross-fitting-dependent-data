library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(patchwork)
library(knitr)
library(kableExtra)
library(stringr)
library(here)

# Cluster results
df_cluster <- read.csv(here::here("data", "results_cluster_hard.csv"))

df_cluster_long <- df_cluster %>%
    pivot_longer(
        cols = matches("(_iid|_2w)$"),
        names_to = c(".value", "type"),
        names_pattern = "(.*)_(iid|2w)$"
    ) %>%
    mutate(N_M = paste0(N, " x ", M)) %>%
    select(-N, -M)
df_cluster_long$type = factor(df_cluster_long$type, levels = c("iid", "2w"), labels = c("As-independent", "Two-way"))
# Create table
df_cluster_tbl = df_cluster %>% 
    mutate(N_M = paste0(N, " x ", M)) %>% select(-N, -M) %>%
    select(N_M, bias_iid, bias_2w, var_iid, var_2w, mse_iid, mse_2w, mean_t_iid, mean_t_2w)

kable(df_cluster_tbl, format = "latex", col.names = c("Clusters", "IID", "Two-Way", "IID", "Two-Way", "IID", "Two-Way", "IID", "Two-Way"), booktabs = TRUE, digits = 3) %>% 
    add_header_above(c(" " = 1, "Bias (%)" = 2, "Variance" = 2, "MSE" = 2, "Mean Time (s)" = 2)) %>%
    kable_styling()

val_colors = c("As-independent" = "#002e5c", "Two-way" = "#C67700")

p_bias <- ggplot(df_cluster_long, aes(x = N_M, y = bias, color = type, group = type)) +
    geom_line() +
    geom_point() +
    geom_hline(aes(yintercept = 0)) +
    scale_color_manual(values = val_colors) +
    labs(y = "Bias (%)", x = "Sample Size") +
    theme_minimal() +
    theme(legend.position = "bottom") + guides(color = guide_legend(title = "Type"))

p_var <- ggplot(df_cluster_long, aes(x = N_M, y = var, color = type, group = type)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = val_colors) +
    labs(y = "Variance", x = "Sample Size") +
    theme_minimal() +
    theme(legend.position = "bottom") + guides(color = guide_legend(title = "Type"))

p_mse <- ggplot(df_cluster_long, aes(x = N_M, y = mse * c(900, 1600, 2500, 3600), color = type, group = type)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = val_colors) +
    labs(y = "Scaled MSE", x = "Sample Size") +
    theme_minimal() +
    theme(legend.position = "bottom") + guides(color = guide_legend(title = "Type"))

# Arrange plots in a grid with shared legend
plot_grid <- (p_bias + p_var + p_mse) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

ggsave(here("figures", "cluster.png"), plot = plot_grid, width = 9, height = 3)

# Time series results
df_ts <- read.csv(here::here("data", "results_ts_hard.csv"))
df_ts$method = factor(df_ts$method, levels = c("nlo", "iid"), labels = c("Neighbors left out", "As-independent"))

# Create table
df_ts_tbl = df_ts %>% pivot_wider(names_from = method, values_from = c(bias, variance, mse, mean_time))

kable(df_ts_tbl, format = "latex", col.names = c("Clusters", "IID", "NLO", "IID", "NLO", "IID", "NLO", "IID", "NLO"), booktabs = TRUE, digits = 3) %>% 
    add_header_above(c(" " = 1, "Bias (%)" = 2, "Variance" = 2, "MSE" = 2, "Mean Time (s)" = 2)) %>%
    kable_styling()

# Create plots

val_colors = c("As-independent" = "#002e5c", "Neighbors left out" = "#C67700")

p_bias <- ggplot(df_ts, aes(x = T, y = bias, color = method, group = method)) +
    geom_line() +
    geom_point() +
    geom_hline(aes(yintercept = 0)) +
    scale_x_continuous(breaks = c(400, 900, 1600, 2500)) +
    scale_color_manual(values = val_colors) +
    labs(y = "Bias (%)", x = "Sample Size") +
    theme_minimal() +
    theme(legend.position = "bottom") + guides(color = guide_legend(title = "Type"))

p_var <- ggplot(df_ts, aes(x = T, y = variance, color = method, group = method)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = c(400, 900, 1600, 2500)) +
    scale_color_manual(values = val_colors) +
    labs(y = "Variance", x = "Sample Size") +
    theme_minimal() +
    theme(legend.position = "bottom") + guides(color = guide_legend(title = "Type"))

p_mse <- ggplot(df_ts, aes(x = T, y = mse * T, color = method, group = method)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = c(400, 900, 1600, 2500)) +
    scale_color_manual(values = val_colors) +
    labs(y = "Scaled MSE", x = "Sample Size") +
    theme_minimal() +
    theme(legend.position = "bottom") + guides(color = guide_legend(title = "Type"))

# Arrange plots in a grid with shared legend
plot_grid <- (p_bias + p_var + p_mse) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

ggsave(here("figures", "time_series.png"), plot = plot_grid, width = 9, height = 3)



##### Network data #####

df_net_full = lapply(c(400, 900, 1600, 2500), function(n){read.csv(str_glue(here::here("data", "sim_net{n}_1.csv")))}) %>% bind_rows()

df_net = df_net_full %>% group_by(n, method) %>% 
    summarize(bias = mean(theta) - mean(truth), variance = var(theta), mse = mean((theta - truth)^2), mean_time = mean(time)) %>%
    mutate(method = factor(method, levels = c("emm", "iid"), labels = c("Emmenegger et al.", "As-independent")))

df_net_tbl =  df_net %>%
    pivot_wider(names_from = method, values_from = c(bias, variance, mse, mean_time))

kable(df_net_tbl, format = "latex", col.names = c("Clusters", "IID", "Net", "IID", "Net", "IID", "Net", "IID", "Net"), booktabs = TRUE, digits = 3) %>% 
    add_header_above(c(" " = 1, "Bias (%)" = 2, "Variance" = 2, "MSE" = 2, "Mean Time (s)" = 2)) %>%
    kable_styling()

# Create plots without legends and titles
val_colors = c("As-independent" = "#002e5c", "Emmenegger et al." = "#C67700")

p_bias <- ggplot(df_net, aes(x = n, y = bias, color = method, group = method)) +
    geom_line() +
    geom_point() +
    geom_hline(aes(yintercept = 0)) +
    scale_x_continuous(breaks = c(400, 900, 1600, 2500)) +
    scale_color_manual(values = val_colors) +
    labs(y = "Bias (%)", x = "Sample Size") +
    theme_minimal() +
    theme(legend.position = "bottom") + guides(color = guide_legend(title = "Type"))

p_var <- ggplot(df_net, aes(x = n, y = variance, color = method, group = method)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = c(400, 900, 1600, 2500)) +
    scale_color_manual(values = val_colors) +
    labs(y = "Variance", x = "Sample Size") +
    theme_minimal() +
    theme(legend.position = "bottom") + guides(color = guide_legend(title = "Type"))

p_mse <- ggplot(df_net, aes(x = n, y = mse * n / log(n), color = method, group = method)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = c(400, 900, 1600, 2500)) +
    scale_color_manual(values = val_colors) +
    labs(y = "Scaled MSE", x = "Sample Size") +
    theme_minimal() +
    theme(legend.position = "bottom") + guides(color = guide_legend(title = "Type"))

p_mean_t <- ggplot(df_net, aes(x = n, y = mean_time, color = method, group = method)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = c(400, 900, 1600, 2500)) +
    scale_color_manual(values = val_colors) +
    labs(y = "Mean Time (s)", x = "Sample Size") +
    theme_minimal() +
    theme(legend.position = "bottom") + guides(color = guide_legend(title = "Type"))

# Arrange plots in a grid with shared legend
plot_grid <- (p_bias + p_var + p_mse) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

ggsave(here("figures", "network.png"), plot = plot_grid, width = 9, height = 3)
