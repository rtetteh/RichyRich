#' @title RichyRich Main Script
#' @description Example script to run the full power analysis workflow and plot results.
#' @examples
#' source('R/main.R')

library(tidyverse)
library(furrr)
library(MASS)

source('R/moment_utils.R')
source('R/test_statistics.R')
source('R/beta_selection.R')
source('R/power_analysis.R')
source('R/plotting.R')

plan(multisession)

results <- run_power_analysis(
  n = 30,
  n_sim = 8,
  n_boot = 30,
  n_rep = 100, # Reduce for demo, increase for real runs
  mu_values = seq(0, 2.1, 0.3),
  c_values = c(50, 20, 10, 1e-1, 1e-3, 1e-6, 1e-9, 0),
  moment_order = 5,
  alpha = 0.05
)

plot_power_results(results, alpha = 0.05)
