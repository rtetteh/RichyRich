#' Main script: Power analysis for moment-based test (tidyverse style)
#' 
#' Sources all modules and runs the simulation and plotting.

library(tidyverse)
library(MASS)
source('R/moments_and_statistics.R')
source('R/regularization_selection.R')
source('R/power_simulation.R')
source('R/plot_power_curves.R')

# Parameters
defaults <- list(
  n = 30,
  n_sim = 8,
  n_boot = 30,
  alpha = 0.05,
  moment_order = 5,
  n_rep = 50,
  mu_grid = seq(0, 2.1, by = 0.3),
  c_grid = c(50, 20, 10, 1e-1, 1e-3, 1e-6, 1e-9, 0)
)

# Run simulation
power_results <- run_power_simulation(
  n = defaults$n,
  n_sim = defaults$n_sim,
  n_boot = defaults$n_boot,
  alpha = defaults$alpha,
  moment_order = defaults$moment_order,
  n_rep = defaults$n_rep,
  mu_grid = defaults$mu_grid,
  c_grid = defaults$c_grid
)

# Print results
tibble::glimpse(power_results)

# Plot
plot_power_curves(power_results, defaults$alpha)
