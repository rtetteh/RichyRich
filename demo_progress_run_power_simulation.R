# Example: Demonstrate progress bars in run_power_simulation

library(progressr)
library(furrr)

source(file.path(getwd(), "R", "moments_and_statistics.R"))
source(file.path(getwd(), "R", "regularization_selection.R"))
source(file.path(getwd(), "R", "power_simulation.R"))

# Set up progressr handler for console
handlers(global = TRUE)
handlers("txtprogressbar")

# Parameters chosen to run for ~2 minutes (adjust n_rep and mu_grid as needed)
n <- 50
n_sim <- 10
n_boot <- 20
alpha <- 0.05
moment_order <- 4
n_rep <- 2000  # Increase for longer runtime
mu_grid <- seq(0, 2, by = 0.2)  # 11 values
c_grid <- c(1, 0.1, 0.01)

progressr::with_progress({
  res <- run_power_simulation(
    n = n,
    n_sim = n_sim,
    n_boot = n_boot,
    alpha = alpha,
    moment_order = moment_order,
    n_rep = n_rep,
    mu_grid = mu_grid,
    c_grid = c_grid
  )
  print(res)
})
