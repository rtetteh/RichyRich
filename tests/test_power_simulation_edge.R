# Additional tests for power simulation: happy path and edge cases

# Source relevant scripts for test execution
testthat::setup({
  source("../R/moments_and_statistics.R")
  source("../R/regularization_selection.R")
  source("../R/power_simulation.R")
  source("../R/plot_power_curves.R")
})

test_that("run_power_simulation handles single mu and single c_grid", {
  n <- 20
  n_sim <- 10
  n_boot <- 10
  alpha <- 0.05
  moment_order <- 2
  n_rep <- 5
  mu_grid <- 0.5
  c_grid <- 0.1

  result <- run_power_simulation(
    n = n,
    n_sim = n_sim,
    n_boot = n_boot,
    alpha = alpha,
    moment_order = moment_order,
    n_rep = n_rep,
    mu_grid = mu_grid,
    c_grid = c_grid
  )

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1)
  expect_true(result$power_unweighted >= 0 && result$power_unweighted <= 1)
  expect_true(result$power_weighted >= 0 && result$power_weighted <= 1)
})

test_that("run_power_simulation errors on invalid input (edge cases)", {
  expect_error(run_power_simulation(
    n = 0, n_sim = 10, n_boot = 10, alpha = 0.05, moment_order = 2, n_rep = 5, mu_grid = 0.5, c_grid = 0.1
  ), "n must be a single positive integer")

  expect_error(run_power_simulation(
    n = 10, n_sim = 0, n_boot = 10, alpha = 0.05, moment_order = 2, n_rep = 5, mu_grid = 0.5, c_grid = 0.1
  ), "n_sim must be a single positive integer")

  expect_error(run_power_simulation(
    n = 10, n_sim = 10, n_boot = 0, alpha = 0.05, moment_order = 2, n_rep = 5, mu_grid = 0.5, c_grid = 0.1
  ), "n_boot must be a single positive integer")

  expect_error(run_power_simulation(
    n = 10, n_sim = 10, n_boot = 10, alpha = 1.5, moment_order = 2, n_rep = 5, mu_grid = 0.5, c_grid = 0.1
  ), "alpha must be a single numeric value between 0 and 1")

  expect_error(run_power_simulation(
    n = 10, n_sim = 10, n_boot = 10, alpha = 0.05, moment_order = 0, n_rep = 5, mu_grid = 0.5, c_grid = 0.1
  ), "moment_order must be a single positive integer")

  expect_error(run_power_simulation(
    n = 10, n_sim = 10, n_boot = 10, alpha = 0.05, moment_order = 2, n_rep = 0, mu_grid = 0.5, c_grid = 0.1
  ), "n_rep must be a single positive integer")

  expect_error(run_power_simulation(
    n = 10, n_sim = 10, n_boot = 10, alpha = 0.05, moment_order = 2, n_rep = 5, mu_grid = numeric(0), c_grid = 0.1
  ), "mu_grid must be a non-empty numeric vector")

  expect_error(run_power_simulation(
    n = 10, n_sim = 10, n_boot = 10, alpha = 0.05, moment_order = 2, n_rep = 5, mu_grid = 0.5, c_grid = numeric(0)
  ), "c_grid must be a non-empty numeric vector")
})
