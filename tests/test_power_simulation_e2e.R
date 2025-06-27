# Test for end-to-end power analysis logic

# Source relevant scripts for test execution
testthat::setup({
  source("../R/moments_and_statistics.R")
  source("../R/regularization_selection.R")
  source("../R/power_simulation.R")
  source("../R/plot_power_curves.R")
})

test_that("run_power_simulation returns correct structure and power increases with mu", {
  # Arrange: Set up parameters for a small, fast test
  n <- 30
  n_sim <- 20
  n_boot <- 20
  alpha <- 0.05
  moment_order <- 2
  n_rep <- 10
  mu_grid <- c(0, 0.5, 1)
  c_grid <- c(0.1, 1)

  # Act: Run the power simulation
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

  # Assert: Structure is correct
  expect_s3_class(result, "tbl_df")
  expect_true(all(c("mu", "power_unweighted", "power_weighted") %in% colnames(result)))
  expect_equal(nrow(result), length(mu_grid))

  # Assert: Power should increase as mu moves away from 0
  expect_true(result$power_unweighted[1] <= result$power_unweighted[2] + 0.1)
  expect_true(result$power_unweighted[2] <= result$power_unweighted[3] + 0.1)
  expect_true(result$power_weighted[1] <= result$power_weighted[2] + 0.1)
  expect_true(result$power_weighted[2] <= result$power_weighted[3] + 0.1)

  # Assert: Power at mu = 0 is near alpha
  expect_true(abs(result$power_unweighted[1] - alpha) < 0.15)
  expect_true(abs(result$power_weighted[1] - alpha) < 0.15)
})
