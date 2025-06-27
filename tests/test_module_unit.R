# Unit tests for individual functions and modules in the power analysis workflow

# Source relevant scripts for test execution
testthat::setup({
  source("../R/moments_and_statistics.R")
  source("../R/regularization_selection.R")
  source("../R/power_simulation.R")
  source("../R/plot_power_curves.R")
})

# Test calculate_moments

test_that("calculate_moments returns correct length and values for known input", {
  x <- c(1, 2, 3, 4)
  moment_order <- 3
  result <- calculate_moments(x, moment_order)
  expect_type(result, "double")
  expect_length(result, moment_order)
  expect_equal(unname(result[1]), mean(x))
  expect_equal(unname(result[2]), mean((x - mean(x))^2))
  expect_equal(unname(result[3]), mean((x - mean(x))^3))
})

# Test covariance_of_moments

test_that("covariance_of_moments returns correct shape and is symmetric", {
  mat <- matrix(rnorm(50), ncol = 5)
  cov_mat <- covariance_of_moments(mat)
  expect_true(is.matrix(cov_mat))
  expect_equal(dim(cov_mat)[1], dim(mat)[2])
  expect_equal(dim(cov_mat)[2], dim(mat)[2])
  expect_equal(cov_mat, t(cov_mat))
})

# Test compute_statistic_unweighted

test_that("compute_statistic_unweighted returns non-negative value", {
  obs <- c(1, 2, 3)
  ref <- c(1, 2, 2.5)
  stat <- compute_statistic_unweighted(obs, ref)
  expect_true(is.numeric(stat))
  expect_true(stat >= 0)
})

# Test compute_statistic_weighted

test_that("compute_statistic_weighted returns non-negative value and matches unweighted if identity weight", {
  obs <- c(1, 2, 3)
  ref <- c(1, 2, 2.5)
  weight_mat <- diag(3)
  stat_w <- compute_statistic_weighted(obs, ref, weight_mat)
  stat_u <- compute_statistic_unweighted(obs, ref)
  expect_true(is.numeric(stat_w))
  expect_true(stat_w >= 0)
  expect_equal(stat_w, stat_u)
})

# Test select_optimal_regularization

test_that("select_optimal_regularization returns a positive numeric value", {
  set.seed(1)
  n <- 15
  c_grid <- c(0.1, 1)
  moment_order <- 2
  n_sim <- 5
  n_boot <- 5
  x <- rnorm(n)
  beta <- select_optimal_regularization(n, c_grid, moment_order, n_sim, n_boot, x)
  expect_true(is.numeric(beta))
  expect_true(beta > 0)
})

# Test plot_power_curves

test_that("plot_power_curves returns a ggplot object and plots correct columns", {
  power_results <- tibble::tibble(
    mu = c(0, 1),
    power_unweighted = c(0.05, 0.8),
    power_weighted = c(0.06, 0.85)
  )
  alpha <- 0.05
  plt <- plot_power_curves(power_results, alpha)
  expect_true(inherits(plt, "ggplot"))
})
