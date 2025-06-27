# Test suite for regularization_selection.R
# These tests check that the R function select_optimal_regularization produces a valid beta and is numerically stable.

testthat::setup({
  source("../R/moments_and_statistics.R")
  source("../R/regularization_selection.R")
  source("../R/power_simulation.R")
  source("../R/plot_power_curves.R")
})

library(testthat)
library(purrr)
library(tibble)
library(MASS)
library(duckplyr)

# Test: select_optimal_regularization returns a numeric value in grid

test_that("select_optimal_regularization returns a value in c_grid", {
  set.seed(123)
  x <- rnorm(30)
  c_grid <- c(10, 1, 0.1, 0.001, 0)
  beta <- select_optimal_regularization(
    n = 30,
    c_grid = c_grid,
    moment_order = 3,
    n_sim = 4,
    n_boot = 10,
    x = x
  )
  expect_true(is.numeric(beta))
  expect_true(any(abs(beta - c_grid / (30^(1/3))) < 1e-8))
})

# Test: select_optimal_regularization is stable for repeated runs

test_that("select_optimal_regularization is stable for repeated runs", {
  set.seed(123)
  x <- rnorm(30)
  c_grid <- c(10, 1, 0.1, 0.001, 0)
  beta1 <- select_optimal_regularization(30, c_grid, 3, 4, 10, x)
  beta2 <- select_optimal_regularization(30, c_grid, 3, 4, 10, x)
  expect_equal(beta1, beta2, tolerance = 1e-10)
})

# Edge case: empty c_grid

test_that("select_optimal_regularization errors for empty c_grid", {
  set.seed(123)
  x <- rnorm(30)
  expect_error(select_optimal_regularization(30, numeric(0), 3, 4, 10, x))
})

# Edge case: single-value c_grid

test_that("select_optimal_regularization works for single-value c_grid", {
  set.seed(123)
  x <- rnorm(30)
  c_grid <- 0.5
  beta <- select_optimal_regularization(30, c_grid, 3, 4, 10, x)
  expect_true(is.numeric(beta))
  expect_equal(beta, c_grid / (30^(1/3)), tolerance = 1e-8)
})

# Edge case: all-equal c_grid

test_that("select_optimal_regularization works for all-equal c_grid", {
  set.seed(123)
  x <- rnorm(30)
  c_grid <- rep(2, 5)
  beta <- select_optimal_regularization(30, c_grid, 3, 4, 10, x)
  expect_true(is.numeric(beta))
  expect_true(any(abs(beta - c_grid[1] / (30^(1/3))) < 1e-8))
})

# Edge case: NA/Inf in c_grid

test_that("select_optimal_regularization errors for NA/Inf in c_grid", {
  set.seed(123)
  x <- rnorm(30)
  expect_error(select_optimal_regularization(30, c(1, NA, 0.1), 3, 4, 10, x))
  expect_error(select_optimal_regularization(30, c(1, Inf, 0.1), 3, 4, 10, x))
})

# Edge case: NA/Inf in x

test_that("select_optimal_regularization errors for NA/Inf in x", {
  c_grid <- c(10, 1, 0.1)
  expect_error(select_optimal_regularization(30, c_grid, 3, 4, 10, c(1, 2, NA)))
  expect_error(select_optimal_regularization(30, c_grid, 3, 4, 10, c(1, 2, Inf)))
})

# Edge case: non-numeric x or c_grid

test_that("select_optimal_regularization errors for non-numeric x or c_grid", {
  expect_error(select_optimal_regularization(30, c("a", "b"), 3, 4, 10, rnorm(30)))
  expect_error(select_optimal_regularization(30, c(1, 2), 3, 4, 10, c("a", "b")))
})

# Edge case: small n, n_sim, n_boot

test_that("select_optimal_regularization works for small n, n_sim, n_boot", {
  set.seed(123)
  x <- rnorm(2)
  c_grid <- c(1, 0.1)
  beta <- select_optimal_regularization(2, c_grid, 2, 1, 1, x)
  expect_true(is.numeric(beta))
})

# Edge case: output is always scalar, finite, non-NA

test_that("select_optimal_regularization output is scalar, finite, non-NA", {
  set.seed(123)
  x <- rnorm(30)
  c_grid <- c(10, 1, 0.1)
  beta <- select_optimal_regularization(30, c_grid, 3, 4, 10, x)
  expect_length(beta, 1)
  expect_true(is.finite(beta))
  expect_false(is.na(beta))
})

# Edge case: negative n or invalid input

test_that("select_optimal_regularization errors for negative n or invalid input", {
  expect_error(select_optimal_regularization(-1, c(1, 2), 3, 4, 10, rnorm(30)))
  expect_error(select_optimal_regularization(30, c(1, 2), 0, 4, 10, rnorm(30)))
})

# Test: duckplyr handles large c_grid
test_that("duckplyr handles large c_grid for select_optimal_regularization", {
  set.seed(123)
  x <- rnorm(30)
  c_grid <- as.numeric(1:1e5)
  c_grid_tbl <- as_duckplyr_df(data.frame(c_grid = c_grid))
  # Just test that function can handle the vector extracted from duckplyr_df
  beta <- select_optimal_regularization(30, c_grid_tbl$c_grid, 3, 4, 10, x)
  expect_true(is.numeric(beta))
})
