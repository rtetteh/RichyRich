library(testthat)
library(tidyverse)

source('../R/moment_utils.R')
source('../R/test_statistics.R')
source('../R/beta_selection.R')

# Test compute_moments

test_that('compute_moments returns correct value', {
  expect_equal(compute_moments(c(1, 2, 3), 2), mean(c(1, 2, 3)^2))
  expect_equal(compute_moments(rep(2, 5), 3), 8)
})

# Test compute_kernel

test_that('compute_kernel returns correct shape and values', {
  mat <- matrix(1:20, nrow = 4, ncol = 5)
  kernel <- compute_kernel(mat, 4)
  expect_equal(dim(kernel), c(5, 5))
  expect_true(is.matrix(kernel))
})

# Test calculate_test_statistic_1

test_that('calculate_test_statistic_1 returns squared distance', {
  expect_equal(calculate_test_statistic_1(c(1, 2), c(1, 1)), 1)
  expect_equal(calculate_test_statistic_1(c(0, 0), c(0, 0)), 0)
})

# Test calculate_test_statistic_1opt

test_that('calculate_test_statistic_1opt returns correct value', {
  psi <- c(1, 2)
  psi_bar <- c(1, 1)
  K_inv_half <- diag(2)
  expect_equal(calculate_test_statistic_1opt(psi, psi_bar, K_inv_half), 1)
})

# Test select_optimal_beta (basic smoke test)
test_that('select_optimal_beta returns a numeric value', {
  set.seed(1)
  beta <- select_optimal_beta(10, c(1, 0.1), 2, 2, 2, rnorm(10))
  expect_true(is.numeric(beta))
})
