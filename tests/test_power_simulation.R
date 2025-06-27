# Test suite for power_simulation.R
# These tests check that the power simulation runs and returns expected structure.

library(testthat)
library(purrr)
library(tibble)
library(MASS)
source("../R/moments_and_statistics.R")
source("../R/regularization_selection.R")
source("../R/power_simulation.R")

# Test: run_power_simulation returns a tibble with correct columns

test_that("run_power_simulation returns tibble with expected columns", {
  set.seed(123)
  res <- run_power_simulation(
    n = 10,
    n_sim = 2,
    n_boot = 5,
    alpha = 0.05,
    moment_order = 2,
    n_rep = 3,
    mu_grid = c(0, 1),
    c_grid = c(1, 0.1)
  )
  expect_s3_class(res, "tbl_df")
  expect_true(all(c("mu", "power_unweighted", "power_weighted") %in% colnames(res)))
  expect_equal(nrow(res), 2)
})

# Test: run_power_simulation power increases with mu (stochastic, but should trend)

test_that("run_power_simulation power increases with mu", {
  set.seed(123)
  res <- run_power_simulation(
    n = 20,
    n_sim = 3,
    n_boot = 10,
    alpha = 0.05,
    moment_order = 2,
    n_rep = 10,
    mu_grid = c(0, 2),
    c_grid = c(1, 0.1)
  )
  expect_true(res$power_unweighted[2] >= res$power_unweighted[1])
  expect_true(res$power_weighted[2] >= res$power_weighted[1])
})
