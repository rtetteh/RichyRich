# Test suite for moments_and_statistics.R
# These tests check that the R functions reproduce the same results as the original Python code.

library(testthat)
library(purrr)
source("../R/moments_and_statistics.R")

# Test: calculate_moments matches numpy mean for powers

test_that("calculate_moments matches numpy mean for powers", {
  set.seed(123)
  x <- rnorm(10)
  # Python: np.mean(x ** 1), np.mean(x ** 2), np.mean(x ** 3)
  expect_equal(
    unname(calculate_moments(x, 3, central = FALSE)),
    unname(c(mean(x), mean(x^2), mean(x^3))),
    tolerance = 1e-10
  )
})

# Test: compute_statistic_unweighted matches numpy dot

test_that("compute_statistic_unweighted matches numpy dot", {
  a <- c(1, 2, 3)
  b <- c(2, 2, 2)
  # Python: np.dot(a-b, a-b)
  expect_equal(
    compute_statistic_unweighted(a, b),
    sum((a-b)^2),
    tolerance = 1e-10
  )
})

# Test: compute_statistic_weighted matches Mahalanobis quadratic form

test_that("compute_statistic_weighted matches Mahalanobis quadratic form", {
  a <- c(1, 2, 3)
  b <- c(2, 2, 2)
  w <- diag(3)
  # Python: np.dot((a-b), np.dot(w, (a-b)))
  expect_equal(
    compute_statistic_weighted(a, b, w),
    as.numeric(t(a-b) %*% w %*% (a-b)),
    tolerance = 1e-10
  )
})

# Test: covariance_of_moments matches numpy cov

test_that("covariance_of_moments matches numpy cov", {
  set.seed(42)
  mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
  # Python: np.cov(mat, rowvar=False)
  # R's cov uses columns as variables, same as rowvar=False
  expect_equal(
    covariance_of_moments(mat),
    cov(mat),
    tolerance = 1e-10
  )
})
