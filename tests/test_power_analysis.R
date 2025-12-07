library(testthat)
library(tidyverse)
source("../R/moment_utils.R")
source("../R/test_statistics.R")
source("../R/beta_selection.R")
source('../R/power_analysis.R')

# Test run_power_analysis returns a tibble with expected columns and values in [0,1]
test_that('run_power_analysis returns valid results', {
  res <- run_power_analysis(n = 10, n_sim = 2, n_boot = 2, n_rep = 2, mu_values = c(0, 1), c_values = c(1, 0.1), moment_order = 2, alpha = 0.05, r_dir = '../R')
  expect_true(is_tibble(res))
  expect_true(all(c('mu', 'power_test_1', 'power_test_1opt') %in% names(res)))
  expect_true(all(res$power_test_1 >= 0 & res$power_test_1 <= 1))
  expect_true(all(res$power_test_1opt >= 0 & res$power_test_1opt <= 1))
})
