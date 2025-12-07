library(testthat)
library(tidyverse)
source('../R/plotting.R')

# Test plot_power_results returns a ggplot object
test_that('plot_power_results returns ggplot', {
  df <- tibble(mu = c(0, 1), power_test_1 = c(0.1, 0.9), power_test_1opt = c(0.2, 0.8))
  plt <- plot_power_results(df, alpha = 0.05)
  expect_true(inherits(plt, 'ggplot'))
})
