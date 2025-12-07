#' @title Power Analysis for Specification Tests
#' @description Main workflow for running power analysis using modularized functions.
#' @param n Integer, sample size.
#' @param n_sim Integer, number of simulations.
#' @param n_boot Integer, number of bootstrap samples.
#' @param n_rep Integer, number of repetitions for empirical size.
#' @param mu_values Numeric vector, means to test.
#' @param c_values Numeric vector, c values for beta selection.
#' @param moment_order Integer, order of moments.
#' @param alpha Numeric, significance level.
#' @param r_dir Character, directory for R scripts (default 'R').
#' @return Data frame of power results.
#' @examples
#' run_power_analysis()
#' @family main

run_power_analysis <- function(
  n = 100,
  n_sim = 8,
  n_boot = 100,
  n_rep = 1000,
  mu_values = seq(0, 2.1, 0.3),
  c_values = c(50, 20, 10, 1e-1, 1e-3, 1e-6, 1e-9, 0),
  moment_order = 5,
  alpha = 0.05,
  r_dir = "R" # allow override for testthat
) {
  library(tidyverse)
  if (!requireNamespace("cli", quietly = TRUE)) stop("cli package required for input validation.")
  cli::cli_h2("Validating inputs for run_power_analysis()")
  if (!is.numeric(n) || length(n) != 1 || n < 1) cli::cli_abort("n must be a positive integer.")
  if (!is.numeric(n_sim) || length(n_sim) != 1 || n_sim < 1) cli::cli_abort("n_sim must be a positive integer.")
  if (!is.numeric(n_boot) || length(n_boot) != 1 || n_boot < 1) cli::cli_abort("n_boot must be a positive integer.")
  if (!is.numeric(n_rep) || length(n_rep) != 1 || n_rep < 1) cli::cli_abort("n_rep must be a positive integer.")
  if (!is.numeric(mu_values) || length(mu_values) < 1) cli::cli_abort("mu_values must be a numeric vector.")
  if (!is.numeric(c_values) || length(c_values) < 1) cli::cli_abort("c_values must be a numeric vector.")
  if (!is.numeric(moment_order) || length(moment_order) != 1 || moment_order < 1) cli::cli_abort("moment_order must be a positive integer.")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) cli::cli_abort("alpha must be a number in (0, 1).")
  if (!is.character(r_dir) || length(r_dir) != 1) cli::cli_abort("r_dir must be a character string.")
  source(file.path(r_dir, "moment_utils.R"))
  source(file.path(r_dir, "test_statistics.R"))
  source(file.path(r_dir, "beta_selection.R"))
  best_beta_results <- map(
    mu_values,
    ~set_names(vector("list", moment_order), as.character(1:moment_order))
  ) %>% set_names(as.character(mu_values))
  power_test_1 <- numeric(length(mu_values))
  power_test_1opt <- numeric(length(mu_values))
  # Parallelize over mu_values
  results <- furrr::future_map2_dfr(
    mu_values, seq_along(mu_values),
    function(mu, idx) {
      rejected_1 <- 0
      rejected_1opt <- 0
      p_values_1 <- numeric(n_rep)
      p_values_1opt <- numeric(n_rep)
      for (rep in seq_len(n_rep)) {
        observed_data <- rnorm(n, mu, 1)
        observed_psi <- map_dbl(1:moment_order, ~compute_moments(observed_data, .x))
        true_data1 <- rnorm(n, 0, 1)
        best_beta <- select_optimal_beta(n, c_values, moment_order, n_sim, n_boot, true_data1)
        best_beta_results[[as.character(mu)]][[as.character(moment_order)]] <-
          c(best_beta_results[[as.character(mu)]][[as.character(moment_order)]], best_beta)
        psi_s <- map(1:n_sim, ~map_dbl(1:moment_order, ~compute_moments(rnorm(n), .x))) %>% do.call(rbind, .)
        bar_psi_S <- colMeans(psi_s)
        K_N <- compute_kernel(psi_s, n)
        K_N_square <- K_N %*% K_N
        K_N_inv_half <- MASS::ginv(K_N_square + best_beta * diag(moment_order)) %*% K_N
        W_N_1 <- calculate_test_statistic_1(observed_psi, bar_psi_S)
        W_N_1opt <- calculate_test_statistic_1opt(observed_psi, bar_psi_S, K_N_inv_half)
        psi_r <- map(1:n_boot, ~map_dbl(1:moment_order, ~compute_moments(rnorm(n), .x))) %>% do.call(rbind, .)
        W_r_1opt <- map_dbl(1:n_boot, ~calculate_test_statistic_1opt(psi_r[.x,], bar_psi_S, K_N_inv_half))
        p_value_1opt <- (sum(W_r_1opt >= W_N_1opt) + 1) / (n_boot + 1)
        W_r_1 <- map_dbl(1:n_boot, ~calculate_test_statistic_1(psi_r[.x,], bar_psi_S))
        p_value_1 <- (sum(W_r_1 >= W_N_1) + 1) / (n_boot + 1)
        p_values_1[rep] <- p_value_1
        p_values_1opt[rep] <- p_value_1opt
        if (p_value_1 < alpha) rejected_1 <- rejected_1 + 1
        if (p_value_1opt < alpha) rejected_1opt <- rejected_1opt + 1
      }
      # Print summary statistics for diagnostics
      print(glue::glue('mu={mu}: Test1 p-value range: {range(p_values_1)}, mean={mean(p_values_1)}'))
      print(glue::glue('mu={mu}: Test1opt p-value range: {range(p_values_1opt)}, mean={mean(p_values_1opt)}'))
      tibble(
        mu = mu,
        power_test_1 = rejected_1 / n_rep,
        power_test_1opt = rejected_1opt / n_rep
      )
    }, .options = furrr::furrr_options(seed = TRUE)
  )
  print(results)
  results
}
