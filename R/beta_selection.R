#' @title Optimal Beta Selection
#' @description Function to select optimal beta for weighting matrix in test statistic.
#' @param n Integer, sample size.
#' @param c_values Numeric vector, candidate c values.
#' @param moment_order Integer, order of moments.
#' @param n_sim Integer, number of simulations.
#' @param n_boot Integer, number of bootstrap samples.
#' @param data_true Numeric vector, true data for kernel estimation.
#' @return Numeric, best beta value.
#' @examples
#' select_optimal_beta(100, c(1,0.1), 5, 8, 100, rnorm(100))
#' @family beta_selection

select_optimal_beta <- function(n, c_values, moment_order, n_sim, n_boot, data_true) {
  if (!requireNamespace("cli", quietly = TRUE)) stop("cli package required for input validation.")
  if (!is.numeric(n) || length(n) != 1 || n < 1) cli::cli_abort("n must be a positive integer.")
  if (!is.numeric(c_values) || length(c_values) < 1) cli::cli_abort("c_values must be a numeric vector.")
  if (!is.numeric(moment_order) || length(moment_order) != 1 || moment_order < 1) cli::cli_abort("moment_order must be a positive integer.")
  if (!is.numeric(n_sim) || length(n_sim) != 1 || n_sim < 1) cli::cli_abort("n_sim must be a positive integer.")
  if (!is.numeric(n_boot) || length(n_boot) != 1 || n_boot < 1) cli::cli_abort("n_boot must be a positive integer.")
  if (!is.numeric(data_true)) cli::cli_abort("data_true must be a numeric vector.")
  library(purrr)
  set.seed(1234)
  idx <- sample(seq_along(data_true))
  tilde_n <- floor(2 * n / 3)
  train_data <- data_true[idx[1:tilde_n]]
  test_data <- data_true[idx[(tilde_n + 1):n]]
  n_test <- n - tilde_n

  candidates <- map(c_values, function(c) {
    beta_train <- c / (tilde_n)^(1/3)
    beta_test <- c / (n_test)^(1/3)
    psi_train <- map_dbl(1:moment_order, ~compute_moments(train_data, .x))
    psi_s_train <- map(1:n_sim, ~map_dbl(1:moment_order, ~compute_moments(rnorm(tilde_n), .x))) %>% do.call(rbind, .)
    psi_bar_s_train <- colMeans(psi_s_train)
    psi_r_train <- map(1:n_boot, ~map_dbl(1:moment_order, ~compute_moments(rnorm(tilde_n), .x))) %>% do.call(rbind, .)
    K_train <- compute_kernel(psi_s_train, tilde_n)
    K_square <- K_train %*% K_train
    K_c_inv_half_train <- MASS::ginv(K_square + beta_train * diag(moment_order)) %*% K_train
    W_train <- calculate_test_statistic_1opt(psi_train, psi_bar_s_train, K_c_inv_half_train)
    W_r_train <- map_dbl(1:n_boot, ~calculate_test_statistic_1opt(psi_r_train[.x,], psi_bar_s_train, K_c_inv_half_train))
    p_value_train <- (sum(W_r_train >= W_train) + 1) / (n_boot + 1)
    psi_test <- map_dbl(1:moment_order, ~compute_moments(test_data, .x))
    psi_s_test <- map(1:n_sim, ~map_dbl(1:moment_order, ~compute_moments(rnorm(n_test), .x))) %>% do.call(rbind, .)
    psi_bar_s_test <- colMeans(psi_s_test)
    psi_r_test <- map(1:n_boot, ~map_dbl(1:moment_order, ~compute_moments(rnorm(n_test), .x))) %>% do.call(rbind, .)
    K_test <- compute_kernel(psi_r_test, n_test)
    K_square_test <- K_test %*% K_test
    K_c_inv_half_test <- MASS::ginv(K_square_test + beta_test * diag(moment_order)) %*% K_test
    W_test <- calculate_test_statistic_1opt(psi_test, psi_bar_s_test, K_c_inv_half_test)
    W_r_test <- map_dbl(1:n_boot, ~calculate_test_statistic_1opt(psi_r_test[.x,], psi_bar_s_test, K_c_inv_half_test))
    p_value_test <- (sum(W_r_test >= W_test) + 1) / (n_boot + 1)
    abs_diff <- abs(p_value_train - p_value_test)
    beta <- c / (n)^(1/3)
    list(c = c, beta = beta, p_train = p_value_train, p_test = p_value_test, abs_diff = abs_diff)
  })
  best <- candidates[[which.min(map_dbl(candidates, "abs_diff"))]]
  best$beta
}
