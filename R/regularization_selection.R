#' Select optimal regularization parameter for weighting matrix
#' @description Cross-validates over candidate values to minimize difference in p-values between train/test splits
#' @param n Integer, sample size
#' @param c_grid Numeric vector, candidate regularization multipliers
#' @param moment_order Integer
#' @param n_sim Integer, number of simulated null samples
#' @param n_boot Integer, number of bootstrap samples
#' @param x Numeric vector, data from null distribution
#' @param verbose Logical, if TRUE print user alerts via cli (default FALSE)
#' @return Numeric, optimal regularization parameter
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
select_optimal_regularization <- function(n, c_grid, moment_order, n_sim, n_boot, x, verbose = FALSE) {
  if (!is.numeric(n) || length(n) != 1 || n < 1 || n != as.integer(n)) {
    stop("n must be a single positive integer.")
  }
  if (!is.numeric(c_grid) || length(c_grid) < 1) {
    stop("c_grid must be a non-empty numeric vector.")
  }
  if (!is.numeric(moment_order) || length(moment_order) != 1 || moment_order < 1 || moment_order != as.integer(moment_order)) {
    stop("moment_order must be a single positive integer.")
  }
  if (!is.numeric(n_sim) || length(n_sim) != 1 || n_sim < 1 || n_sim != as.integer(n_sim)) {
    stop("n_sim must be a single positive integer.")
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1 || n_boot < 1 || n_boot != as.integer(n_boot)) {
    stop("n_boot must be a single positive integer.")
  }
  if (!is.numeric(x) || !is.vector(x) || length(x) < n) {
    stop("x must be a numeric vector of length at least n.")
  }
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value.")
  }
  if (verbose) cli::cli_alert_info("Selecting optimal regularization parameter (cross-validation over c_grid of length {length(c_grid)})...")
  set.seed(1234)
  idx <- sample(seq_along(x))
  n_train <- floor(2 * n / 3)
  train <- x[idx[1:n_train]]
  test <- x[idx[(n_train + 1):n]]
  n_test <- n - n_train
  results <- purrr::map(c_grid, function(c_val) {
    beta_train <- c_val / (n_train ^ (1/3))
    beta_test <- c_val / (n_test ^ (1/3))
    # Train moments
    psi_train <- calculate_moments(train, moment_order)
    sim_train <- replicate(n_sim, calculate_moments(rnorm(n_train), moment_order))
    sim_train <- t(sim_train)
    mean_sim_train <- colMeans(sim_train)
    boot_train <- replicate(n_boot, calculate_moments(rnorm(n_train), moment_order))
    boot_train <- t(boot_train)
    cov_train <- covariance_of_moments(sim_train)
    weight_train <- MASS::ginv(cov_train + beta_train * diag(moment_order))
    stat_train <- compute_statistic_weighted(psi_train, mean_sim_train, weight_train)
    stat_boot_train <- purrr::map_dbl(
      seq_len(n_boot),
      ~compute_statistic_weighted(boot_train[., ], mean_sim_train, weight_train)
    )
    p_train <- (sum(stat_boot_train >= stat_train) + 1) / (n_boot + 1)
    # Test moments
    psi_test <- calculate_moments(test, moment_order)
    sim_test <- replicate(n_sim, calculate_moments(rnorm(n_test), moment_order))
    sim_test <- t(sim_test)
    mean_sim_test <- colMeans(sim_test)
    boot_test <- replicate(n_boot, calculate_moments(rnorm(n_test), moment_order))
    boot_test <- t(boot_test)
    cov_test <- covariance_of_moments(sim_test)
    weight_test <- MASS::ginv(cov_test + beta_test * diag(moment_order))
    stat_test <- compute_statistic_weighted(psi_test, mean_sim_test, weight_test)
    stat_boot_test <- purrr::map_dbl(
      seq_len(n_boot),
      ~compute_statistic_weighted(boot_test[., ], mean_sim_test, weight_test)
    )
    p_test <- (sum(stat_boot_test >= stat_test) + 1) / (n_boot + 1)
    tibble::tibble(c = c_val, beta = c_val / (n ^ (1/3)), p_train, p_test, diff = abs(p_train - p_test))
  }) |> dplyr::bind_rows()
  best <- results[which.min(results$diff), ]
  if (verbose) cli::cli_alert_success("Optimal regularization parameter selected: beta = {best$beta}")
  best$beta
}
