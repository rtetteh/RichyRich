#' Select optimal c* and beta*
#' @param N Integer, sample size
#' @param c_values Numeric vector, candidate c values
#' @param moment_order Integer
#' @param S Integer, number of simulations
#' @param R Integer, number of bootstrap samples
#' @param true_data Numeric vector
#' @param verbose Logical, if TRUE print user alerts via cli (default FALSE)
#' @return Numeric, best beta
select_optimal_beta <- function(N, c_values, moment_order, S, R, true_data, verbose = FALSE) {
  # Input validation
  if (!is.numeric(N) || length(N) != 1 || N < 2 || N != as.integer(N)) {
    stop("N must be a single integer >= 2.")
  }
  if (!is.numeric(c_values) || length(c_values) < 1 || anyNA(c_values) || any(!is.finite(c_values))) {
    stop("c_values must be a non-empty numeric vector with no NA/Inf.")
  }
  if (!is.numeric(moment_order) || length(moment_order) != 1 || moment_order < 1 || moment_order != as.integer(moment_order)) {
    stop("moment_order must be a single positive integer.")
  }
  if (!is.numeric(S) || length(S) != 1 || S < 1 || S != as.integer(S)) {
    stop("S must be a single positive integer.")
  }
  if (!is.numeric(R) || length(R) != 1 || R < 1 || R != as.integer(R)) {
    stop("R must be a single positive integer.")
  }
  if (!is.numeric(true_data) || !is.vector(true_data) || length(true_data) < N || anyNA(true_data) || any(!is.finite(true_data))) {
    stop("true_data must be a numeric vector of length at least N, with no NA/Inf.")
  }
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value.")
  }
  if (verbose) cli::cli_alert_info("Selecting optimal beta (cross-validation over c_values of length {length(c_values)})...")
  set.seed(1234)
  true_data <- sample(true_data)
  tilde_N <- as.integer(2 * N / 3)
  train_data <- true_data[1:tilde_N]
  test_data <- true_data[(tilde_N + 1):N]
  N_test <- N - tilde_N
  all_candidates <- list()
  for (c in c_values) {
    beta_train <- c / (tilde_N ^ (1/3))
    beta_test <- c / (N_test ^ (1/3))
    psi_train <- sapply(1:moment_order, function(h) compute_moments(train_data, h))
    psi_s_train <- replicate(S, sapply(1:moment_order, function(h) compute_moments(rnorm(tilde_N), h)))
    psi_bar_s_train <- rowMeans(psi_s_train)
    psi_r_train <- replicate(R, sapply(1:moment_order, function(h) compute_moments(rnorm(tilde_N), h)))
    K_train <- compute_kernel(t(psi_s_train), tilde_N)
    K_square <- K_train %*% K_train
    K_c_inv_half_train <- MASS::ginv(K_square + beta_train * diag(moment_order)) %*% K_train
    W_train <- calculate_test_statistic_1opt(psi_train, psi_bar_s_train, K_c_inv_half_train)
    W_r_train <- sapply(1:R, function(i) calculate_test_statistic_1opt(psi_r_train[,i], psi_bar_s_train, K_c_inv_half_train))
    p_value_train <- (sum(W_r_train >= W_train) + 1) / (R + 1)
    psi_test <- sapply(1:moment_order, function(h) compute_moments(test_data, h))
    psi_s_test <- replicate(S, sapply(1:moment_order, function(h) compute_moments(rnorm(N_test), h)))
    psi_bar_s_test <- rowMeans(psi_s_test)
    psi_r_test <- replicate(R, sapply(1:moment_order, function(h) compute_moments(rnorm(N_test), h)))
    K_test <- compute_kernel(t(psi_r_test), N_test)
    K_square_test <- K_test %*% K_test
    K_c_inv_half_test <- MASS::ginv(K_square_test + beta_test * diag(moment_order)) %*% K_test
    W_test <- calculate_test_statistic_1opt(psi_test, psi_bar_s_test, K_c_inv_half_test)
    W_r_test <- sapply(1:R, function(i) calculate_test_statistic_1opt(psi_r_test[,i], psi_bar_s_test, K_c_inv_half_test))
    p_value_test <- (sum(W_r_test >= W_test) + 1) / (R + 1)
    abs_diff <- abs(p_value_train - p_value_test)
    beta <- c / (N ^ (1/3))
    all_candidates[[length(all_candidates) + 1]] <- list(c = c, beta = beta, p_train = p_value_train, p_test = p_value_test, diff = abs_diff)
  }
  if (length(all_candidates) == 0) {
    stop("No valid candidates found in c_values.")
  }
  best_candidate <- all_candidates[[which.min(sapply(all_candidates, function(x) c(x$diff, x$beta)))[1]]]
  if (verbose) cli::cli_alert_success("Optimal beta selected: beta = {best_candidate$beta}")
  best_candidate$beta
}
