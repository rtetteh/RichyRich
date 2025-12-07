#' @title Moment Calculation Utilities
#' @description Functions for computing moments and kernels for statistical tests.
#' @param data Numeric vector.
#' @param order Integer, order of the moment.
#' @return Numeric value or matrix.
#' @examples
#' compute_moments(rnorm(100), 3)
#' compute_kernel(matrix(rnorm(100*5), ncol=5), 100)
#' @family moment_utils

compute_moments <- function(data, order) {
  if (!requireNamespace("cli", quietly = TRUE)) stop("cli package required for input validation.")
  if (!is.numeric(data)) cli::cli_abort("data must be numeric.")
  if (!is.numeric(order) || length(order) != 1 || order < 1) cli::cli_abort("order must be a positive integer.")
  mean(data ^ order)
}

compute_kernel <- function(matrix, sample_size) {
  if (!requireNamespace("cli", quietly = TRUE)) stop("cli package required for input validation.")
  if (!is.matrix(matrix)) cli::cli_abort("matrix must be a matrix.")
  if (!is.numeric(sample_size) || length(sample_size) != 1 || sample_size < 1) cli::cli_abort("sample_size must be a positive integer.")
  moment_order <- ncol(matrix)
  kernel <- matrix(0, nrow = moment_order, ncol = moment_order)
  means <- colMeans(matrix)
  for (h1 in seq_len(moment_order)) {
    for (h2 in seq_len(moment_order)) {
      diff1 <- matrix[, h1] - means[h1]
      diff2 <- matrix[, h2] - means[h2]
      kernel[h1, h2] <- sum(diff1 * diff2) / sample_size
    }
  }
  kernel
}
