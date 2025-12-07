#' @title Test Statistic Calculations
#' @description Functions for calculating test statistics for specification tests.
#' @param psi Numeric vector of observed moments.
#' @param psi_bar Numeric vector of mean simulated moments.
#' @param K_inv_half Optional, matrix for optimal weighting.
#' @return Numeric value (test statistic).
#' @examples
#' calculate_test_statistic_1(1:5, 1:5)
#' calculate_test_statistic_1opt(1:5, 1:5, diag(5))
#' @family test_statistics

calculate_test_statistic_1 <- function(psi, psi_bar) {
  if (!requireNamespace("cli", quietly = TRUE)) stop("cli package required for input validation.")
  if (!is.numeric(psi)) cli::cli_abort("psi must be numeric.")
  if (!is.numeric(psi_bar)) cli::cli_abort("psi_bar must be numeric.")
  if (length(psi) != length(psi_bar)) cli::cli_abort("psi and psi_bar must have the same length.")
  diff <- psi - psi_bar
  sum(diff^2)
}

calculate_test_statistic_1opt <- function(psi, psi_bar, K_inv_half) {
  if (!requireNamespace("cli", quietly = TRUE)) stop("cli package required for input validation.")
  if (!is.numeric(psi)) cli::cli_abort("psi must be numeric.")
  if (!is.numeric(psi_bar)) cli::cli_abort("psi_bar must be numeric.")
  if (!is.matrix(K_inv_half)) cli::cli_abort("K_inv_half must be a matrix.")
  if (length(psi) != length(psi_bar)) cli::cli_abort("psi and psi_bar must have the same length.")
  if (ncol(K_inv_half) != length(psi)) cli::cli_abort("K_inv_half must have number of columns equal to length of psi.")
  diff <- psi - psi_bar
  transformed <- K_inv_half %*% diff
  as.numeric(crossprod(transformed))
}
