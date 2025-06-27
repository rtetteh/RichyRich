#' Calculate sample moments up to a specified order
#' @description Computes either raw or central moments up to `max_order`.
#' @param x Numeric vector
#' @param max_order Integer, maximum moment order
#' @param central Logical, if TRUE (default), compute central moments; if FALSE, compute raw moments
#' @return Named numeric vector of moments
#' @examples
#' calculate_moments(1:5, 3, central = TRUE)
#' calculate_moments(1:5, 3, central = FALSE)
calculate_moments <- function(x, max_order, central = TRUE) {
  if (!is.numeric(x) || !is.vector(x)) {
    stop("x must be a numeric vector.")
  }
  if (!is.numeric(max_order) || length(max_order) != 1 || max_order < 1 || max_order != as.integer(max_order)) {
    stop("max_order must be a single positive integer.")
  }
  if (!is.logical(central) || length(central) != 1) {
    stop("central must be a single logical value.")
  }
  names_out <- paste0("moment_", seq_len(max_order))
  if (central) {
    mu <- mean(x)
    out <- purrr::map_dbl(seq_len(max_order), function(k) {
      if (k == 1) {
        mu
      } else {
        mean((x - mu)^k)
      }
    })
  } else {
    out <- purrr::map_dbl(seq_len(max_order), ~mean(x ^ .x))
  }
  names(out) <- names_out
  out
}

#' Mahalanobis-style test statistic (unweighted)
#' @description Squared Euclidean distance between two vectors
#' @param observed Numeric vector
#' @param reference Numeric vector
#' @return Numeric, test statistic
compute_statistic_unweighted <- function(observed, reference) {
  if (!is.numeric(observed) || !is.numeric(reference)) {
    stop("Both observed and reference must be numeric vectors.")
  }
  if (length(observed) != length(reference)) {
    stop("observed and reference must have the same length.")
  }
  sum((observed - reference) ^ 2)
}

#' Mahalanobis-style test statistic (weighted)
#' @description Quadratic form using a weighting matrix
#' @param observed Numeric vector
#' @param reference Numeric vector
#' @param weight_mat Numeric matrix
#' @return Numeric, test statistic
compute_statistic_weighted <- function(observed, reference, weight_mat) {
  if (!is.numeric(observed) || !is.numeric(reference)) {
    stop("Both observed and reference must be numeric vectors.")
  }
  if (length(observed) != length(reference)) {
    stop("observed and reference must have the same length.")
  }
  if (!is.matrix(weight_mat) || !is.numeric(weight_mat)) {
    stop("weight_mat must be a numeric matrix.")
  }
  if (nrow(weight_mat) != length(observed) || ncol(weight_mat) != length(observed)) {
    stop("weight_mat must be a square matrix with dimension equal to the length of observed/reference.")
  }
  diff <- observed - reference
  as.numeric(t(diff) %*% weight_mat %*% diff)
}

#' Covariance matrix for simulated moments
#' @description Computes the sample covariance matrix of simulated moment vectors
#' @param sim_moments Matrix (simulations x moments)
#' @return Covariance matrix
covariance_of_moments <- function(sim_moments) {
  if (!is.matrix(sim_moments) || !is.numeric(sim_moments)) {
    stop("sim_moments must be a numeric matrix.")
  }
  stats::cov(sim_moments)
}
