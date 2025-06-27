#' Plot power curves for moment-based tests
#' @param power_results Tibble with columns mu, power_unweighted, power_weighted
#' @param alpha Numeric, significance level
#' @param verbose Logical, if TRUE print user alerts via cli (default FALSE)
#' @return ggplot object
#' @import ggplot2
plot_power_curves <- function(power_results, alpha, verbose = FALSE) {
  if (!is.data.frame(power_results) ||
      !all(c("mu", "power_unweighted", "power_weighted") %in% colnames(power_results))) {
    stop("power_results must be a data frame or tibble with columns: mu, power_unweighted, power_weighted.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single numeric value between 0 and 1.")
  }
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value.")
  }
  if (verbose) cli::cli_alert_info("Plotting power curves for {nrow(power_results)} mu values...")
  p <- ggplot2::ggplot(power_results, ggplot2::aes(x = mu)) +
    ggplot2::geom_line(ggplot2::aes(y = power_unweighted, color = 'Unweighted'), linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(y = power_unweighted, color = 'Unweighted')) +
    ggplot2::geom_line(ggplot2::aes(y = power_weighted, color = 'Weighted'), linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(y = power_weighted, color = 'Weighted')) +
    ggplot2::geom_hline(yintercept = alpha, linetype = 'dashed', color = 'red') +
    ggplot2::labs(title = 'Power Analysis for Deviation in Mean',
         x = 'Mean (μ)',
         y = 'Empirical Power',
         color = 'Test') +
    ggplot2::theme_minimal()
  if (verbose) cli::cli_alert_success("Power curves plot created.")
  p
}
