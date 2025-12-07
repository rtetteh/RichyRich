#' @title Power Analysis Plotting
#' @description Plot power analysis results for specification tests.
#' @param results Data frame from run_power_analysis.
#' @param alpha Numeric, significance level.
#' @return ggplot object.
#' @examples
#' plot_power_results(results, alpha = 0.05)
#' @family plotting

plot_power_results <- function(results, alpha = 0.05) {
  if (!requireNamespace("cli", quietly = TRUE)) stop("cli package required for input validation.")
  if (!is.data.frame(results) && !tibble::is_tibble(results)) cli::cli_abort("results must be a data frame or tibble.")
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) cli::cli_abort("alpha must be a number in (0, 1).")
  library(ggplot2)
  ggplot(results, aes(x = mu)) +
    geom_line(aes(y = power_test_1, color = 'Test 1'), size = 1) +
    geom_line(aes(y = power_test_1opt, color = 'Test 1opt'), size = 1) +
    geom_hline(yintercept = alpha, linetype = 'dashed', color = 'red') +
    labs(
      title = 'Power Analysis for Deviation in Mean Values',
      x = 'Mean (μ)',
      y = 'Empirical Power',
      color = 'Test'
    ) +
    theme_minimal()
}
