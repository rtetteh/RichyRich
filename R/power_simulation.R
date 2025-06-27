#' Power simulation for moment-based test
#' @description Simulates power for both unweighted and optimally weighted tests across a grid of means, writing results to DuckDB for memory efficiency
#' @param n Integer, sample size
#' @param n_sim Integer, number of simulated null samples
#' @param n_boot Integer, number of bootstrap samples
#' @param alpha Numeric, significance level
#' @param moment_order Integer
#' @param n_rep Integer, number of repetitions per mean
#' @param mu_grid Numeric vector, means to test
#' @param c_grid Numeric vector, regularization candidates
#' @param verbose Logical, if TRUE print user alerts via cli (default FALSE)
#' @return Tibble with columns mu, power_unweighted, power_weighted
#' @importFrom purrr map map_dbl map_int
#' @importFrom dplyr bind_rows
#' @importFrom furrr future_map
#' @importFrom duckdb duckdb, dbConnect, dbDisconnect
#' @import progressr
run_power_simulation <- function(n, n_sim, n_boot, alpha, moment_order, n_rep, mu_grid, c_grid, verbose = FALSE) {
  if (!is.numeric(n) || length(n) != 1 || n < 1 || n != as.integer(n)) {
    stop("n must be a single positive integer.")
  }
  if (!is.numeric(n_sim) || length(n_sim) != 1 || n_sim < 1 || n_sim != as.integer(n_sim)) {
    stop("n_sim must be a single positive integer.")
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1 || n_boot < 1 || n_boot != as.integer(n_boot)) {
    stop("n_boot must be a single positive integer.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1 || is.na(alpha) || !is.finite(alpha)) {
    stop("alpha must be a single numeric value between 0 and 1 (not NA/Inf).")
  }
  if (!is.numeric(moment_order) || length(moment_order) != 1 || moment_order < 1 || moment_order != as.integer(moment_order)) {
    stop("moment_order must be a single positive integer.")
  }
  if (!is.numeric(n_rep) || length(n_rep) != 1 || n_rep < 1 || n_rep != as.integer(n_rep)) {
    stop("n_rep must be a single positive integer.")
  }
  if (!is.numeric(mu_grid) || length(mu_grid) < 1 || anyNA(mu_grid) || any(!is.finite(mu_grid))) {
    stop("mu_grid must be a non-empty numeric vector with no NA/Inf.")
  }
  if (!is.numeric(c_grid) || length(c_grid) < 1 || anyNA(c_grid) || any(!is.finite(c_grid))) {
    stop("c_grid must be a non-empty numeric vector with no NA/Inf.")
  }
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value.")
  }
  if (verbose) cli::cli_alert_info("Starting power simulation with {length(mu_grid)} mu values and {n_rep} repetitions per mu...")
  # Set up parallel plan: use one less than the number of logical cores
  n_workers <- max(1, parallel::detectCores() - 1)
  future::plan(future::multisession, workers = n_workers)
  if (verbose) cli::cli_alert_info("Parallel plan set with {n_workers} workers.")
  # Set up a temporary DuckDB database
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  DBI::dbExecute(con, "CREATE TABLE sim_results (mu DOUBLE, reject_unweighted INTEGER, reject_weighted INTEGER)")
  if (verbose) cli::cli_alert_info("DuckDB in-memory database initialized.")
  # Set up progress bar for mu_grid (outer parallel loop only)
  pb_mu <- progressr::progressor(along = mu_grid)
  all_results <- furrr::future_map(
    mu_grid,
    function(mu) {
      pb_mu(sprintf("mu = %s", mu))
      results <- purrr::map(
        seq_len(n_rep),
        function(j) {
          obs <- rnorm(n, mean = mu, sd = 1)
          obs_moments <- calculate_moments(obs, moment_order, central = TRUE)
          null_data <- rnorm(n)
          beta <- select_optimal_regularization(n, c_grid, moment_order, n_sim, n_boot, null_data)
          sim_moments <- purrr::map(seq_len(n_sim), ~calculate_moments(rnorm(n), moment_order, central = TRUE))
          sim_moments <- purrr::map(sim_moments, function(x) if (is.null(x) || !is.numeric(x)) tibble::tibble() else tibble::as_tibble_row(x))
          sim_moments <- dplyr::bind_rows(sim_moments)
          if (nrow(sim_moments) == 0) return(tibble::tibble(reject_unweighted = NA_real_, reject_weighted = NA_real_))
          sim_moments_mat <- as.matrix(sim_moments)
          mean_sim <- colMeans(sim_moments_mat)
          cov_sim <- covariance_of_moments(sim_moments_mat)
          weight_mat <- MASS::ginv(cov_sim + beta * diag(moment_order))
          stat_unweighted <- compute_statistic_unweighted(obs_moments, mean_sim)
          stat_weighted <- compute_statistic_weighted(obs_moments, mean_sim, weight_mat)
          boot_moments <- purrr::map(seq_len(n_boot), ~calculate_moments(rnorm(n), moment_order, central = TRUE))
          boot_moments <- purrr::map(boot_moments, function(x) if (is.null(x) || !is.numeric(x)) tibble::tibble() else tibble::as_tibble_row(x))
          boot_moments <- dplyr::bind_rows(boot_moments)
          if (nrow(boot_moments) == 0) return(tibble::tibble(reject_unweighted = NA_real_, reject_weighted = NA_real_))
          boot_moments_mat <- as.matrix(boot_moments)
          boot_stats_unweighted <- purrr::map_dbl(
            seq_len(n_boot),
            ~compute_statistic_unweighted(boot_moments_mat[., ], mean_sim)
          )
          boot_stats_weighted <- purrr::map_dbl(
            seq_len(n_boot),
            ~compute_statistic_weighted(boot_moments_mat[., ], mean_sim, weight_mat)
          )
          pval_unweighted <- (sum(boot_stats_unweighted >= stat_unweighted) + 1) / (n_boot + 1)
          pval_weighted <- (sum(boot_stats_weighted >= stat_weighted) + 1) / (n_boot + 1)
          tibble::tibble(
            reject_unweighted = as.integer(pval_unweighted < alpha),
            reject_weighted = as.integer(pval_weighted < alpha)
          )
        }
      )
      results <- purrr::map(results, function(x) if (is.null(x) || !is.data.frame(x)) tibble::tibble(reject_unweighted = NA_real_, reject_weighted = NA_real_) else x)
      results <- dplyr::bind_rows(results)
      # Always add mu column, even if results is empty, and ensure correct type
      results <- tibble::add_column(results, mu = as.numeric(mu), .before = 1)
      results
    },
    .options = furrr::furrr_options(seed = TRUE)
  )
  # Bind all results and write to DuckDB in the main process
  all_results <- dplyr::bind_rows(all_results)
  # Ensure 'mu' column always exists and is numeric, even if all_results is empty
  if (!"mu" %in% names(all_results)) {
    all_results <- tibble::add_column(all_results, mu = numeric(0), .before = 1)
  }
  if (nrow(all_results) > 0 && "mu" %in% names(all_results)) {
    DBI::dbAppendTable(con, "sim_results", all_results[, c("mu", "reject_unweighted", "reject_weighted")])
    sim_tbl <- duckplyr::as_duckdb_tibble(dplyr::tbl(con, "sim_results"))
    summary <- sim_tbl |
      dplyr::group_by(mu) |
      dplyr::summarise(
        power_unweighted = mean(reject_unweighted, na.rm = TRUE),
        power_weighted = mean(reject_weighted, na.rm = TRUE),
        .groups = "drop"
      ) |
      dplyr::collect()
    if (verbose) cli::cli_alert_success("Power simulation complete. Returning summary tibble.")
    return(summary)
  } else {
    if (verbose) cli::cli_alert_warning("No valid results to summarize: returning empty summary tibble.")
    return(tibble::tibble(mu = numeric(0), power_unweighted = numeric(0), power_weighted = numeric(0)))
  }
}
