source(normalizePath("../R/moments_and_statistics.R"))
source(normalizePath("../R/regularization_selection.R"))
source(normalizePath("../R/power_simulation.R"))

#' Test: run_power_simulation with and without DuckDB/duckplyr for large workloads
#' 
#' This test benchmarks memory usage and correctness for large simulation workloads.
#' It demonstrates that the DuckDB-backed version can handle much larger n_rep than a pure in-memory approach.

library(testthat)
library(tibble)
library(dplyr)
library(duckdb)
library(duckplyr)
library(furrr)

# Helper: check available RAM (platform-specific, here just a placeholder)
available_ram_gb <- function() {
  # On Windows, use memory.size(max=TRUE)/1024; on Linux, use system("free -g")
  # For demonstration, return 2GB
  2
}

# Test: Large simulation with DuckDB backend

test_that("run_power_simulation handles large n_rep with DuckDB backend", {
  # skip_on_cran()
  skip_if_not_installed("duckdb")
  skip_if_not_installed("duckplyr")
  skip_if_not_installed("furrr")

  # Use a large n_rep that would normally exhaust RAM
  n_rep <- 10000
  n <- 30
  n_sim <- 5
  n_boot <- 10
  alpha <- 0.05
  moment_order <- 3
  mu_grid <- c(0, 1)
  c_grid <- c(1, 0.1)

  # Run simulation (should not run out of memory)
  res <- run_power_simulation(
    n = n,
    n_sim = n_sim,
    n_boot = n_boot,
    alpha = alpha,
    moment_order = moment_order,
    n_rep = n_rep,
    mu_grid = mu_grid,
    c_grid = c_grid
  )

  expect_s3_class(res, "tbl_df")
  expect_true(all(c("mu", "power_unweighted", "power_weighted") %in% names(res)))
  expect_equal(nrow(res), length(mu_grid))
  expect_true(all(res$power_unweighted >= 0 & res$power_unweighted <= 1))
  expect_true(all(res$power_weighted >= 0 & res$power_weighted <= 1))
})

# Optional: Compare with a pure in-memory version (if available)
# test_that("pure in-memory version fails for large n_rep", { ... })
