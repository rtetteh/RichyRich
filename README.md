# RichyRich: Moment-Based Hypothesis Testing and Simulation

## Overview

This project provides modular R scripts for conducting moment-based hypothesis tests, optimal weighting via regularization, and power simulations. All code follows the tidyverse style guide, uses snake_case naming, and is pipe-friendly. The workflow is suitable for advanced statistical simulation and robust hypothesis testing.

## Core Functionalities

### 1. Sample Moments and Test Statistics (`R/moments_and_statistics.R`)
- **calculate_moments(x, max_order):**
  Computes the mean of a numeric vector raised to each integer power up to `max_order`.
- **compute_statistic_unweighted(observed, reference):**
  Calculates the squared Euclidean distance between two vectors (unweighted test statistic).
- **compute_statistic_weighted(observed, reference, weight_mat):**
  Computes a quadratic form using a weighting matrix (e.g., covariance matrix) for Mahalanobis-style statistics.
- **covariance_of_moments(sim_moments):**
  Computes the sample covariance matrix of simulated moment vectors.

### 2. Regularization Selection (`R/regularization_selection.R`)
- **select_optimal_regularization(n, c_grid, moment_order, n_sim, n_boot, x):**
  Cross-validates over candidate regularization values to minimize the difference in p-values between train/test splits, returning the optimal regularization parameter for weighting matrices.

### 3. Power Simulation (`R/power_simulation.R`)
- **run_power_simulation(n, n_sim, n_boot, alpha, moment_order, n_rep, mu_grid, c_grid):**
  Simulates the empirical power of both unweighted and optimally weighted moment-based tests across a grid of means. Returns a tibble with power results for each test.

### 4. Visualization (`R/plot_power_curves.R`)
- **plot_power_curves(power_results, alpha):**
  Plots power curves for both unweighted and weighted tests using ggplot2, with a reference line for the significance level.

### 5. Main Workflow (`main_power_analysis.R`)
- Sources all modules, sets up parameters, runs the power simulation, and generates plots. This script demonstrates the full workflow and can be adapted for custom analyses.

## Documentation and Examples

- **Vignettes:**
  - `vignettes/moments_and_statistics.qmd`: How to use the core moment/statistic functions.
  - `vignettes/regularization_selection.qmd`: How to select optimal regularization for weighting matrices.
  - `vignettes/power_simulation.qmd`: How to run and visualize power simulations.

Each vignette provides setup, code examples, and interpretation.

## Style and Best Practices
- All code uses snake_case for functions and variables (tidyverse standard).
- Only native R pipe (`|>`) is used.
- All functions are pipe-friendly and return tibbles or vectors.
- Roxygen2 documentation is provided for all public functions.

## Getting Started
1. Source the relevant R scripts in your session or use as a package.
2. See the vignettes in the `vignettes/` folder for usage examples and workflows.
3. Run `main_power_analysis.R` for a full demonstration.

## Requirements
- R (>= 4.1.0)
- Packages: `purrr`, `tibble`, `ggplot2`, `MASS`, `dplyr` (all tidyverse)

## License
MIT
