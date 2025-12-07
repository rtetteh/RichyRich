# RichyRich Power Analysis Toolkit

This toolkit provides modular R functions for conducting power analysis of specification tests using higher-order moments and optimal weighting. The code is organized for reproducibility, parallelization, and tidyverse compatibility.

## Directory Structure

- `R/moment_utils.R`: Functions for computing moments and kernel matrices.
- `R/test_statistics.R`: Functions for calculating test statistics.
- `R/beta_selection.R`: Function to select the optimal beta parameter for weighting.
- `R/power_analysis.R`: Main workflow for running the power analysis simulation.
- `R/plotting.R`: Function to plot power analysis results.
- `R/main.R`: Example script to run the full workflow.

## Function Overview

### 1. `compute_moments(data, order)`
Computes the sample moment of a given order for a numeric vector.

### 2. `compute_kernel(matrix, sample_size)`
Calculates the kernel (covariance) matrix of simulated moments.

### 3. `calculate_test_statistic_1(psi, psi_bar)`
Computes the basic test statistic as the squared Euclidean distance between observed and simulated moments.

### 4. `calculate_test_statistic_1opt(psi, psi_bar, K_inv_half)`
Computes the optimally weighted test statistic using a kernel-based weighting matrix.

### 5. `select_optimal_beta(n, c_values, moment_order, n_sim, n_boot, data_true)`
Selects the optimal beta parameter for the weighting matrix by minimizing the difference in p-values between training and test splits.

### 6. `run_power_analysis(...)`
Runs the full power analysis simulation for a range of means (`mu_values`), returning a tibble of empirical power for each test.

### 7. `plot_power_results(results, alpha = 0.05)`
Plots the empirical power curves for both test statistics.

## How to Use

1. **Install dependencies:**
   - `tidyverse`, `furrr`, `MASS`, `glue`
2. **Source the main script:**
   ```r
   source('R/main.R')
   ```
   This will run the power analysis and plot the results.
3. **Customize parameters:**
   - Edit `R/main.R` to change sample size, number of simulations, bootstrap samples, repetitions, or the range of means.

## Example
```r
# In R console or script:
source('R/main.R')
```

## Notes
- All code follows the tidyverse style guide and uses the native R pipe (`|>`).
- Parallelization is handled via `furrr::future_map2_dfr`.
- All functions are documented with roxygen2 comments.
- For further customization, edit the relevant R scripts in the `R/` directory.
