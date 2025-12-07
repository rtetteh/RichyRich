# Empirical Size Analysis of Hybrid LogNormal-Pareto Distributions

This project implements a modular, high-performance Python system for performing empirical size analysis of statistical tests on Hybrid LogNormal-Pareto distributions. It is designed to validate statistical methods used for analyzing heavy-tailed data (e.g., wealth or income distributions).

## Project Structure

The project is organized into a modular `src` package and execution scripts:

```text
RichyRich/
├── src/                        # Core source code
│   ├── distribution.py         # Data generation (Hybrid LogNormal-Pareto)
│   ├── estimation.py           # Parameter estimation (MLE, Hill, MoM)
│   ├── statistics.py           # Auxiliary statistics (Percentiles, Gini, Theil)
│   ├── testing.py              # Hypothesis testing & Beta selection logic
│   ├── simulation.py           # Parallel simulation runner
│   └── utils.py                # Utility functions (Matrix operations)
├── main.py                     # Command-line script to run the full simulation
├── demo_analysis.ipynb         # Jupyter notebook demonstrating system usage
├── requirements.txt            # Python dependencies
└── README.md                   # This file
```

## Installation

1. **Prerequisites**: Ensure you have Python 3.8+ installed.

2. **Install Dependencies**:
   It is recommended to use a virtual environment.

   ```bash
   pip install -r requirements.txt
   ```

   *Dependencies include: `numpy`, `pandas`, `scipy`, `joblib`, `matplotlib`, `seaborn`.*

## Usage

### 1. Running the Full Simulation (Command Line)

To run the complete Monte Carlo simulation as a background process or script:

```bash
python main.py
```

This will:

* Initialize the simulation with default parameters (defined in `main.py`).
* Run the analysis in parallel using all available CPU cores.
* Save the results to `empirical_size_results.csv`.
* Print progress and summary statistics to the console.

You can modify `main.py` to adjust parameters like sample size (`n_data`), repetitions (`M`), or moment orders (`orders`).

### 2. Interactive Analysis (Jupyter Notebook)

For a step-by-step walkthrough and visualization of the system components:

1. Open `demo_analysis.ipynb` in VS Code or Jupyter Lab.
2. Run the cells sequentially.

The notebook demonstrates:

* **Data Generation**: Visualizing the hybrid distribution and verifying percentile placement.
* **Estimation**: Recovering parameters from synthetic data.
* **Testing**: Running the hypothesis tests and optimizing parameters ($\tau$ and $\beta$).
* **Simulation**: Running a small-scale simulation directly in the notebook.

## Key Features

* **Hybrid Distribution**: Generates data with a LogNormal body and Pareto tail, ensuring the tail starts exactly at a specified percentile.
* **Robust Estimation**: Uses a combination of MLE, Hill Estimator, and Method of Moments to robustly estimate tail parameters.
* **Parallel Processing**: Uses `joblib` for efficient parallel execution of Monte Carlo simulations.
* **Modular Design**: Follows SOLID principles, making it easy to extend or modify individual components (e.g., adding a new test statistic) without breaking the rest of the system.
