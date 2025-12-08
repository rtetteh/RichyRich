import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from src2.api import ExperimentConfig, run_empirical_power_experiment

def main():
    print("Starting Power Analysis for Multiple Deviations...")

    # Vector of deviations (delta values)
    deltas = np.array([0.01, 0.02, 0.05, 0.10, 0.15, 0.20])

    results_list = []

    for delta in deltas:
        print(f"\nRunning for delta = {delta}...")

        config = ExperimentConfig(
            n_data=2000,
            M=10,  # Reduced for faster testing; increase for production
            S=200,
            R=299,
            orders=[4],
            delta=delta,
            n_jobs=-1,
            verbose=False,  # Set to False to reduce output
        )

        result = run_empirical_power_experiment(config)
        df = result.to_dataframe()

        # Extract the results for order 4
        row = df.iloc[0]
        results_list.append({
            'Delta': delta,
            'Test_1_Power': row['Test 1'],
            'Test_1opt_Power': row['Test 1opt'],
            'Rejections_Test1': row['Rejections_Test1'],
            'Rejections_Test1opt': row['Rejections_Test1opt'],
            'Best_Tau': row['Best_Tau'],
            'mu_est': row['mu_est'],
            'sigma_est': row['sigma_est'],
            'alpha_est': row['alpha_est'],
        })

    # Create results DataFrame
    results_df = pd.DataFrame(results_list)

    # Save to CSV
    results_df.to_csv("power_analysis_results.csv", index=False)
    print("\nResults saved to power_analysis_results.csv")

    # Print table
    print("\nPower Analysis Results:")
    print(results_df.to_string(index=False, float_format="%.4f"))

    # Create plot
    plt.figure(figsize=(10, 6))
    plt.plot(results_df['Delta'], results_df['Test_1_Power'], 'o-', label='Test 1', color='blue', linewidth=2, markersize=8)
    plt.plot(results_df['Delta'], results_df['Test_1opt_Power'], 's-', label='Test 1opt', color='red', linewidth=2, markersize=8)
    plt.xlabel('Deviation Fraction (δ)', fontsize=12)
    plt.ylabel('Empirical Power', fontsize=12)
    plt.title('Empirical Power vs Deviation in Right Tail (Log-Cauchy)', fontsize=14)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('power_analysis_plot.png', dpi=300, bbox_inches='tight')
    print("Plot saved to power_analysis_plot.png")
    plt.show()

if __name__ == "__main__":
    main()