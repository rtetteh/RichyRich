import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from src2.api import ExperimentConfig, run_empirical_power_experiment

def main():
    # ===========================================
    # CONFIGURABLE PARAMETERS
    # ===========================================

    # Experiment parameters
    M = 10          # Monte Carlo repetitions
    S = 200         # Simulated auxiliary statistics
    R = 299         # Bootstrap replications
    orders = [4]    # Moment orders to test
    alpha = 0.05    # Nominal significance level

    # Analysis parameters
    deltas = np.array([0.01, 0.02, 0.05, 0.10, 0.15, 0.20])  # Deviation fractions
    use_mock_data = True  # Set to False to run actual experiments (slower)

    # ===========================================

    print("Starting Power Analysis for Multiple Deviations...")
    print(f"Parameters: M={M}, S={S}, R={R}, orders={orders}, alpha={alpha}")
    print(f"Using {'mock data' if use_mock_data else 'actual experiments'}")

    if use_mock_data:
        # Mock results for demonstration (based on previous runs)
        results_list = [
            {'Delta': 0.01, 'Test_1_Power': 0.4, 'Test_1opt_Power': 0.1, 'Rejections_Test1': 4, 'Rejections_Test1opt': 1, 'Best_Tau': 32524.5, 'mu_est': 9.0247, 'sigma_est': 1.0615, 'alpha_est': 3.8574},
            {'Delta': 0.02, 'Test_1_Power': 0.6, 'Test_1opt_Power': 0.1, 'Rejections_Test1': 6, 'Rejections_Test1opt': 1, 'Best_Tau': 33024.1, 'mu_est': 9.0645, 'sigma_est': 1.0856, 'alpha_est': 3.8203},
            {'Delta': 0.05, 'Test_1_Power': 0.9, 'Test_1opt_Power': 0.2, 'Rejections_Test1': 9, 'Rejections_Test1opt': 2, 'Best_Tau': 33384.8, 'mu_est': 9.1078, 'sigma_est': 1.1110, 'alpha_est': 3.6636},
            {'Delta': 0.10, 'Test_1_Power': 0.8, 'Test_1opt_Power': 0.2, 'Rejections_Test1': 8, 'Rejections_Test1opt': 2, 'Best_Tau': 32941.6, 'mu_est': 9.1191, 'sigma_est': 1.1171, 'alpha_est': 3.3341},
            {'Delta': 0.15, 'Test_1_Power': 0.7, 'Test_1opt_Power': 0.2, 'Rejections_Test1': 7, 'Rejections_Test1opt': 2, 'Best_Tau': 32212.2, 'mu_est': 9.1420, 'sigma_est': 1.1300, 'alpha_est': 3.1007},
            {'Delta': 0.20, 'Test_1_Power': 0.9, 'Test_1opt_Power': 0.2, 'Rejections_Test1': 9, 'Rejections_Test1opt': 2, 'Best_Tau': 34312.2, 'mu_est': 9.1687, 'sigma_est': 1.1444, 'alpha_est': 3.1392},
        ]
    else:
        # Run actual experiments
        results_list = []

        for delta in deltas:
            print(f"\nRunning for delta = {delta}...")

            config = ExperimentConfig(
                n_data=2000,
                M=M,
                S=S,
                R=R,
                orders=orders,
                delta=delta,
                alpha=alpha,
                n_jobs=-1,
                verbose=False,
            )

            result = run_empirical_power_experiment(config)
            df = result.to_dataframe()

            # Extract the results for the first order
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

    # Create formatted table similar to empirical size results
    print("\n" + "="*120)
    print("EMPIRICAL POWER ANALYSIS RESULTS")
    print("="*120)
    print(f"Sample Size: 2000")
    print(f"Monte Carlo Repetitions: {M}")
    print(f"Simulated Statistics (S): {S}")
    print(f"Bootstrap Replications (R): {R}")
    print(f"Moment Orders: {orders}")
    print(f"Nominal Significance Level: {alpha}")
    print(f"Delta Values Tested: {', '.join([f'{d:.2f}' for d in deltas])}")
    print(f"Data Source: {'Mock Data (Demo)' if use_mock_data else 'Actual Experiments'}")
    print("="*120)

    # Print header
    header = f"{'Delta':<8} {'Test 1':<8} {'Test 1opt':<10} {'Rejections Test1':<16} {'Rejections Test1opt':<18} {'Best Tau':<10} {'mu_est':<8} {'sigma_est':<10} {'alpha_est':<10}"
    print(header)
    print("-" * 120)

    # Print each row
    for result in results_list:
        row = (f"{result['Delta']:<8.3f} "
               f"{result['Test_1_Power']:<8.4f} "
               f"{result['Test_1opt_Power']:<10.4f} "
               f"{result['Rejections_Test1']:>2d}/{M:<13} "
               f"{result['Rejections_Test1opt']:>2d}/{M:<15} "
               f"{result['Best_Tau']:<10.1f} "
               f"{result['mu_est']:<8.4f} "
               f"{result['sigma_est']:<10.4f} "
               f"{result['alpha_est']:<10.4f}")
        print(row)

    print("="*120)

    # Create plot
    plt.figure(figsize=(12, 8))
    plt.plot(results_df['Delta'], results_df['Test_1_Power'], 'o-', label='Test 1', color='blue', linewidth=2, markersize=8)
    plt.plot(results_df['Delta'], results_df['Test_1opt_Power'], 's-', label='Test 1opt', color='red', linewidth=2, markersize=8)
    plt.axhline(y=alpha, color='black', linestyle='--', alpha=0.7, label=f'Nominal Size ({alpha})')
    plt.xlabel('Deviation Fraction (δ)', fontsize=12)
    plt.ylabel('Empirical Power', fontsize=12)
    plt.title('Empirical Power vs Deviation in Right Tail (Log-Cauchy)', fontsize=14, fontweight='bold')
    plt.legend(loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('power_analysis_plot.png', dpi=300, bbox_inches='tight')
    print("Plot saved to power_analysis_plot.png")
    plt.show()

if __name__ == "__main__":
    main()