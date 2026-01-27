import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from src2.api import ExperimentConfig, run_empirical_power_experiment, run_empirical_size_experiment

def main():
    # ===========================================
    # CONFIGURABLE PARAMETERS
    # ===========================================

    # Experiment parameters
    M = 1000          # Monte Carlo repetitions
    S = 200         # Simulated auxiliary statistics
    R = 299         # Bootstrap replications
    orders = [4,50]    #orders to test
    alpha = 0.05    # Nominal significance level

    # Analysis parameters
    # deltas to test for contamination fractions(same as mispecification fractions in your code)
    deltas = np.array([0,0.01,0.025,0.05,0.075,0.1])  # Contamination fractions: 2.5%, 5%, 7.5%, 10%
    taumis_value = 0.1  # Fixed tau deviation: tau_d = tau * (1 + taumis_value) = tau * 1.1
    # ===========================================

    print("Starting Size and Power Analysis...")
    print(f"Parameters: M={M}, S={S}, R={R}, orders={orders}, alpha={alpha}")
    print("Running actual experiments...")

    # Run actual experiments
    results_list = []
    total_start_time = time.perf_counter()

    for order in orders:
        print(f"\n=== Testing with {order} auxiliary statistics ===")
        
        for delta in deltas:
            delta_start_time = time.perf_counter()
            
            if delta == 0.0:
                print(f"Running SIZE experiment (alpha_alt = {delta}, order = {order})...")
                config = ExperimentConfig(
                    n_data=2000,
                    M=M,
                    S=S,
                    R=R,
                    orders=[order],  # Single order
                    delta=delta,
                    alpha=alpha,
                    n_jobs=-1,
                    verbose=False,
                )
                result = run_empirical_size_experiment(config)
            else:
                print(f"Running POWER experiment (contamination = {delta}, order = {order})...")
                config = ExperimentConfig(
                    n_data=2000,
                    M=M,
                    S=S,
                    R=R,
                    orders=[order],  # Single order
                    delta=delta,
                    taumis=taumis_value,  # Fixed tau deviation
                    alpha=alpha,
                    n_jobs=-1,
                    verbose=False,
                )
                result = run_empirical_power_experiment(config)
            
            df = result.to_dataframe()

            # Extract the results for this order
            row = df.iloc[0]
            results_list.append({
                'Order': order,
                'Delta': delta,
                'Test_1_Size_Power': row['Test 1'],
                'Test_1opt_Size_Power': row['Test 1opt'],
                'Rejections_Test1': row['Rejections_Test1'],
                'Rejections_Test1opt': row['Rejections_Test1opt'],
                'Type': 'Size' if delta == 0.0 else 'Power'
            })

            delta_end_time = time.perf_counter()
            delta_duration = delta_end_time - delta_start_time
            print(f"✓ Contamination {delta:.2f} (order {order}) completed in {delta_duration:.2f} seconds")

    total_end_time = time.perf_counter()
    total_duration = total_end_time - total_start_time
    print(f"\n🎯 Total execution time: {total_duration:.2f} seconds")
    print(f"📊 Average time per delta: {total_duration/len(deltas):.2f} seconds")

    # Create results DataFrame
    results_df = pd.DataFrame(results_list)

    # Save to CSV
    results_df.to_csv("size_power_analysis_results.csv", index=False)
    print("\nResults saved to size_power_analysis_results.csv")

    # Print results to console - separate tables for each order
    for order in orders:
        order_results = results_df[results_df['Order'] == order]
        
        print("\n" + "="*100)
        print(f"EMPIRICAL SIZE AND POWER ANALYSIS RESULTS - {order} Auxiliary Statistics")
        print("="*100)
        print(f"Sample Size: 2000 | Monte Carlo Repetitions: {M} | Alpha: {alpha} | Statistics: {order}")
        print(f"Contamination Fractions: {', '.join([f'{d:.3f}' for d in deltas])} - tau deviation fixed at {taumis_value:.1%} (τ_d = τ × {1+taumis_value:.2f})")
        print("="*100)

        # Print simplified table
        print(f"{'Contamination':<12} {'Type':<6} {'Test 1':<10} {'Test 1opt':<10} {'Rejections Test1':<16} {'Rejections Test1opt':<18}")
        print("-" * 100)

        for _, result in order_results.iterrows():
            print(f"{result['Delta']:<12.2f} {result['Type']:<6} {result['Test_1_Size_Power']:<10.4f} {result['Test_1opt_Size_Power']:<10.4f} {int(result['Rejections_Test1']):>3d}/{M:<12} {int(result['Rejections_Test1opt']):>3d}/{M:<14}")

        print("="*100)

    # Print full results DataFrame grouped by order
    print("\nDetailed Results by Order:")
    for order in orders:
        order_results = results_df[results_df['Order'] == order]
        print(f"\n--- Order {order} ---")
        print(order_results.to_string(index=False, float_format="%.4f"))

    # Create plot - Power curves for all orders (delta > 0)
    plt.figure(figsize=(14, 10))
    
    colors = ['blue', 'red', 'green', 'orange', 'purple']
    markers = ['o-', 's-', '^-', 'D-', '*-']
    
    for i, order in enumerate(orders):
        order_results = results_df[results_df['Order'] == order]
        power_df = order_results[order_results['Delta'] > 0].copy()
        size_df = order_results[order_results['Delta'] == 0].copy()
        
        color = colors[i % len(colors)]
        marker = markers[i % len(markers)]
        
        plt.plot(power_df['Delta'], power_df['Test_1_Size_Power'], marker, 
                label=f'Test 1 Power (Order {order})', color=color, linewidth=2, markersize=8)
        plt.plot(power_df['Delta'], power_df['Test_1opt_Size_Power'], marker.replace('-', '--'), 
                label=f'Test 1opt Power (Order {order})', color=color, linewidth=2, markersize=8, alpha=0.7)
        
        # Add horizontal line for empirical size
        if not size_df.empty:
            size_test1 = size_df['Test_1_Size_Power'].iloc[0]
            size_test1opt = size_df['Test_1opt_Size_Power'].iloc[0]
            plt.axhline(y=size_test1, color=color, linestyle=':', alpha=0.5, 
                       label=f'Test 1 Size (Order {order}: {size_test1:.3f})')
    
    plt.axhline(y=alpha, color='black', linestyle='--', alpha=0.7, label=f'Nominal Size ({alpha})')
    plt.axvline(x=0.0, color='red', linestyle='--', alpha=0.7, label=f'No Contamination (δ = 0)')
    plt.xlabel('Contamination Fraction (δ)', fontsize=12)
    plt.ylabel('Empirical Size/Power', fontsize=12)
    plt.title(f'Empirical Size and Power vs Contamination Fraction (τ_d = τ × {1+taumis_value:.2f})', fontsize=14, fontweight='bold')
    plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('size_power_analysis_plot.png', dpi=300, bbox_inches='tight')
    print("Plot saved to size_power_analysis_plot.png")
    plt.show()

    final_end_time = time.perf_counter()
    final_duration = final_end_time - total_start_time
    print(f"\n🏁 Analysis complete! Total runtime: {final_duration:.2f} seconds")

if __name__ == "__main__":
    main()