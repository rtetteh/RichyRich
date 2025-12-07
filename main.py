from src.api import ExperimentConfig, run_empirical_size_experiment

def main():
    print("Starting Empirical Size Analysis...")

    config = ExperimentConfig(
        n_data=2000,
        M=10,
        S=200,
        R=299,
        orders=[4, 8],
        n_jobs=-1,
        verbose=True,
    )

    result = run_empirical_size_experiment(config)
    df = result.to_dataframe()

    print("\nFinal Results:")
    print(result.summary())

    df.to_csv("empirical_size_results.csv", index=False)
    print("Results saved to empirical_size_results.csv")

if __name__ == "__main__":
    main()
