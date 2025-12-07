from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np
import os
import pandas as pd

from .distribution import HybridLogNormalPareto, HybridParams
from .estimation import ParameterEstimator
from .testing import HypothesisTest


@dataclass(frozen=True)
class ExperimentConfig:
    """Configuration for an empirical size experiment.

    This object captures all tunable knobs of the simulation in a
    single, well-documented structure so that the public API is
    explicit and stable.
    """

    n_data: int = 2000
    mu_true: float = 9.0
    s2_true: float = 1.17
    alpha_true: float = 3.05
    tau_true: float = 32000.0
    tau_perc_true: float = 0.82
    M: int = 10
    S: int = 200
    R: int = 299
    alpha: float = 0.05
    orders: List[int] = None
    c_values: Optional[np.ndarray] = None
    tau_grid: Optional[np.ndarray] = None
    seed: int = 1234
    n_jobs: int = -1
    use_rich: bool = True
    verbose: bool = True

    def __post_init__(self) -> None:
        if self.orders is None:
            object.__setattr__(self, "orders", [4, 8])


@dataclass(frozen=True)
class OrderResult:
    """Summary results for a single moment order."""

    order: int
    best_tau: float
    mu_est: float
    sigma_est: float
    alpha_est: float
    empirical_size_test1: float
    empirical_size_test1opt: float
    rejections_test1: int
    rejections_test1opt: int
    M_repetitions: int
    nominal_size: float


@dataclass(frozen=True)
class EmpiricalSizeResult:
    """Full result of an empirical size experiment.

    Provides convenience helpers to work with the results as a
    pandas DataFrame or formatted text summary.
    """

    config: ExperimentConfig
    order_results: List[OrderResult]
    true_data: np.ndarray

    def to_dataframe(self) -> pd.DataFrame:
        rows = [
            {
                "Order": r.order,
                "Best_Tau": r.best_tau,
                "mu_est": r.mu_est,
                "sigma_est": r.sigma_est,
                "alpha_est": r.alpha_est,
                "Nominal_Size": r.nominal_size,
                "Test 1": r.empirical_size_test1,
                "Test 1opt": r.empirical_size_test1opt,
                "Rejections_Test1": r.rejections_test1,
                "Rejections_Test1opt": r.rejections_test1opt,
                "M_Repetitions": r.M_repetitions,
            }
            for r in self.order_results
        ]
        return pd.DataFrame(rows)

    def summary(self) -> str:
        df = self.to_dataframe()
        if df.empty:
            return "No results available."
        with pd.option_context("display.max_columns", None, "display.width", 120):
            return df.to_string(index=False, float_format="%.4f")


def run_empirical_size_experiment(config: ExperimentConfig) -> EmpiricalSizeResult:
    """Run the full empirical size experiment.

    This is the main public API. It mirrors the logic of the
    original notebook and of :class:`SimulationRunner`, but
    exposes a clean functional interface and typed result object.
    """

    true_params = HybridParams(
        mu=config.mu_true,
        sigma=float(np.sqrt(config.s2_true)),
        alpha=config.alpha_true,
        tau=config.tau_true,
        tau_perc=config.tau_perc_true,
    )
    true_data = HybridLogNormalPareto.generate_sample(
        config.n_data, true_params, seed=config.seed
    )

    tau_grid = (
        config.tau_grid
        if config.tau_grid is not None
        else np.arange(70, 91, 1)
    )
    tau_values = np.percentile(true_data, tau_grid)

    if config.c_values is None:
        c_values = np.concatenate(
            ([50, 40, 30, 20, 15, 10], np.logspace(1, -7, num=20))
        )
    else:
        c_values = config.c_values

    order_results: List[OrderResult] = []

    for order in config.orders:
        if config.verbose:
            print(f"Processing Order {order}...")

        best_tau = HypothesisTest.evaluate_tau(true_data, tau_values, order, S=config.S)
        if best_tau is None:
            if config.verbose:
                print(f"Could not find optimal tau for order {order}")
            continue

        mu_hat, sigma2_hat, alpha_hat = ParameterEstimator.estimate_parameters(
            best_tau, true_data
        )
        sigma_hat = float(np.sqrt(sigma2_hat))

        est_params = HybridParams(
            mu=mu_hat,
            sigma=sigma_hat,
            alpha=alpha_hat,
            tau=best_tau,
        )

        from joblib import Parallel, delayed
        from .statistics import AuxiliaryStatistics

        # Determine effective number of workers
        if config.n_jobs == -1:
            cpu_count = os.cpu_count() or 1
            effective_n_jobs = max(1, cpu_count - 1)
        else:
            effective_n_jobs = config.n_jobs

        if config.verbose:
            print(
                f"Running {config.M} Monte Carlo repetitions with n_jobs="
                f"{effective_n_jobs} (requested {config.n_jobs})..."
            )

        def _single(iteration: int) -> tuple[bool, bool]:
            iter_seed = config.seed + iteration * 1000

            observed_data = HybridLogNormalPareto.generate_sample(
                config.n_data, est_params, seed=iter_seed
            )
            observed_psi = AuxiliaryStatistics.compute_vector(observed_data, order)

            psi_s = [
                AuxiliaryStatistics.compute_vector(
                    HybridLogNormalPareto.generate_sample(
                        config.n_data, est_params, seed=iter_seed + s
                    ),
                    order,
                )
                for s in range(config.S)
            ]
            psi_s_matrix = np.array(psi_s)
            bar_psi = np.mean(psi_s_matrix, axis=0)

            W_N_1 = HypothesisTest.calculate_test_statistic_1(observed_psi, bar_psi)

            # best_tau is guaranteed non-None here because of the early
            # continue above when evaluate_tau returns None.
            best_beta = HypothesisTest.select_optimal_beta(
                config.n_data,
                c_values,
                order,
                config.S,
                config.R,
                true_data,
                best_tau,  # type: ignore[arg-type]
                tau_grid,
                seed=iter_seed,
            )

            K_T = HypothesisTest.compute_kernel(psi_s_matrix, config.S)
            squared_integral_operator = np.dot(K_T, K_T)
            pinv = np.linalg.pinv(
                squared_integral_operator + best_beta * np.eye(order)
            )
            optimal_weight = np.dot(pinv, K_T)

            W_N_1opt = HypothesisTest.calculate_test_statistic_1opt(
                observed_psi, bar_psi, optimal_weight
            )

            stats_1 = []
            stats_1opt = []

            for b in range(config.R):
                boot_data = HybridLogNormalPareto.generate_sample(
                    config.n_data, est_params, seed=iter_seed + config.S + b
                )
                boot_psi = AuxiliaryStatistics.compute_vector(boot_data, order)
                stats_1.append(
                    HypothesisTest.calculate_test_statistic_1(bar_psi, boot_psi)
                )
                stats_1opt.append(
                    HypothesisTest.calculate_test_statistic_1opt(
                        bar_psi, boot_psi, optimal_weight
                    )
                )

            p_val_1 = float(np.mean(np.array(stats_1) >= W_N_1))
            p_val_1opt = float(np.mean(np.array(stats_1opt) >= W_N_1opt))

            return p_val_1 < config.alpha, p_val_1opt < config.alpha
        # Optional warm-up: run a few repetitions sequentially to estimate
        # per-repetition cost and provide a rough ETA.
        warmup_reps = min(3, max(1, config.M // 10))
        warmup_time: float | None = None

        if config.verbose and warmup_reps > 0:
            from time import perf_counter

            print(f"Warm-up: running {warmup_reps} sequential repetitions for timing...")
            t0 = perf_counter()
            _warmup_results = [_single(i) for i in range(warmup_reps)]
            t1 = perf_counter()
            warmup_time = t1 - t0
            per_rep = warmup_time / warmup_reps
            est_total = per_rep * config.M / max(1, effective_n_jobs)
            print(
                f"Estimated wall-clock per order ≈ {est_total:0.1f} s "
                f"({per_rep:0.2f} s per repetition, {effective_n_jobs} workers)"
            )

        # Execute remaining Monte Carlo repetitions with optional Rich progress bar
        if config.use_rich and config.verbose:
            from rich.console import Console
            from rich.progress import (
                Progress,
                BarColumn,
                TimeElapsedColumn,
                TimeRemainingColumn,
                TextColumn,
            )

            console = Console()
            with Progress(
                TextColumn("[bold blue]Order {task.fields[order]}[/]"),
                BarColumn(),
                "[progress.percentage]{task.percentage:>3.0f}%",
                "•",
                TimeElapsedColumn(),
                "remaining",
                TimeRemainingColumn(),
                console=console,
            ) as progress:
                task_id = progress.add_task("mc", total=config.M, order=order)

                # Use Parallel to compute all results, but advance the
                # progress bar from the main process as each result is
                # consumed from the iterator.
                parallel = Parallel(n_jobs=effective_n_jobs)

                results = []

                # Start with any warm-up results we may have computed
                completed = 0
                for _ in range(warmup_reps):
                    progress.update(task_id, advance=1)
                    completed += 1

                remaining_indices = range(warmup_reps, config.M)
                iterator = parallel(delayed(_single)(i) for i in remaining_indices)

                for res in iterator:
                    results.append(res)
                    completed += 1
                    progress.update(task_id, advance=1)
        else:
            # If no Rich progress, just run all repetitions in parallel,
            # including any warm-up already done above (their cost is minor).
            results = Parallel(n_jobs=effective_n_jobs)(
                delayed(_single)(i) for i in range(config.M)
            )

        # Filter out any unexpected None entries for safety
        filtered_results = [r for r in results if r is not None]

        rejections_1 = sum(int(r[0]) for r in filtered_results)
        rejections_1opt = sum(int(r[1]) for r in filtered_results)

        empirical_size_1 = rejections_1 / config.M
        empirical_size_1opt = rejections_1opt / config.M

        order_results.append(
            OrderResult(
                order=order,
                best_tau=float(best_tau),
                mu_est=float(mu_hat),
                sigma_est=float(sigma_hat),
                alpha_est=float(alpha_hat),
                empirical_size_test1=float(empirical_size_1),
                empirical_size_test1opt=float(empirical_size_1opt),
                rejections_test1=int(rejections_1),
                rejections_test1opt=int(rejections_1opt),
                M_repetitions=config.M,
                nominal_size=config.alpha,
            )
        )

        if config.verbose:
            print(f"Order {order} Results:")
            print(f"  Test 1: {empirical_size_1:.4f}")
            print(f"  Test 1opt: {empirical_size_1opt:.4f}")

    return EmpiricalSizeResult(config=config, order_results=order_results, true_data=true_data)
