#!/usr/bin/env python3
"""Demo: fork-join queue approximation vs simulation."""

import pathlib
import numpy as np

from forkjoin import (
    mean_response_time,
    mean_response_time_lh,
    mean_response_time_lh_enhanced,
    upper_bound_independent,
    lower_bound_bottleneck,
    nelson_tantawi,
    simulate,
)
from forkjoin.plotting import plot_vs_load, plot_vs_heterogeneity

OUT_DIR = pathlib.Path(__file__).parent


def print_comparison_table():
    """Print approximation vs simulation for several parameter sets."""
    scenarios = [
        (1.0, 1.0, 0.3),
        (1.0, 1.0, 0.6),
        (1.0, 1.0, 0.9),
        (1.0, 1.5, 0.3),
        (1.0, 1.5, 0.6),
        (1.0, 1.5, 0.9),
        (1.0, 2.0, 0.3),
        (1.0, 2.0, 0.6),
        (1.0, 2.0, 0.9),
        (1.0, 3.0, 0.6),
        (1.0, 5.0, 0.6),
    ]

    header = f"{'mu1':>5} {'mu2':>5} {'lam':>5} {'rho1':>5} | {'T_bot':>7} {'T_sim':>7} {'T_approx':>8} {'Err%':>6} {'T_LH':>7} {'ErrLH%':>7} {'T_UB':>7}"
    print(header)
    print("-" * len(header))

    for mu1, mu2, lam in scenarios:
        rho1 = lam / mu1
        t_bot = lower_bound_bottleneck(lam, mu1, mu2)
        t_ub = upper_bound_independent(lam, mu1, mu2)
        t_approx = mean_response_time(lam, mu1, mu2)
        t_lh = mean_response_time_lh(lam, mu1, mu2)
        res = simulate(lam, mu1, mu2, n_jobs=2_000_000, warmup=100_000, seed=42)
        t_sim = res.mean_response_time
        error = (t_approx - t_sim) / t_sim * 100
        error_lh = (t_lh - t_sim) / t_sim * 100
        print(
            f"{mu1:5.1f} {mu2:5.1f} {lam:5.1f} {rho1:5.2f} | "
            f"{t_bot:7.3f} {t_sim:7.3f} {t_approx:8.3f} {error:+6.2f}% "
            f"{t_lh:7.3f} {error_lh:+7.2f}% {t_ub:7.3f}"
        )


def print_lh_comparison_table():
    """Print LH vs enhanced LH vs simulation for the same parameter sets."""
    scenarios = [
        (1.0, 1.0, 0.3),
        (1.0, 1.0, 0.6),
        (1.0, 1.0, 0.9),
        (1.0, 1.5, 0.3),
        (1.0, 1.5, 0.6),
        (1.0, 1.5, 0.9),
        (1.0, 2.0, 0.3),
        (1.0, 2.0, 0.6),
        (1.0, 2.0, 0.9),
        (1.0, 3.0, 0.6),
        (1.0, 5.0, 0.6),
    ]

    header = f"{'mu1':>5} {'mu2':>5} {'lam':>5} {'rho1':>5} | {'T_sim':>7} {'T_LH':>7} {'ErrLH%':>7} {'T_LHenh':>8} {'ErrEnh%':>8}"
    print(header)
    print("-" * len(header))

    for mu1, mu2, lam in scenarios:
        rho1 = lam / mu1
        t_lh = mean_response_time_lh(lam, mu1, mu2)
        t_enh = mean_response_time_lh_enhanced(lam, mu1, mu2)
        res = simulate(lam, mu1, mu2, n_jobs=2_000_000, warmup=100_000, seed=42)
        t_sim = res.mean_response_time
        err_lh = (t_lh - t_sim) / t_sim * 100
        err_enh = (t_enh - t_sim) / t_sim * 100
        print(
            f"{mu1:5.1f} {mu2:5.1f} {lam:5.1f} {rho1:5.2f} | "
            f"{t_sim:7.3f} {t_lh:7.3f} {err_lh:+7.2f}% "
            f"{t_enh:8.3f} {err_enh:+8.2f}%"
        )


def main():
    print("=" * 60)
    print("Fork-Join Queue: Approximation vs Simulation")
    print("=" * 60)
    print()
    print_comparison_table()

    print()
    print("LH vs Enhanced LH vs Simulation")
    print("-" * 60)
    print_lh_comparison_table()

    print("\nGenerating plots...")

    fig1, _ = plot_vs_load(1.0, 2.0, lam_points=50, run_simulation=True)
    fig1.savefig(OUT_DIR / "response_time_vs_load.png", dpi=150)
    print(f"  Saved {OUT_DIR / 'response_time_vs_load.png'}")

    mu2_range = np.linspace(1.1, 5.0, 40)
    fig2, _ = plot_vs_heterogeneity(1.0, 0.6, mu2_range, run_simulation=True)
    fig2.savefig(OUT_DIR / "response_time_vs_heterogeneity.png", dpi=150)
    print(f"  Saved {OUT_DIR / 'response_time_vs_heterogeneity.png'}")

    print("\nDone.")


if __name__ == "__main__":
    main()
