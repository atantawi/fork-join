"""Visualization for fork-join queue analysis."""

import numpy as np
import matplotlib.pyplot as plt

from . import analytical
from .simulation import simulate


def plot_vs_load(mu1, mu2, lam_points=50, run_simulation=True, n_jobs=500_000):
    """Plot response time quantities vs load (lambda / min(mu1, mu2)).

    Returns (fig, ax).
    """
    mu_min = min(mu1, mu2)
    loads = np.linspace(0.05, 0.95, lam_points)
    lams = loads * mu_min

    t_approx = [analytical.mean_response_time(l, mu1, mu2) for l in lams]
    t_lh = [analytical.mean_response_time_lh(l, mu1, mu2) for l in lams]
    t_ub = [analytical.upper_bound_independent(l, mu1, mu2) for l in lams]
    t_bot = [analytical.lower_bound_bottleneck(l, mu1, mu2) for l in lams]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(loads, t_approx, "b-", linewidth=2, label="Approximation")
    ax.plot(loads, t_lh, "m-", linewidth=2, label="LH interpolation")
    ax.plot(loads, t_ub, "r--", linewidth=1.5, label="Upper bound (indep.)")
    ax.plot(loads, t_bot, "g--", linewidth=1.5, label="Lower bound (bottleneck)")

    if run_simulation:
        sim_loads = loads[::max(1, lam_points // 10)]
        sim_lams = sim_loads * mu_min
        t_sim = []
        for l in sim_lams:
            res = simulate(l, mu1, mu2, n_jobs=n_jobs, seed=42)
            t_sim.append(res.mean_response_time)
        ax.plot(sim_loads, t_sim, "ko", markersize=5, label="Simulation")

    ax.set_xlabel(r"Load $\lambda / \min(\mu_1, \mu_2)$")
    ax.set_ylabel("Mean response time $T$")
    ax.set_title(f"Fork-Join Response Time ($\\mu_1={mu1}$, $\\mu_2={mu2}$)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig, ax


def plot_vs_heterogeneity(mu1, lam, mu2_range, run_simulation=True, n_jobs=500_000):
    """Plot response time quantities vs heterogeneity ratio mu2/mu1.

    mu2_range: array-like of mu2 values to plot.
    Returns (fig, ax).
    """
    mu2_arr = np.asarray(mu2_range, dtype=float)
    ratios = mu2_arr / mu1

    t_approx = [analytical.mean_response_time(lam, mu1, m) for m in mu2_arr]
    t_lh = [analytical.mean_response_time_lh(lam, mu1, m) for m in mu2_arr]
    t_ub = [analytical.upper_bound_independent(lam, mu1, m) for m in mu2_arr]
    t_bot = [analytical.lower_bound_bottleneck(lam, mu1, m) for m in mu2_arr]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(ratios, t_approx, "b-", linewidth=2, label="Approximation")
    ax.plot(ratios, t_lh, "m-", linewidth=2, label="LH interpolation")
    ax.plot(ratios, t_ub, "r--", linewidth=1.5, label="Upper bound (indep.)")
    ax.plot(ratios, t_bot, "g--", linewidth=1.5, label="Lower bound (bottleneck)")

    if run_simulation:
        step = max(1, len(mu2_arr) // 10)
        sim_mu2 = mu2_arr[::step]
        sim_ratios = sim_mu2 / mu1
        t_sim = []
        for m in sim_mu2:
            res = simulate(lam, mu1, m, n_jobs=n_jobs, seed=42)
            t_sim.append(res.mean_response_time)
        ax.plot(sim_ratios, t_sim, "ko", markersize=5, label="Simulation")

    ax.set_xlabel(r"Heterogeneity ratio $\mu_2 / \mu_1$")
    ax.set_ylabel("Mean response time $T$")
    ax.set_title(f"Fork-Join Response Time ($\\mu_1={mu1}$, $\\lambda={lam}$)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig, ax
