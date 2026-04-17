"""Visualization for fork-join queue analysis."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

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
    t_lhe = [analytical.mean_response_time_lh_enhanced(l, mu1, mu2) for l in lams]
    t_ub = [analytical.upper_bound_independent(l, mu1, mu2) for l in lams]
    t_bot = [analytical.lower_bound_bottleneck(l, mu1, mu2) for l in lams]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(loads, t_approx, "b-", linewidth=2, label=r"$T_{\mathrm{UL}}$")
    ax.plot(loads, t_lh, "m-", linewidth=2, label=r"$T_{\mathrm{LH}}$")
    ax.plot(loads, t_lhe, "c-", linewidth=2, label=r"$T_{\mathrm{LH}}^{\mathrm{enh}}$")
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
    t_lhe = [analytical.mean_response_time_lh_enhanced(lam, mu1, m) for m in mu2_arr]
    t_ub = [analytical.upper_bound_independent(lam, mu1, m) for m in mu2_arr]
    t_bot = [analytical.lower_bound_bottleneck(lam, mu1, m) for m in mu2_arr]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(ratios, t_approx, "b-", linewidth=2, label=r"$T_{\mathrm{UL}}$")
    ax.plot(ratios, t_lh, "m-", linewidth=2, label=r"$T_{\mathrm{LH}}$")
    ax.plot(ratios, t_lhe, "c-", linewidth=2, label=r"$T_{\mathrm{LH}}^{\mathrm{enh}}$")
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


def plot_t_ul_vs_heterogeneity_panel(
    rho_values=(0.52, 0.76, 0.88, 0.94),
    r_values=(1.0, 1.15, 1.3, 1.5, 2.0, 3.0, 5.0, 8.0),
    *,
    mu1=1.0,
    n_analytical=100,
    sim_results=None,
    n_jobs=1_000_000,
    warmup=100_000,
    seed=42,
    show_bounds=True,
    err_ylim=(-1.5, 1.5),
    figsize=(10, 8),
):
    """2x2 panel: T_UL vs simulation across heterogeneity ratio r for 4 load levels.

    Each subplot shows:
      - Left y-axis: smooth T_UL curve + simulation markers with 95% CI bars
      - Right y-axis (red): relative error (T_UL - T_sim) / T_sim * 100%
      - x-axis: r = mu_2 / mu_1 on log scale

    Args:
        rho_values: 4 load levels (rho = lambda / mu1). One subplot each.
        r_values: Heterogeneity ratios at which simulation is run.
        mu1: Bottleneck service rate. mu2 = r * mu1, lam = rho * mu1.
        n_analytical: Number of r points for smooth analytical curves.
        sim_results: Optional dict mapping (rho, r) -> SimResult. If None,
            simulations are run using n_jobs/warmup/seed.
        n_jobs: Jobs per simulation run (when sim_results is None).
        warmup: Warmup jobs discarded from statistics.
        seed: Random seed for reproducibility.
        show_bounds: Whether to overlay T_UB and T_bot as light dashed lines.
        figsize: Figure size in inches.

    Returns:
        (fig, axes) where axes is the 2x2 ndarray of Axes.
    """
    rho_values = list(rho_values)
    r_values = list(r_values)
    r_smooth = np.geomspace(r_values[0], r_values[-1], n_analytical)

    fig, axes = plt.subplots(2, 2, figsize=figsize, constrained_layout=True)

    # Colors
    col_tul = "#1f77b4"       # blue for T_UL curve
    col_sim = "#0d47a1"       # dark blue for simulation markers
    col_err = "#e53935"       # red for error axis
    col_bound = "#90a4ae"     # gray for bounds

    legend_handles = []
    legend_labels = []

    for idx, (ax, rho) in enumerate(zip(axes.flat, rho_values)):
        lam = rho * mu1

        # --- Smooth analytical curves ---
        t_ul_smooth = [
            analytical.mean_response_time(lam, mu1, r * mu1) for r in r_smooth
        ]
        h_tul, = ax.plot(
            r_smooth, t_ul_smooth, "-", color=col_tul, linewidth=2,
            label=r"$T_{\mathrm{UL}}$",
        )

        if show_bounds:
            t_ub_smooth = [
                analytical.upper_bound_independent(lam, mu1, r * mu1) for r in r_smooth
            ]
            t_bot_smooth = [
                analytical.lower_bound_bottleneck(lam, mu1, r * mu1) for r in r_smooth
            ]
            h_ub, = ax.plot(
                r_smooth, t_ub_smooth, "--", color=col_bound, linewidth=1, alpha=0.5,
                label=r"$T_{\mathrm{UB}}$",
            )
            h_bot, = ax.plot(
                r_smooth, t_bot_smooth, ":", color=col_bound, linewidth=1, alpha=0.5,
                label=r"$T_{\mathrm{bot}}$",
            )

        # --- Simulation points ---
        t_ul_pts = []
        t_sim_pts = []
        ci_lo = []
        ci_hi = []

        for r in r_values:
            mu2 = r * mu1
            t_ul_pts.append(analytical.mean_response_time(lam, mu1, mu2))
            if sim_results is not None and (rho, r) in sim_results:
                res = sim_results[(rho, r)]
            else:
                res = simulate(lam, mu1, mu2, n_jobs=n_jobs, warmup=warmup, seed=seed)
            t_sim_pts.append(res.mean_response_time)
            ci_lo.append(res.mean_response_time - res.ci_95[0])
            ci_hi.append(res.ci_95[1] - res.mean_response_time)

        h_sim = ax.errorbar(
            r_values, t_sim_pts, yerr=[ci_lo, ci_hi],
            fmt="o", color=col_sim, markersize=6, capsize=3, linewidth=1,
            label="Simulation (95% CI)",
        )

        # --- Right y-axis: relative error ---
        ax2 = ax.twinx()
        rel_err = [
            (ul - sim) / sim * 100
            for ul, sim in zip(t_ul_pts, t_sim_pts)
        ]
        h_err, = ax2.plot(
            r_values, rel_err, "^:", color=col_err, markersize=5, linewidth=1,
            alpha=0.85, label="Rel. error (%)",
        )
        ax2.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.4)
        if err_ylim is not None:
            ax2.set_ylim(err_ylim)
        ax2.yaxis.label.set_color(col_err)
        ax2.tick_params(axis="y", colors=col_err)
        ax2.spines["right"].set_color(col_err)
        if idx % 2 == 1:  # right column only
            ax2.set_ylabel("Relative error (%)", color=col_err)
        else:
            ax2.set_ylabel("")

        # --- x-axis formatting ---
        ax.set_xscale("log")
        ax.set_xticks(r_values)
        ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.xaxis.set_minor_formatter(ticker.NullFormatter())
        if idx >= 2:  # bottom row only
            ax.set_xlabel(r"$r = \mu_2 / \mu_1$")

        # --- Left y-axis ---
        if idx % 2 == 0:  # left column only
            ax.set_ylabel("Mean response time $T$")
        ax.grid(True, alpha=0.2, which="both")

        # --- Panel title ---
        ax.set_title(fr"$\rho = {rho}$", fontsize=11)

        # Collect legend items from first panel only
        if idx == 0:
            legend_handles.append(h_tul)
            legend_labels.append(r"$T_{\mathrm{UL}}$")
            legend_handles.append(h_sim)
            legend_labels.append("Simulation (95% CI)")
            legend_handles.append(h_err)
            legend_labels.append("Rel. error (%)")
            if show_bounds:
                legend_handles.append(h_ub)
                legend_labels.append(r"$T_{\mathrm{UB}}$")
                legend_handles.append(h_bot)
                legend_labels.append(r"$T_{\mathrm{bot}}$")

    fig.legend(
        legend_handles, legend_labels,
        loc="lower center", ncol=len(legend_handles),
        bbox_to_anchor=(0.5, -0.04), frameon=True, fontsize=9,
    )

    return fig, axes
