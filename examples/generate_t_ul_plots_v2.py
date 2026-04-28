#!/usr/bin/env python3
"""Generate single-panel T_UL validation plot: four rho curves on one axes.

Loads simulation results from the existing cache (t_ul_sim_cache.json).
Only T_UL and simulation points are shown (no bounds, no relative error).

Usage:
    python generate_t_ul_plots_v2.py
"""

import json
import pathlib

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
import numpy as np

from forkjoin import analytical
from forkjoin.simulation import SimResult

SCRIPT_DIR = pathlib.Path(__file__).parent
CACHE_FILE = SCRIPT_DIR / "t_ul_sim_cache.json"

RHO_VALUES = (0.52, 0.76, 0.88, 0.94)
R_VALUES   = (1.0, 1.15, 1.3, 1.5, 2.0, 3.0, 5.0, 8.0)
MU1        = 1.0

COLORS  = ["steelblue", "darkorange", "forestgreen", "crimson"]
MARKERS = ["o", "D", "s", "^"]
R_SMOOTH = np.geomspace(R_VALUES[0], R_VALUES[-1], 300)


def load_cache(path):
    if not path.exists():
        raise FileNotFoundError(f"Cache not found: {path}. Run generate_t_ul_plots.py first.")
    with open(path) as f:
        raw = json.load(f)
    result = {}
    for key_str, v in raw.items():
        rho_s, r_s = key_str.split(",")
        result[(float(rho_s), float(r_s))] = SimResult(
            mean_response_time=v["mean_response_time"],
            std_response_time=v["std_response_time"],
            n_samples=v["n_samples"],
            ci_95=(v["ci_95"][0], v["ci_95"][1]),
        )
    return result


def main():
    sim_results = load_cache(CACHE_FILE)

    fig, ax = plt.subplots(figsize=(7, 5))

    for rho, color, marker in zip(RHO_VALUES, COLORS, MARKERS):
        lam = rho * MU1

        # Smooth analytical curve
        t_ul_smooth = [analytical.mean_response_time(lam, MU1, r * MU1) for r in R_SMOOTH]
        ax.plot(R_SMOOTH, t_ul_smooth, "-", color=color, linewidth=2, zorder=2)

        # Simulation markers with 95% CI
        t_sim = []
        ci_lo = []
        ci_hi = []
        for r in R_VALUES:
            res = sim_results[(rho, r)]
            t_sim.append(res.mean_response_time)
            ci_lo.append(res.mean_response_time - res.ci_95[0])
            ci_hi.append(res.ci_95[1] - res.mean_response_time)

        ax.errorbar(
            R_VALUES, t_sim, yerr=[ci_lo, ci_hi],
            fmt=marker, color=color, markersize=6, capsize=3, linewidth=1, zorder=3,
        )

    ax.set_xscale("log")
    ax.set_xticks(R_VALUES)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_ylim(0, 24)
    ax.grid(True, alpha=0.3, which="both")
    ax.set_xlabel(r"$r = \hat{\mu}_2 / \hat{\mu}_1$", fontsize=11)
    ax.set_ylabel(r"$\mathbb{E}[T_{FJ}]$", fontsize=11)

    legend_handles = [
        Line2D([0], [0], color=c, linewidth=2, marker=m, markersize=6,
               label=rf"$\rho = {rho}$")
        for rho, c, m in zip(RHO_VALUES, COLORS, MARKERS)
    ]
    ax.legend(handles=legend_handles, fontsize=10, loc="upper right")

    for out in [
        SCRIPT_DIR / "t_ul_vs_heterogeneity_v2.png",
        SCRIPT_DIR / "t_ul_vs_heterogeneity_v2.pdf",
    ]:
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Saved {out}")
    plt.show()


if __name__ == "__main__":
    main()
