#!/usr/bin/env python3
"""Plot the eq.23 vs eq.14 comparison across the Table 1 grid.

Companion figure for eq23_vs_eq14_comparison.md (issue #6, Finding 2).

Two panels, both analytical (no simulation needed for the expressions):
  (A) Pairwise relative gap  (T_LHe[eq.23] - T_FJ1[eq.14]) / T_FJ1  vs rho,
      one line per heterogeneity ratio r. Shows the gap = c2*rho^2/(mu_min(1-rho))
      growing toward saturation and peaking at intermediate r (=4).
  (B) Signed relative error vs simulation for both expressions: eq.23 (solid)
      and eq.14 (dashed), colored by r. Shows eq.14 systematically
      under-predicting at moderate-to-heavy load.

Color encodes the entity (r), never rank; line style encodes the expression.
T_sim is the cached multi-replica / 20M-job grand mean (correct t-CI protocol;
currently 10 replicas -- whatever reproduce_table1.py has cached is used).

Usage:  python generate_eq23_vs_eq14_plot.py
Outputs: eq23_vs_eq14.png, eq23_vs_eq14.pdf
"""

import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from forkjoin import mean_response_time_lh_enhanced  # eq.23 (quadratic)

MU1 = 1.0
BETA = 10.0
RHO_VALUES = [0.4, 0.8, 0.9, 0.95]
R_VALUES = [2, 4, 8]
# Fixed categorical order for r (matches sibling figure table1_errors.png).
R_COLORS = {2: "steelblue", 4: "darkorange", 8: "forestgreen"}

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_FILE = os.path.join(SCRIPT_DIR, "table1_results.json")


def eq14_first_order(lam, mu1, mu2, beta=BETA):
    """eq.14 == eq.23 with c2 set to zero (verified algebraically)."""
    mn, mx = min(mu1, mu2), max(mu1, mu2)
    t0 = 1 / mu1 + 1 / mu2 - 1 / (mu1 + mu2)
    f1 = 1 / mn**2 + 1 / mx**2 - 2 / (mn + mx) ** 2 - 2 * mn * mx / (mn + mx) ** 4
    c0 = mn * t0
    c1 = mn**2 * f1 - mn * t0
    rho = lam / mn
    return (c1 * rho + c0) / (mn * (1 - rho))


def load_tsim():
    """Grand mean over every cached 20M-job replication of each cell.

    reproduce_table1.py caches one entry per replication, keyed
    "<rho>,<r>|<n_jobs>|<warmup>|seed=<s>", so this picks up however many
    replicas have been run (currently 10) without a hard-coded seed list.
    """
    cache = json.load(open(CACHE_FILE))

    def tsim(rho, r):
        prefix = f"{rho},{r}|20000000|500000|seed="
        means = [v["mean"] for k, v in cache.items() if k.startswith(prefix)]
        if not means:
            raise KeyError(f"no cached 20M replications for rho={rho}, r={r}; "
                           "run reproduce_table1.py first")
        return sum(means) / len(means)

    return tsim


def n_replicas():
    cache = json.load(open(CACHE_FILE))
    return sum(1 for k in cache
               if k.startswith(f"{RHO_VALUES[0]},{R_VALUES[0]}|20000000|500000|seed="))


def main():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tsim = load_tsim()

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.2))

    # ---- Panel A: pairwise relative gap between the two expressions --------
    for r in R_VALUES:
        gaps = []
        for rho in RHO_VALUES:
            lam, mu2 = rho * MU1, r * MU1
            a = mean_response_time_lh_enhanced(lam, MU1, mu2, BETA)
            b = eq14_first_order(lam, MU1, mu2, BETA)
            gaps.append((a - b) / b * 100)
        axA.plot(RHO_VALUES, gaps, marker="o", markersize=6, linewidth=2,
                 color=R_COLORS[r], label=f"$r={r}$")
    axA.axhline(0, color="0.6", linewidth=0.8, zorder=0)
    axA.set_xlabel(r"System load  $\rho = \lambda/\mu_{\min}$")
    axA.set_ylabel("Relative gap  (eq.23 $-$ eq.14) / eq.14  [\\%]")
    axA.set_title("(A) Gap between the two expressions")
    axA.grid(True, linestyle="--", linewidth=0.5, color="0.85", zorder=0)
    axA.legend(title="heterogeneity", frameon=False)

    # ---- Panel B: signed error vs simulation for each expression -----------
    for r in R_VALUES:
        err23, err14 = [], []
        for rho in RHO_VALUES:
            lam, mu2 = rho * MU1, r * MU1
            ts = tsim(rho, r)
            a = mean_response_time_lh_enhanced(lam, MU1, mu2, BETA)
            b = eq14_first_order(lam, MU1, mu2, BETA)
            err23.append((a - ts) / ts * 100)
            err14.append((b - ts) / ts * 100)
        axB.plot(RHO_VALUES, err23, marker="o", markersize=6, linewidth=2,
                 color=R_COLORS[r], linestyle="-", label=f"$r={r}$")
        axB.plot(RHO_VALUES, err14, marker="s", markersize=6, linewidth=2,
                 color=R_COLORS[r], linestyle="--")
    axB.axhline(0, color="0.6", linewidth=0.8, zorder=0)
    axB.set_xlabel(r"System load  $\rho = \lambda/\mu_{\min}$")
    axB.set_ylabel("Relative error vs simulation  [\\%]")
    axB.set_title("(B) Accuracy vs simulation")
    axB.grid(True, linestyle="--", linewidth=0.5, color="0.85", zorder=0)
    # Two legends: color = r (entity), line style = expression.
    color_handles = [plt.Line2D([], [], color=R_COLORS[r], linewidth=2,
                                label=f"$r={r}$") for r in R_VALUES]
    style_handles = [
        plt.Line2D([], [], color="0.3", linewidth=2, linestyle="-",
                   marker="o", label=r"$T_{LHe}$ (eq.23, quadratic)"),
        plt.Line2D([], [], color="0.3", linewidth=2, linestyle="--",
                   marker="s", label=r"$T_{FJ}^{(1)}$ (eq.14, $c_2{=}0$)"),
    ]
    leg1 = axB.legend(handles=color_handles, title="heterogeneity",
                      frameon=False, loc="lower left")
    axB.add_artist(leg1)
    axB.legend(handles=style_handles, frameon=False, loc="upper right",
               fontsize=8)

    fig.suptitle(r"$T_{LHe}$ (eq.23) vs $T_{FJ}^{(1)}$ (eq.14) "
                 r"across the Table 1 grid  ($\mu_1{=}1,\ \mu_2{=}r,\ \beta{=}10$)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    for ext in ("png", "pdf"):
        out = os.path.join(SCRIPT_DIR, f"eq23_vs_eq14.{ext}")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
