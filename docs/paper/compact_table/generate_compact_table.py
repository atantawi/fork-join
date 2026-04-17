#!/usr/bin/env python3
"""Generate compact LaTeX comparison table for T_UL, T_LH, T_LHe vs simulation.

Grid: mu1=1, r in {1,2,4,8}, rho in {0.4,0.8,0.9,0.95}.
Simulation: 10M jobs, warmup 500K, seed 42, normal-approximation 95% CI.
Output: compact_table.tex in this directory.
Cache: sim_cache.json avoids re-running completed simulations.
"""

import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

from forkjoin.analytical import (
    mean_response_time,
    mean_response_time_lh,
    mean_response_time_lh_enhanced,
)
from forkjoin.simulation import simulate

MU1 = 1.0
R_VALUES = [1, 2, 4, 8]
RHO_VALUES = [0.4, 0.8, 0.9, 0.95]
N_JOBS = 10_000_000
WARMUP = 500_000
SEED = 42
BETA = 10.0

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_FILE = os.path.join(SCRIPT_DIR, "sim_cache.json")
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "compact_table.tex")


def load_cache():
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE) as f:
            return json.load(f)
    return {}


def save_cache(cache):
    with open(CACHE_FILE, "w") as f:
        json.dump(cache, f, indent=2)


def run_or_load(rho, r, cache):
    key = f"{rho},{r}"
    if key in cache:
        d = cache[key]
        ci_hw = (d["ci_95"][1] - d["ci_95"][0]) / 2
        return d["mean_response_time"], ci_hw
    lam = rho * MU1
    mu2 = r * MU1
    print(f"  Simulating rho={rho}, r={r}  (lam={lam}, mu2={mu2}) ...")
    res = simulate(lam, MU1, mu2, n_jobs=N_JOBS, warmup=WARMUP, seed=SEED)
    cache[key] = {
        "mean_response_time": res.mean_response_time,
        "std_response_time": res.std_response_time,
        "n_samples": res.n_samples,
        "ci_95": list(res.ci_95),
    }
    ci_hw = (res.ci_95[1] - res.ci_95[0]) / 2
    return res.mean_response_time, ci_hw


def fmt_t(val):
    return f"{val:.3f}"


def fmt_err(val):
    sign = "+" if val >= 0 else ""
    return f"${sign}{val:.2f}\\%$"


def fmt_sim(t, hw):
    return f"${t:.3f} \\pm {hw:.3f}$"


def generate_table(data):
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        (
            r"\caption{Approximations $T_\mathrm{UL}$, $T_\mathrm{LH}$, and "
            r"$T_\mathrm{LH}^\mathrm{enh}$ vs.\ simulation for the heterogeneous "
            r"two-queue fork-join ($\mu_1 = 1$, $\beta = 10$). "
            r"$T_\mathrm{sim} \pm 95\%\,\mathrm{CI}$ from $10\,\mathrm{M}$ jobs "
            r"(seed\,42, normal approximation). "
            r"Errors are $(T_\mathrm{approx} - T_\mathrm{sim})/T_\mathrm{sim} \times 100\%$.}"
        ),
        r"\label{tab:compact}",
        r"\resizebox{\textwidth}{!}{%",
        r"\begin{tabular}{@{}ccrcccccc@{}}",
        r"\toprule",
        (
            r"& & \multicolumn{1}{c}{Simulation}"
            r" & \multicolumn{2}{c}{$T_\mathrm{UL}$}"
            r" & \multicolumn{2}{c}{$T_\mathrm{LH}$}"
            r" & \multicolumn{2}{c}{$T_\mathrm{LH}^\mathrm{enh}$} \\"
        ),
        r"\cmidrule(lr){3-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}\cmidrule(lr){8-9}",
        (
            r"$\rho$ & $r$"
            r" & $T_\mathrm{sim} \pm 95\%\,\mathrm{CI}$"
            r" & $T$ & Err"
            r" & $T$ & Err"
            r" & $T$ & Err \\"
        ),
        r"\midrule",
    ]

    for i, rho in enumerate(RHO_VALUES):
        if i > 0:
            lines.append(r"\midrule")
        for j, r in enumerate(R_VALUES):
            d = data[(rho, r)]
            rho_cell = rf"\multirow{{4}}{{*}}{{${rho}$}}" if j == 0 else ""
            row = (
                f"{rho_cell} & ${r}$"
                f" & {fmt_sim(d['t_sim'], d['ci_hw'])}"
                f" & {fmt_t(d['t_ul'])} & {fmt_err(d['err_ul'])}"
                f" & {fmt_t(d['t_lh'])} & {fmt_err(d['err_lh'])}"
                f" & {fmt_t(d['t_lhe'])} & {fmt_err(d['err_lhe'])} \\\\"
            )
            lines.append(row)

    lines += [
        r"\bottomrule",
        r"\end{tabular}%",
        r"}",
        r"\end{table}",
    ]
    return "\n".join(lines)


def main():
    cache = load_cache()
    data = {}

    print("Running simulations (cached results reused where available)...")
    for rho in RHO_VALUES:
        for r in R_VALUES:
            lam = rho * MU1
            mu2 = r * MU1
            t_sim, ci_hw = run_or_load(rho, r, cache)
            t_ul = mean_response_time(lam, MU1, mu2)
            t_lh = mean_response_time_lh(lam, MU1, mu2, beta=BETA)
            t_lhe = mean_response_time_lh_enhanced(lam, MU1, mu2, beta=BETA)
            data[(rho, r)] = {
                "t_sim": t_sim,
                "ci_hw": ci_hw,
                "t_ul": t_ul,
                "err_ul": (t_ul - t_sim) / t_sim * 100,
                "t_lh": t_lh,
                "err_lh": (t_lh - t_sim) / t_sim * 100,
                "t_lhe": t_lhe,
                "err_lhe": (t_lhe - t_sim) / t_sim * 100,
            }

    save_cache(cache)

    tex = generate_table(data)
    with open(OUTPUT_FILE, "w") as f:
        f.write(tex + "\n")
    print(f"Table written to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
