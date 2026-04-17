#!/usr/bin/env python3
"""Generate 2x2 T_UL validation panel: T_UL vs heterogeneity ratio r.

Usage:
    python generate_t_ul_plots.py               # 10M jobs (production)
    python generate_t_ul_plots.py --quick        # 100K jobs (quick visual check)
    python generate_t_ul_plots.py --clear-cache  # delete cache and re-run
"""

import argparse
import json
import pathlib
import time

from forkjoin import simulate
from forkjoin.plotting import plot_t_ul_vs_heterogeneity_panel
from forkjoin.simulation import SimResult

SCRIPT_DIR = pathlib.Path(__file__).parent
PAPER_FIG_DIR = SCRIPT_DIR.parent / "docs" / "paper" / "figures"
CACHE_FILE = SCRIPT_DIR / "t_ul_sim_cache.json"

RHO_VALUES = (0.52, 0.76, 0.88, 0.94)
R_VALUES = (1.0, 1.15, 1.3, 1.5, 2.0, 3.0, 5.0, 8.0)
MU1 = 1.0


def _load_cache(path):
    if not path.exists():
        return {}
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


def _save_cache(path, results):
    raw = {}
    for (rho, r), res in results.items():
        raw[f"{rho},{r}"] = {
            "mean_response_time": res.mean_response_time,
            "std_response_time": res.std_response_time,
            "n_samples": res.n_samples,
            "ci_95": list(res.ci_95),
        }
    with open(path, "w") as f:
        json.dump(raw, f, indent=2)


def run_simulations(rho_values, r_values, mu1, n_jobs, warmup, seed, cache_path):
    cache = _load_cache(cache_path)
    results = dict(cache)
    total = len(rho_values) * len(r_values)
    done = sum(1 for k in results if k[0] in rho_values and k[1] in r_values)
    skipped = 0

    for i, rho in enumerate(rho_values):
        lam = rho * mu1
        for j, r in enumerate(r_values):
            key = (rho, r)
            if key in cache:
                skipped += 1
                continue
            mu2 = r * mu1
            count = i * len(r_values) + j + 1
            print(
                f"  [{count}/{total}] rho={rho}, r={r} (lam={lam:.4f}, mu2={mu2:.4f})"
                f" ... ",
                end="", flush=True,
            )
            t0 = time.time()
            res = simulate(lam, mu1, mu2, n_jobs=n_jobs, warmup=warmup, seed=seed)
            elapsed = time.time() - t0
            ci_half = (res.ci_95[1] - res.ci_95[0]) / 2
            print(
                f"{elapsed:.1f}s  T={res.mean_response_time:.4f}"
                f"  ±{ci_half:.5f} (95% CI)"
            )
            results[key] = res

    if skipped:
        print(f"  ({skipped}/{total} loaded from cache)")

    _save_cache(cache_path, results)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quick", action="store_true",
        help="Use 100K jobs for a fast visual check",
    )
    parser.add_argument(
        "--clear-cache", action="store_true",
        help="Delete cached simulation results before running",
    )
    parser.add_argument(
        "--n-jobs", type=int, default=None,
        help="Override number of simulation jobs",
    )
    args = parser.parse_args()

    if args.clear_cache and CACHE_FILE.exists():
        CACHE_FILE.unlink()
        print(f"Cleared cache: {CACHE_FILE}")

    n_jobs = args.n_jobs or (100_000 if args.quick else 10_000_000)
    warmup = max(10_000, n_jobs // 20)
    print(f"Simulation: n_jobs={n_jobs:,}, warmup={warmup:,}, seed=42")
    print(f"Grid: {len(RHO_VALUES)} rho values x {len(R_VALUES)} r values"
          f" = {len(RHO_VALUES) * len(R_VALUES)} simulations")
    print()

    sim_results = run_simulations(
        RHO_VALUES, R_VALUES, MU1, n_jobs, warmup, seed=42, cache_path=CACHE_FILE,
    )

    print("\nGenerating panel plot...")
    fig, _ = plot_t_ul_vs_heterogeneity_panel(
        rho_values=RHO_VALUES,
        r_values=R_VALUES,
        mu1=MU1,
        sim_results=sim_results,
        show_bounds=True,
    )

    dpi = 150 if args.quick else 300
    outputs = [SCRIPT_DIR / "t_ul_vs_heterogeneity.png"]
    if not args.quick:
        outputs += [
            PAPER_FIG_DIR / "t_ul_vs_heterogeneity.png",
            PAPER_FIG_DIR / "t_ul_vs_heterogeneity.pdf",
        ]

    for dest in outputs:
        fig.savefig(dest, dpi=dpi, bbox_inches="tight")
        print(f"  Saved {dest}")

    print("\nDone.")


if __name__ == "__main__":
    main()
