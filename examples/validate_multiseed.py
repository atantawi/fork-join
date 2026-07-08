"""Multi-seed validation harness for the heterogeneous 2-queue fork-join.

Reproduces the protocol described in docs/paper/sections/validation.tex
(Tables 1 and 2): for each scenario, run N_JOBS jobs after a WARMUP-job
warmup, repeated over N_SEEDS independent seeds. The reported T_sim is the
grand mean over seeds; the 95% CI half-width is t_{0.975, df} * s / sqrt(k),
where s is the sample std (ddof=1) of the per-seed means and df = k - 1.

This is the "independent replications" method: the between-seed variance is
the correct basis for a CI, unlike a single run's normal approximation, which
ignores the strong autocorrelation between consecutive response times near
saturation.

Each replication is a single long run of the library simulate() (one continuous
Lindley recursion, warmup discarded) -- the same routine used everywhere else in
the repo. We loop it over seeds here; we do NOT reimplement the simulator.

Usage:
    python examples/validate_multiseed.py            # 20M jobs, 5 seeds (paper protocol; ~15 min)
    python examples/validate_multiseed.py --quick    # 200K jobs, 5 seeds (smoke test)
"""

import argparse
import json
import math
import os
import time

from forkjoin import (
    simulate,
    mean_response_time,           # T_UL
    mean_response_time_lh,        # T_LH
    mean_response_time_lh_enhanced,  # T_LHe
    lower_bound_bottleneck,       # T_bot
)

# Scenarios from Table 1 / Table 2 of validation.tex: (mu1, mu2, lambda).
SCENARIOS = [
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

# Paper's current numbers (T_sim, CI half-width) for side-by-side comparison.
# Keyed by (mu1, mu2, lam). Taken from Table 2 (more decimals) where available.
PAPER = {
    (1.0, 1.0, 0.3): (2.089, 0.001),
    (1.0, 1.0, 0.6): (3.561, 0.002),
    (1.0, 1.0, 0.9): (13.855, 0.061),
    (1.0, 1.5, 0.3): (1.707, 0.001),
    (1.0, 1.5, 0.6): (2.767, 0.001),
    (1.0, 1.5, 0.9): (10.120, 0.024),
    (1.0, 2.0, 0.3): (1.583, 0.001),
    (1.0, 2.0, 0.6): (2.623, 0.001),
    (1.0, 2.0, 0.9): (10.028, 0.024),
    (1.0, 3.0, 0.6): (2.547, 0.001),
    (1.0, 5.0, 0.6): (2.515, 0.001),
}

# Student-t 0.975 quantiles by degrees of freedom (avoids a scipy dependency).
T_975 = {1: 12.706205, 2: 4.302653, 3: 3.182446, 4: 2.776445, 5: 2.570582,
         6: 2.446912, 7: 2.364624, 8: 2.306004, 9: 2.262157}

CACHE_FILE = os.path.join(os.path.dirname(__file__), "validate_multiseed_cache.json")


def grand_mean_and_ci(per_seed_means):
    """Grand mean and 95% t-CI half-width across independent replications."""
    k = len(per_seed_means)
    gm = sum(per_seed_means) / k
    var = sum((m - gm) ** 2 for m in per_seed_means) / (k - 1)
    s = math.sqrt(var)
    half = T_975[k - 1] * s / math.sqrt(k)
    return gm, half


def run(n_jobs, warmup, seeds, cache):
    rows = []
    for (mu1, mu2, lam) in SCENARIOS:
        key = f"{mu1},{mu2},{lam}"
        if key in cache and cache[key]["n_jobs"] == n_jobs \
                and cache[key]["seeds"] == list(seeds):
            per_seed = cache[key]["per_seed_means"]
            print(f"[cache] mu2={mu2} lam={lam}")
        else:
            print(f"[run]   mu2={mu2} lam={lam}: {len(seeds)} seeds x {n_jobs:,} jobs ...",
                  flush=True)
            t0 = time.time()
            per_seed = [
                simulate(lam, mu1, mu2, n_jobs=n_jobs, warmup=warmup, seed=s).mean_response_time
                for s in seeds
            ]
            cache[key] = {"n_jobs": n_jobs, "warmup": warmup,
                          "seeds": list(seeds), "per_seed_means": per_seed}
            with open(CACHE_FILE, "w") as f:
                json.dump(cache, f, indent=2)
            print(f"        done in {time.time() - t0:.0f}s  per-seed={['%.4f' % m for m in per_seed]}",
                  flush=True)

        t_sim, ci = grand_mean_and_ci(per_seed)
        t_ul = mean_response_time(lam, mu1, mu2)
        t_lh = mean_response_time_lh(lam, mu1, mu2)
        t_lhe = mean_response_time_lh_enhanced(lam, mu1, mu2)
        t_bot = lower_bound_bottleneck(lam, mu1, mu2)
        rows.append({
            "mu1": mu1, "mu2": mu2, "lam": lam, "t_bot": t_bot,
            "t_sim": t_sim, "ci": ci,
            "t_ul": t_ul, "err_ul": (t_ul - t_sim) / t_sim * 100,
            "t_lh": t_lh, "err_lh": (t_lh - t_sim) / t_sim * 100,
            "t_lhe": t_lhe, "err_lhe": (t_lhe - t_sim) / t_sim * 100,
        })
    return rows


def print_report(rows, n_jobs, seeds):
    df = len(seeds) - 1
    print("\n" + "=" * 108)
    print(f"REGENERATED  ({len(seeds)} seeds x {n_jobs:,} jobs, t-CI df={df})   "
          f"vs  PAPER (validation.tex)")
    print("=" * 108)
    hdr = (f"{'mu1':>3} {'mu2':>4} {'lam':>4} | "
           f"{'T_sim(new)':>10} {'+-CI':>7} | {'T_sim(paper)':>12} {'+-CI':>6} | "
           f"{'dT_sim':>7} | {'errUL%':>7} {'errLH%':>7} {'errLHe%':>7}")
    print(hdr)
    print("-" * 108)
    for r in rows:
        p_ts, p_ci = PAPER[(r["mu1"], r["mu2"], r["lam"])]
        d = r["t_sim"] - p_ts
        flag = "  <-- shifted" if abs(d) > max(r["ci"], p_ci) else ""
        print(f"{r['mu1']:>3} {r['mu2']:>4} {r['lam']:>4} | "
              f"{r['t_sim']:>10.4f} {r['ci']:>7.4f} | {p_ts:>12.3f} {p_ci:>6.3f} | "
              f"{d:>+7.4f} | {r['err_ul']:>+7.2f} {r['err_lh']:>+7.2f} {r['err_lhe']:>+7.2f}{flag}")
    print("=" * 108)
    print("dT_sim = new grand mean - paper value. '<-- shifted' = gap exceeds both CIs.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="200K jobs (smoke test)")
    ap.add_argument("--n-jobs", type=int, default=None)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    args = ap.parse_args()

    n_jobs = args.n_jobs or (200_000 if args.quick else 20_000_000)
    warmup = 500_000 if not args.quick else 10_000
    seeds = args.seeds

    cache = {}
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE) as f:
            cache = json.load(f)

    print(f"Protocol: {len(seeds)} seeds {seeds}, n_jobs={n_jobs:,}, warmup={warmup:,}")
    rows = run(n_jobs, warmup, seeds, cache)
    print_report(rows, n_jobs, seeds)


if __name__ == "__main__":
    main()
