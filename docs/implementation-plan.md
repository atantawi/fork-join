# Implementation Plan: Python Package for Fork-Join Queue Approximations

## Context

This research project has derived an approximate closed-form expression for the mean response time of a heterogeneous 2-queue M/M/1 fork-join system (documented in `docs/heterogeneous-fj-approximation.md`). The formula has been validated via simulation with <2% error. We now need a Python implementation providing the analytical formulas, a discrete-event simulator for validation, and visualization tools.

The repo currently has no Python code — only `docs/` and `CLAUDE.md`.

## Structure

```
fork-join/
├── forkjoin/
│   ├── __init__.py          # Package exports
│   ├── analytical.py        # All closed-form formulas (approximation + bounds)
│   ├── simulation.py        # Discrete-event FJ queue simulator
│   └── plotting.py          # Visualization: approximation vs simulation vs bounds
├── examples/
│   └── demo.py              # Runnable example showing all capabilities
├── pyproject.toml            # Package metadata + dependencies
└── (existing docs/, CLAUDE.md)
```

## Files to Create

### 1. `pyproject.toml`
- Minimal config: name `forkjoin`, Python >=3.10
- Dependencies: `numpy`, `matplotlib`
- No build system complexity — just a simple research package

### 2. `forkjoin/__init__.py`
- Re-export the key public functions from analytical and simulation modules

### 3. `forkjoin/analytical.py`
Core formulas, all taking `(lam, mu1, mu2)` as arguments and returning floats:
- `mean_response_time(lam, mu1, mu2)` — the main approximation (Section 4.2 of the doc)
- `upper_bound_independent(lam, mu1, mu2)` — T_UB
- `lower_bound_bottleneck(lam, mu1, mu2)` — T_bot
- `upper_bound_split_merge(lam, mu1, mu2)` — T_SM (P-K formula)
- `nelson_tantawi(lam, mu)` — homogeneous exact result for reference
- Input validation: check stability condition, raise `ValueError` if violated

### 4. `forkjoin/simulation.py`
- `simulate(lam, mu1, mu2, n_jobs=1_000_000, warmup=100_000, seed=None)` — discrete-event simulation
  - Returns a results dataclass with: `mean_response_time`, `std_response_time`, `n_samples`, `ci_95` (confidence interval)
  - Uses simple event-driven logic: track `last_departure` per server, compute response per job as `max(depart1, depart2) - arrival`
  - Uses `numpy.random.Generator` for performance

### 5. `forkjoin/plotting.py`
Two main plot functions:
- `plot_vs_load(mu1, mu2, lam_points=50, run_simulation=True, n_jobs=500_000)` — plots T_approx, T_UB, T_bot, (optionally T_sim) vs lambda/min(mu1,mu2) from 0.05 to 0.95
- `plot_vs_heterogeneity(mu1, lam, mu2_range, run_simulation=True, n_jobs=500_000)` — plots all quantities vs mu2/mu1 ratio
- Both return `(fig, ax)` for further customization
- Clean academic style: labeled axes, legend, gridlines

### 6. `examples/demo.py`
- Runnable script that:
  1. Prints a table comparing approximation vs simulation for several parameter sets (like Section 4.5 of the doc)
  2. Generates two plots: T vs load, and T vs heterogeneity ratio
  3. Saves plots to `examples/` as PNGs

## Verification

1. `python -c "from forkjoin import mean_response_time; print(mean_response_time(0.6, 1.0, 1.0))"` should output ~3.5625
2. `python examples/demo.py` should produce the comparison table and two plot files
3. Homogeneous case: `mean_response_time(lam, mu, mu)` must equal `nelson_tantawi(lam, mu)` for all valid inputs
