# AGENTS.md

This file provides guidance to agents when working with code in this repository.

## Project Overview

This is a **research software project** implementing closed-form approximations for the mean response time of heterogeneous 2-queue fork-join systems. The project provides:

1. **Analytical formulas** for fork-join queue response time approximations
2. **Discrete-event simulation** for validation
3. **Visualization tools** for comparing approximations against simulation and theoretical bounds

### Mathematical Context

The core problem: Given a fork-join system with 2 parallel servers (service rates μ₁, μ₂) and Poisson arrivals (rate λ), approximate the mean job response time T = E[max(R₁, R₂)].

**Key Result**: The main approximation uses a convex combination interpolation:
```
T ≈ (1 - α)·T_UB + α·T_bot
where α = (ρ₁ + ρ₂)/8
```

This formula is:
- **Exact** for the homogeneous case (μ₁ = μ₂)
- Within **2% error** for heterogeneous cases (validated via simulation)
- Bounded by T_bot ≤ T ≤ T_UB

### Technology Stack

- **Language**: Python 3.10+
- **Core Dependencies**: NumPy (numerical computation), Matplotlib (visualization)
- **Package Manager**: pip (via pyproject.toml)

## Building and Running

### Installation

```bash
# Install in development mode
pip install -e .

# Or install dependencies only
pip install numpy matplotlib
```

### Running Examples

```bash
# Run the main demonstration
python examples/demo.py

# This will:
# 1. Print comparison table (approximation vs simulation)
# 2. Generate plots: response_time_vs_load.png, response_time_vs_heterogeneity.png
```

### Using the Package

```python
from forkjoin import mean_response_time, simulate

# Compute approximation
lam, mu1, mu2 = 0.6, 1.0, 1.5
t_approx = mean_response_time(lam, mu1, mu2)

# Validate with simulation
result = simulate(lam, mu1, mu2, n_jobs=1_000_000)
print(f"Approximation: {t_approx:.4f}")
print(f"Simulation: {result.mean_response_time:.4f}")
print(f"95% CI: {result.ci_95}")
```

## Code Organization

```
forkjoin/
├── __init__.py       # Public API exports
├── analytical.py     # Closed-form formulas (approximation + bounds)
├── simulation.py     # Discrete-event simulator
└── plotting.py       # Visualization functions

examples/
└── demo.py           # Demonstration script

docs/
├── implementation-plan.md  # Development roadmap
└── paper/                  # LaTeX research paper
```

### Module Responsibilities

**`analytical.py`**: All closed-form formulas
- `mean_response_time(lam, mu1, mu2)` - Main approximation
- `mean_response_time_lh(lam, mu1, mu2)` - Light-heavy traffic interpolation
- `upper_bound_independent(lam, mu1, mu2)` - T_UB (independence assumption)
- `lower_bound_bottleneck(lam, mu1, mu2)` - T_bot (bottleneck bound)
- `upper_bound_split_merge(lam, mu1, mu2)` - T_SM (Pollaczek-Khinchine)
- `nelson_tantawi(lam, mu)` - Exact homogeneous case result

**`simulation.py`**: Discrete-event simulation
- `simulate(lam, mu1, mu2, n_jobs, warmup, seed)` - Returns `SimResult` dataclass
- Uses vectorized NumPy operations for performance
- Tracks per-server departure times, computes max(depart1, depart2) - arrival

**`plotting.py`**: Visualization
- `plot_vs_load(mu1, mu2, ...)` - Response time vs load (λ/min(μ₁,μ₂))
- `plot_vs_heterogeneity(mu1, lam, mu2_range, ...)` - Response time vs μ₂/μ₁ ratio
- Both functions optionally run simulation for validation points

## Development Conventions

### Input Validation

All analytical functions validate:
1. **Positivity**: λ, μ₁, μ₂ > 0
2. **Stability**: λ < min(μ₁, μ₂)

Violations raise `ValueError` with descriptive messages.

### Numerical Precision

- Use float64 (NumPy default) for all computations
- Simulation uses large sample sizes (default: 1M jobs, 100K warmup) for statistical accuracy
- Confidence intervals computed as mean ± 1.96·(std/√n)

### Code Style

- **Docstrings**: All public functions have concise docstrings explaining parameters and return values
- **Type hints**: Used where beneficial (e.g., `SimResult` dataclass)
- **Naming**: Mathematical notation preserved where clear (e.g., `lam` for λ, `mu1`/`mu2` for μ₁/μ₂)
- **Imports**: Standard library → third-party → local, grouped and sorted

### Testing Strategy

Validation is primarily **simulation-based**:
1. **Homogeneous case**: Verify `mean_response_time(lam, mu, mu) ≈ nelson_tantawi(lam, mu)` (should be exact)
2. **Bound compliance**: Verify T_bot ≤ T_approx ≤ T_UB for all parameter sets
3. **Simulation comparison**: Run `examples/demo.py` to generate comparison tables

### Performance Considerations

- **Simulation**: Vectorized NumPy operations (avoid Python loops)
- **Plotting**: Simulation points are subsampled (default: 10 points) to balance accuracy vs runtime
- **Typical runtime**: `simulate()` with 1M jobs takes ~1-2 seconds on modern hardware

## Research Context

### Key References

The approximation builds on:
- **Nelson & Tantawi (1988)**: Exact result for homogeneous 2-queue fork-join
- **Flatto & Hahn (1984, 1985)**: Exact generating function (no closed form for mean)
- **Baccelli et al. (1989)**: Stochastic ordering (correlation reduces max)
- **Varma & Makowski (1994)**: Light/heavy traffic interpolation for symmetric systems

### Validation Approach

The approximation was validated against discrete-event simulation across:
- **Heterogeneity ratios**: μ₂/μ₁ ∈ [1, 5]
- **Load levels**: ρ₁ ∈ [0.3, 0.9]
- **Sample size**: 2M jobs per scenario (100K warmup)

Results: Errors consistently < 2%, with largest errors (~1.8%) at moderate heterogeneity (μ₂/μ₁ ≈ 1.5) under heavy load (ρ₁ = 0.9).

## Common Tasks

### Adding a New Approximation Formula

1. Add function to `analytical.py` with signature `(lam, mu1, mu2) -> float`
2. Include input validation via `_validate(lam, mu1, mu2)`
3. Export from `__init__.py`
4. Add to comparison in `examples/demo.py`
5. Update plots in `plotting.py` if desired

### Modifying Simulation Parameters

Default parameters in `simulation.py`:
- `n_jobs=1_000_000` - Number of jobs to simulate (after warmup)
- `warmup=100_000` - Warmup jobs (discarded from statistics)
- `seed=None` - Random seed (None = non-reproducible)

Increase `n_jobs` for tighter confidence intervals; decrease for faster iteration.

### Generating Publication-Quality Plots

```python
from forkjoin.plotting import plot_vs_load

fig, ax = plot_vs_load(1.0, 2.0, lam_points=100, run_simulation=True, n_jobs=2_000_000)
ax.set_ylim([0, 10])  # Customize as needed
fig.savefig('figure.pdf', dpi=300, bbox_inches='tight')
```

### Verifying Correctness

Run the demo and check:
1. **Homogeneous case** (μ₁ = μ₂): Error should be < 0.5%
2. **Bound compliance**: T_bot < T_approx < T_UB in all rows
3. **Visual inspection**: Plots should show approximation tracking simulation closely

## Notes for AI Assistants

- **Mathematical notation**: Use LaTeX in docstrings/comments when helpful (e.g., `μ₁`, `ρ = λ/μ`)
- **Stability conditions**: Always check λ < min(μ₁, μ₂) before computation
- **Simulation runtime**: Be mindful of `n_jobs` parameter; 1M jobs is reasonable for development
- **Research focus**: This is a research prototype, not production software—prioritize correctness and clarity over optimization
- **Documentation**: The README.md contains the full mathematical derivation and literature review
