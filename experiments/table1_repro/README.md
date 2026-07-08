# Reproducing Table 1 (Performance2026 paper)

Reproduces **Table 1** (`\label{tab:approx-results}`, p.17) of the Squillante &
Tantawi *Performance2026* submission: the expected sojourn time of the **core
parallel computational phase** (a heterogeneous two-queue fork-join) from the
three interpolation approximations vs. discrete-event simulation, over
`rho ∈ {0.4, 0.8, 0.9, 0.95}` and `r ∈ {2, 4, 8}`.

## Notation mapping (paper ⇄ this repo)

| Paper symbol | Paper eq. | `docs/paper` eq. | Repo function | Nickname |
|---|---|---|---|---|
| `E[T_FJ]`      | 6  | 8  | `mean_response_time`             | UL  |
| `E[T_FJ^(0)]`  | 13 | 16 | `mean_response_time_lh` (β=10)   | LH  |
| `E[T_FJ^(1)]`  | 14 | 23 | `mean_response_time_lh_enhanced` (β=10) | LHe |

## Parameterization

With `mu^* = mu_1` the slower/bottleneck server:

```
mu_1 = 1,   mu_2 = r,   gamma = lambda = rho * mu_1
```

so `rho = gamma / mu^*` is the bottleneck utilization and `r = mu_2 / mu_1 > 1`.
Verified directly against the paper: e.g. `rho=0.4, r=2` gives `E[T_FJ] = 1.824`,
matching the paper's Table 1 cell exactly.

## Simulation protocol (from the paper)

> "The simulation of each instance … consisted of 20-million jobs with a warm-up
> period of 500-thousand jobs repeated over five independent replicas, where the
> 95% confidence interval … is computed via the t-distribution across the
> different replicas."

So: **20 M jobs, 500 K warmup, 5 independent seeds**; `T_sim` is the grand mean
over seeds and the 95% CI half-width is `t_{0.975, k-1} · s / sqrt(k)` with `s`
the sample std (ddof=1) of the per-seed means (independent-replications method).
This drives the library `forkjoin.simulate()` once per seed — the simulator is
not reimplemented here.

## Running

```bash
python reproduce_table1.py           # correct protocol: 20M jobs × 5 seeds, t-CI (~12 min)
python reproduce_table1.py --paper-exact   # reproduce PUBLISHED numbers: seed 42, 10M, normal CI (~8 min)
python reproduce_table1.py --quick   # 200K jobs × 5 seeds (smoke test, seconds)
```

Results cache to `table1_results.json` (keyed by rho, r, n_jobs, warmup, seeds),
so re-runs are instant. Outputs:

- `table1_results.json` — full numbers + per-seed means
- `table1.tex` — correct-protocol LaTeX table (5-seed t-CI)
- `table1_paperexact.tex` — exact reproduction of the published numbers
- `table1_errors.{png,pdf}` — relative-error bar chart

An extra homogeneous row (`r=1`) is simulated as a sanity check (recovers
Nelson–Tantawi to the third decimal) but is not part of the main `{2,4,8}` table.

## Two findings worth reconciling before publication

### 1. The published Table 1 does not use the protocol the text describes

The paper text states the simulation used **20M jobs, 500K warmup, 5 independent
replicas, t-distribution CI**. But the published numbers (both `T_sim` and the
CIs) are reproduced *exactly* — to the third decimal — by a **single seed (42),
10M jobs, normal-approximation CI** (the `docs/paper/compact_table.py` protocol):

| rho, r | published `T_sim ± CI` | seed 42, 10M, normal CI |
|---|---|---|
| 0.4, 2  | 1.816 ± 0.001  | 1.816 ± 0.001  |
| 0.8, 2  | 5.109 ± 0.003  | 5.109 ± 0.003  |
| 0.9, 2  | 10.133 ± 0.006 | 10.133 ± 0.006 |
| 0.95, 2 | 20.260 ± 0.012 | 20.260 ± 0.012 |

The normal-approximation CI ignores the strong autocorrelation of consecutive
response times and is far too tight near saturation. A **correct 5-seed t-CI**
(what the text describes, produced by the default run) is much wider and honest:

| rho | published CI | correct 5-seed t-CI |
|---|---|---|
| 0.4  | ±0.001 | ±0.001 |
| 0.8  | ±0.003 | ±0.019 |
| 0.9  | ±0.006 | ±0.083 |
| 0.95 | ±0.012 | ±0.295 |

Under the correct protocol the grand-mean `T_sim` also shifts slightly (single
seed 42 runs a touch high), which nudges a couple of approximation errors above
the "< 1.2%" claim (e.g. LH at rho=0.8, r=2: +1.71% vs the published +1.16%),
though every approximation still lies within the honest simulation CI at
moderate-to-high load. **Either regenerate Table 1 with the 5-seed t-CI protocol
the text describes, or change the text to match the single-seed method used.**

### 2. `E[T_FJ^(1)]` column vs. eq. 14 as typeset

## Note: `E[T_FJ^(1)]` and eq. 14 as typeset

The paper's **Table 1 `E[T_FJ^(1)]` column matches the repo's quadratic enhanced-LH
formula** (`mean_response_time_lh_enhanced`, i.e. `docs/paper` eq. 23):

```
E[T_FJ^(1)] = (c2·rho² + c1·rho + c0) / (mu_min·(1 − rho))
```

which anchors the heavy-traffic limit `h(r)` (the `c2` term) in addition to the
0th/1st light-traffic derivatives.

**Equation 14 as typeset in Performance2026** is a purely first-order form,

```
E[T_FJ^(1)] = gamma/(mu_1 − gamma) · ( 1/mu_1 + mu_1/mu_2² − 2mu_1/(mu_1+mu_2)²
              − 2mu_1²mu_2/(mu_1+mu_2)⁴ ) + (mu_1²+mu_1mu_2+mu_2²)/(mu_1mu_2(mu_1+mu_2))
```

which **does not reproduce the table** (it drops the heavy-traffic/quadratic
structure). Examples:

| rho, r | table `T_FJ^(1)` (repo quadratic) | eq. 14 as typeset |
|---|---|---|
| 0.4, 2  | 1.825  | 1.819  |
| 0.8, 2  | 5.151  | 5.080  |
| 0.95, 8 | 20.003 | 19.795 |

The script computes the literal eq. 14 too (`t_fj1_literal_eq14` in the JSON) to
document the gap. This looks like an equation-typesetting issue in eq. 14 rather
than a problem with the table numbers — worth reconciling before publication.
