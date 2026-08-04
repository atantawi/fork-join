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

So: **20 M jobs, 500 K warmup**, independent seeds; `T_sim` is the grand mean
over seeds and the 95% CI half-width is `t_{0.975, k-1} · s / sqrt(k)` with `s`
the sample std (ddof=1) of the per-seed means (independent-replications method).
This drives the library `forkjoin.simulate()` once per seed — the simulator is
not reimplemented here.

**The default is now `k=10` replicas, not the 5 the paper text states**, to
tighten the interval. Two effects compound: the *t* quantile drops
(`t_{.975,4}=2.776` → `t_{.975,9}=2.262`) and `sqrt(k)` grows, so at fixed
replica spread the half-width shrinks ~42%; in practice it roughly halved (see
the CI table below). **The paper text must therefore say ten replicas** — the
generated caption reports the actual `k`.

## Running

```bash
python reproduce_table1.py           # 20M jobs × 10 seeds, t-CI (~25 min cold, ~12 min from the 5-seed cache)
python reproduce_table1.py --paper-exact   # reproduce PUBLISHED numbers: seed 42, 10M, normal CI (~8 min)
python reproduce_table1.py --quick   # 200K jobs × 10 seeds (smoke test, seconds)
python reproduce_table1.py --seeds 0 1 2 3 4   # the earlier 5-replica table
```

Results cache to `table1_results.json`, **one entry per replication** (keyed by
rho, r, n_jobs, warmup, *seed*), so raising the replica count re-runs only the
new seeds and re-runs are otherwise instant. `simulate()` is deterministic in
the seed, so a cached replica mean is identical to re-running it. (Legacy
per-seed-*list* cache entries are migrated to per-seed keys automatically on
load.) Outputs:

- `table1_results.json` — per-replication cache (mean + normal-CI half-width)
- `table1.tex` — the table (10-seed t-CI)
- `table1_paperexact.tex` — exact reproduction of the published numbers
- `table1_errors.{png,pdf}` — relative-error bar chart
- `eq23_vs_eq14_comparison.md` — full 12-cell eq. 23 vs eq. 14 write-up (Finding 2)
- `eq23_vs_eq14.{png,pdf}` — eq. 23 vs eq. 14 comparison figure
  (`generate_eq23_vs_eq14_plot.py`)

The superseded **5-replica** artifacts are kept alongside for comparison:
`table1_5rep.tex`, `table1_errors_5rep.{png,pdf}`, `table1_preview_5rep.pdf`,
`eq23_vs_eq14_5rep.{png,pdf}`. (The console logs `run_full.log` /
`run_full_10rep.log` are local only — `*.log` is gitignored.) Only the
approximation columns are identical between the two — every `T_sim`, CI, and
error differs, and the bolding differs in one cell.

An extra homogeneous row (`r=1`) is simulated as a sanity check but is not part
of the main `{2,4,8}` table. Over 10 replicas it agrees with the exact
Nelson–Tantawi `T_2 = (12-rho)/(8(mu-lambda))` **within the 95% CI at all four
loads** — exactly to the third decimal at `rho=0.4` (2.4166 vs 2.4167), but only
to ~0.3% at `rho=0.95` (27.715 ± 0.188 vs 27.625), since the CI itself is ±0.19
there:

| rho | `T_sim(r=1)` ± CI | Nelson–Tantawi exact | err |
|---|---|---|---|
| 0.4  | 2.4166 ± 0.0006  | 2.4167  | −0.00% |
| 0.8  | 7.0020 ± 0.0076  | 7.0000  | +0.03% |
| 0.9  | 13.8881 ± 0.0408 | 13.8750 | +0.09% |
| 0.95 | 27.7154 ± 0.1875 | 27.6250 | +0.33% |

## Rendering a table to PDF

`table1.tex` / `table1_paperexact.tex` are table *fragments* (they use
`\toprule`, `\multirow`, the `\ex[...]` operator, and `\eqref{...}` refs to the
paper), so they need a small wrapper before they compile. `table1_preview.tex`
is that wrapper (defines `\ex`, drops the paper cross-refs, loads
`booktabs`/`multirow`):

```bash
latexmk -pdf table1_preview.tex        # -> table1_preview.pdf
```

`latexmk` runs `pdflatex` the required 2× (for `\multirow`). To preview the
published-numbers table instead, edit the `\input{...}` line in
`table1_preview.tex` to point at `table1_paperexact.tex`.

## The final Table 1 (updates A + B applied)

**`table1.tex` is the final table.** It applies both corrections:

- **A — multi-replica simulation.** `T_sim` is the grand mean over 10 independent
  replications of 20 M jobs (500 K warm-up each) and the CI is the Student-$t$
  interval across the 10 replica means (the paper text says five — it needs
  updating to ten). The caption
  now states this explicitly, because the CI no longer means what it did in the
  published table (a within-run normal approximation).
- **B — correct `E[T_FJ^(1)]`.** The column is the quadratic enhanced-LH form
  (`mean_response_time_lh_enhanced`, `docs/paper` eq. 23), which retains the
  heavy-traffic anchoring term `c_2 rho^2`. This is what the published column
  already contained, so **the table numbers were right and eq. 14 is what needs
  fixing**: it must be typeset *with* the `c_2 rho^2` term. A LaTeX comment at
  the top of `table1.tex` records this.

All 12 cells were verified by recomputing `T_sim`, the CI, the three
approximations, every error, and the bold placement directly from the
per-replication cache. `table1_preview.pdf` is the rendered result.

### Two consequences for the paper text

1. **The "< 1.2%" accuracy claim no longer holds.** Under the honest protocol the
   largest error is **+1.73%** (`E[T_FJ^(0)]` at `rho=0.8, r=2`); `E[T_FJ^(1)]`
   peaks at +1.38% (same cell) and `E[T_FJ]` at +0.57%. In the published
   single-replica table the max was 1.16%, which is what made "< 1.2%" true. The
   claim needs to be relaxed (e.g. "< 2%") or restated per-approximation — mean
   |error| is **0.17% / 0.61% / 0.36%** for `E[T_FJ]` / `E[T_FJ^(0)]` /
   `E[T_FJ^(1)]`.
2. **Which approximation "wins" changes.** Bold counts shift from
   `E[T_FJ]` 3 / `E[T_FJ^(0)]` 6 / `E[T_FJ^(1)]` 3 in the published table to
   **`E[T_FJ]` 9 / `E[T_FJ^(0)]` 2 / `E[T_FJ^(1)]` 1**. Under the multi-replica
   reference, the *simplest* interpolation (eq. 6) is most accurate in 9 of 12
   cells. Any text asserting the enhanced approximations are the most accurate
   needs revisiting.

   Caveat worth stating in the paper: in **4 of 12 cells** (`rho=0.9, r=8` and
   all of `rho=0.95`) the spread among the three approximations is *smaller than
   the CI half-width* and all three lie inside the simulation CI — so at high
   load the bolding ranks differences the simulation cannot resolve. (With 5
   replicas this was 6 of 12; the tighter 10-replica CI resolves two more cells,
   `rho=0.8, r=8` and `rho=0.9, r=4`.)

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
response times and is far too tight near saturation. A **correct multi-seed
t-CI** is much wider and honest — and going from 5 to 10 replicas roughly halves
it:

| rho | published CI | 5-seed t-CI | **10-seed t-CI (current)** |
|---|---|---|---|
| 0.4  | ±0.001 | ±0.001 | **±0.001** |
| 0.8  | ±0.003 | ±0.019 | **±0.008** |
| 0.9  | ±0.006 | ±0.083 | **±0.039** |
| 0.95 | ±0.012 | ±0.295 | **±0.149** |

The 10-seed reduction (~50–57%) beats the ~42% that the `t`/`sqrt(k)` change
alone predicts, because the five added replicas also happened to lower the
replica-to-replica sample std. The interval is still an order of magnitude wider
than the published normal-approximation one at `rho=0.95` — that gap is the
autocorrelation the normal CI ignores, and no number of replicas removes it.

Under the correct protocol the grand-mean `T_sim` also shifts slightly (single
seed 42 runs a touch high), which nudges a couple of approximation errors above
the "< 1.2%" claim (e.g. LH at rho=0.8, r=2: +1.73% vs the published +1.16%),
though every approximation still lies within the honest simulation CI at
`rho=0.95`. **Either regenerate Table 1 with the multi-seed t-CI protocol the
text describes, or change the text to match the single-seed method used.**

> **Resolved:** regenerated with the multi-seed t-CI protocol, now at 10 replicas
> — this is update A in
> [The final Table 1](#the-final-table-1-updates-a--b-applied) above.

### 2. `E[T_FJ^(1)]` column vs. eq. 14 as typeset

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

This is **exactly eq. 23 with `c2` set to 0** (verified algebraically), so it
drops the heavy-traffic/quadratic anchoring and **does not reproduce the table**.

**Full-grid comparison (all 12 Table 1 cells), not just the 3 originally
flagged:** eq. 23 is uniformly larger than eq. 14 — by the dropped term
`c2·rho²/(mu_min·(1−rho))` — with the relative gap ranging **+0.18% to +2.81%**
(mean +1.35%). The gap grows toward saturation and is *non-monotonic in `r`*,
peaking at `r=4` (because `c2` peaks there and → 0 as `r → ∞`). Against
simulation, the quadratic eq. 23 has mean |error| **0.36%** (more accurate in
**9/12** cells) vs **1.07%** for eq. 14, which under-predicts systematically at
moderate-to-heavy load (worst −2.88% at `rho=0.95, r=4`).

| rho, r | table `T_FJ^(1)` (repo quadratic, eq. 23) | eq. 14 as typeset (c2=0) |
|---|---|---|
| 0.4, 2  | 1.825  | 1.819  |
| 0.8, 2  | 5.151  | 5.080  |
| 0.95, 8 | 20.003 | 19.795 |

The script computes the literal eq. 14 too (`t_fj1_literal_eq14` in the JSON) to
document the gap. This looks like an equation-typesetting issue in eq. 14 rather
than a problem with the table numbers — worth reconciling before publication.

> **Resolved:** the table keeps the quadratic eq. 23 column (update B); the fix
> belongs in the paper, where eq. 14 must be typeset *with* the `c_2 rho^2` term.
> This is a `docs/paper`/manuscript edit and is **not yet applied**.

See **`eq23_vs_eq14_comparison.md`** (+ figure `eq23_vs_eq14.png`, script
`generate_eq23_vs_eq14_plot.py`) for the full 12-cell table, the algebraic
`c2=0` equivalence, and the accuracy analysis.
