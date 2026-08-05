# Validating Table 1 with qsim-service

Re-measures the 16 simulated cells of **Table 1** (`\label{tab:approx-results}`) with an
**independent** discrete-event simulator — [qsim-service](https://github.com/modeling-analysis/qsim-service),
an HTTP/JSON wrapper around the headless JMT engine — and reports whether the two
simulators agree, cell by cell, within their confidence intervals.

This does **not** replace Table 1's numbers and does not re-derive the approximation
error columns. Design: `docs/superpowers/specs/2026-07-31-qsim-table1-validation-design.md`.

> **Result: all 16 cells agree.** The full grid ran on qsim-service `0.2.0` — 16 cells x
> 5 seeds x 20 M fixed-length jobs, 8,548 engine-seconds (2.37 h). Every qsim mean falls
> inside the combined 95% CI of the reference simulator's; largest relative gap **0.342%**
> (`rho=0.95, r=4`), median **0.091%**, largest gap/combined-half-width **0.56**. No cell
> disagrees. All result artifacts are committed; see
> [Artifact provenance](#artifact-provenance) for exactly what each one does and does not
> reflect.

## Parameterization

`mu_1 = 1`, `mu_2 = r`, `lambda = rho * mu_1` — identical to `../table1_repro/`, so
`rho` is bottleneck utilization and `r = mu_2 / mu_1`.

Grid: `rho in {0.4, 0.8, 0.9, 0.95}`, `r in {1, 2, 4, 8}` (`r = 1` is the closed-form gate).

## Protocol, and why it is not the paper's

Per seed: a **fixed-length** run of 20 M jobs (`disableStatisticStop: true`,
`maxSamples: 20000000`) — the paper's own job count. Five seeds (`0..4`), aggregated by
the independent-replications method: grand mean plus `t_{0.975,4} * s / sqrt(5)`.

Two deliberate departures, both forced by measured service behaviour:

1. **No convergence stopping.** Two upstream bugs originally forced this and are now
   **fixed** in qsim-service 0.2.0: `minSamples` never reached JMT
   ([#10](https://github.com/modeling-analysis/qsim-service/issues/10), fixed in `395ecfc`) and the
   measure `alpha` was written as its complement, which made JMT return a negated
   one-sided quantile — so the stopping rule passed unconditionally and every interval came
   back ~17% too narrow (#12, fixed in `df45cd1`). Together they produced a run of **10,880
   samples** reporting `success: true` against a `precision: 0.005` target it missed by 16x.

   Fixed-length runs remain the choice **on their own merits**, not because of the bugs:
   equal-length replications keep the five per-seed means i.i.d., which is what makes the
   cross-seed Student-t interval valid, and 20 M matches the reference protocol's job
   count. Run length is verified by `samplesAnalyzed == N` — never by `success` or
   `completed`, both of which report "CI targets unmet" under `disableStatisticStop` and
   are recorded as diagnostics. On 0.2.0, `completed: false` is the *normal* result for a
   high-load fixed-length cell: `rho = 0.9, r = 1` returns it with the full sample count
   and 1.154 s against an 1800 s watchdog.
2. **No warm-up *knob*.** qsim has no caller-settable warm-up. It does remove a
   transient — JMT detects one per run — but the amount is neither controllable nor
   comparable to the reference's. Measured over the 80 grid runs: **75 discarded a
   non-zero transient**, median 6,490 samples and max 58,255, i.e. **0.03% and 0.29% of a
   run** against the reference's caller-fixed 500 K = **2.44%**. So the difference is one
   of magnitude and control, not presence versus absence. Residual bias is therefore
   measured, not assumed: `--transient-check` compares N against 2N per cell.

   *An earlier version of this file claimed `samplesDiscarded` is always 0 and that qsim
   performs no transient removal. That was generalised from a 2 M-job probe — where it is
   0 — and is false at the 20 M protocol length. Corrected here and in the generated
   LaTeX.*

Run length matters more at short lengths than the above might suggest. Measured against
the **exact closed form** at both lengths, so the baseline does not shift: at 200 K
jobs/seed JMT discards nothing and `rho = 0.95, r = 1` came back **25.70** vs exact
**27.63** = **−6.95%**, biased low as an un-removed empty-system transient predicts; at
20 M the same cell is **27.81 = +0.69%**. A ~10x reduction, and the sign flips. This is
one reason `--quick` is a pipeline check and not evidence of agreement.

(Against the *reference simulator* rather than the exact form, the 20 M figure is +0.21% —
a different and smaller number. Quoting the two baselines interchangeably understates
residual deviation by 3.2x, which an earlier version of this file did.)

## The two CIs are not the same kind of object

Both columns in the comparison table are cross-seed Student-t intervals, built by one
code path, so they are comparable. qsim *also* reports a per-run interval, and that one
is a different animal — recorded as a diagnostic only, never as the reported CI:

|  | Simple simulator | qsim per-run CI |
|---|---|---|
| Unit of observation | 5 whole-run means | batch means **within one run** |
| Method | Student-t on 5 values | spectral: batch means + polynomial fit to the log-spectrum |
| Autocorrelation | independent replications | spectral variance-of-the-mean estimate |
| Warm-up | 500 K jobs (2.44%), caller-set | auto-detected per run; median 6.5 K (0.03%), max 58 K (0.29%) |
| Count | 1 per cell | 1 per seed |

Determined by decompiling `jmt.engine.dataAnalysis.NewDynamicDataAnalyzer` (fields
`batches`, `batchMean`, `weightBatchMean`, `batchLen`, `polyOrder`, `transientLen`,
`crossesNum`, plus a `PolyFit` class). `../table1_repro/reproduce_table1.py` documents
the same distinction for its own `--paper-exact` mode.

## Correctness gates

Run in order, each gating the next:

1. `GET /health`.
2. **Determinism** — the same seed twice must give the same mean.
3. **The `r = 1` closed form** — for equal rates the exact result is known
   (Nelson & Tantawi 1988, `T = (12 - rho) / (8 (mu - lambda))`), making those four
   cells an analytic oracle for qsim itself. If they fail, the heterogeneous cells do
   not run. A missing `r = 1` cell also fails the gate: a gate that inspects zero cells
   is not a passing gate.

   Its tolerance is `hw_qsim + 1% * exact`, not `hw_qsim` alone, because **both**
   simulators deviate from the exact form at high load by more than pure sampling noise
   suggests. The reference sits +0.01% above exact at `rho = 0.4` rising to +0.47% at
   `rho = 0.95`, the latter consuming 39% of its own CI; qsim sits **+0.69%** above exact
   at `rho = 0.95`, which a CI-only tolerance scores at **0.968** — a 3% margin from
   failing the gate and skipping all 12 heterogeneous cells. So the floor is load-bearing,
   and the reason is finite-run behaviour as `rho → 1`, *not* (as an earlier version of
   this section claimed) that qsim removes no warm-up. Note the deviation is **upward**,
   the opposite sign to what an un-removed transient would cause. The gate exists to catch
   a *broken* fork-join implementation, wrong by tens of percent; the measured gap is
   printed either way.

Additionally, **every run** asserts `system-response-time == response-time`. This network
has one station, so they must agree; the service returned them identical to the last
digit, confirming the fork-join `response-time` carries whole-region semantics rather
than the join station's own residence time.

The hand-transcribed model JSON was also checked against `qopt`'s emitter
(`PYTHONPATH=~/Projects/quantum/optimizer python -c "import qsim_fj;
print(qsim_fj.cross_check_qopt(0.8, 1.0, 4.0, 'x'))"` → `ok: identical to qopt's
emitter`), so the duplication in §4.2 of the design is verified, not merely guarded.

## Running

Start the service (never started or stopped by these scripts):

**Requires qsim-service at `df45cd1` or later** (tagged `0.2.0` locally). Earlier builds
carry the `alpha` inversion, which makes every per-run interval ~17% too narrow — that
interval is what `--transient-check` tests against.

```bash
# podman/docker -- the authoritative runtime (temurin 17-jre, the tested JVM)
cd ~/Projects/quantum/qsim-service
podman build -t qsim-service:0.2.0 .
podman run -d --name qsim -p 8080:8080 qsim-service:0.2.0

# or a local JVM (JDK 26 verified bit-identical to temurin 17 -- see the probe below)
mvn package
java -cp "target/qsim-service.jar:target/dependency/*:lib/JMT-singlejar-1.4.0.jar" \
     qsim.http.App
```

Then:

```bash
python qsim_fj.py --rho 0.8 --r 4 --cross-check-qopt   # one run, ~11 s
python validate_table1_qsim.py --quick                 # smoke test, ~2 min
python validate_table1_qsim.py                         # full grid, 2.37 h measured
python validate_table1_qsim.py --time-baseline --transient-check   # + ~1.9 h
python -m pytest                                       # 42 unit tests, no network
```

`QSIM_URL` overrides the default `http://localhost:8080`. Every run is cached in
`qsim_results.json` immediately — flushed per run, not per cell — so an interrupted
grid resumes and loses at most one run.

**Exit codes:** disagreeing cells exit **0** — that is the experiment's finding, not its
failure. Non-zero means the run itself is untrustworthy: a failed gate, an engine error,
or an unreachable service.

## Cost, measured

Measured at the **20 M protocol length** (work-normalized medians): **simple 1,611,809
jobs/s vs qsim 187,241 jobs/s → qsim is 8.61x slower**, or **8.42x** net of the podman VM
overhead quantified below. The grid itself took 8,548 engine-seconds (2.37 h);
`--transient-check` added 3,423 s.

Read the `jobs/s` columns, not the raw `ratio` column — see the equal-work note below. (An
earlier version of this section quoted 9.8x from the 200 K smoke run, and a headline 8.4x
taken from the raw-seconds column this file elsewhere calls dishonest. Both corrected.)

**The timing comparison is not equal work, despite what the design's §10 claims.**
`forkjoin.simulate` computes `total = warmup + n_jobs`, so its 500 K warm-up jobs are
simulated *in addition* to `n_jobs` — 20.5 M against qsim's 20 M. The `jobs/s` columns
are normalized by each arm's true total and are the honest comparison; the raw
seconds-ratio column is not. One further caveat remains: the two arms do not remove the
same transient (§6.5).

### Podman VM overhead: measured, +2.2%

qsim's seconds are engine-side inside a Linux VM, so they are not natively comparable to
`forkjoin.simulate` on the host. Quantified by running two cells on both runtimes,
interleaved, at the full 20 M jobs (seed 0) — `podman_overhead_probe.json`:

| cell | podman (temurin 17-jre) | local JVM (JDK 26) | overhead |
|---|---|---|---|
| `rho=0.8, r=4` | 107.7 s | 104.0 s | +3.5% |
| `rho=0.95, r=8` | 110.8 s | 109.8 s | +0.9% |
| | | | **mean +2.2%** |

So the 8.61x work-normalized figure is essentially unaffected by virtualization;
correcting for the VM gives 8.42x. (At 2 M jobs the overhead measures 7.0%, since fixed startup cost weighs
more at short run lengths — another reason to read timings at the measurement scale.)

**This also settles the JVM-compatibility question the design left open.** qsim-service
targets Java 17 and its spec listed "JMT 1.4.0 JVM compatibility" as unresolved, which is
why the authoritative grid ran on temurin 17-jre. The two runtimes return
**bit-identical** results — mean and CI agreeing to 17 significant digits on a shared
2 M-job request, and identical means on both 20 M probe cells. JMT 1.4.0 on JDK 26 is
therefore not a source of numerical difference.

To reproduce:

```bash
cd ~/Projects/quantum/qsim-service     # must be at df45cd1 or later
mvn -DskipTests package                # only if target/qsim-service.jar predates that
QSIM_PORT=8081 java -cp "target/qsim-service.jar:target/dependency/*:lib/JMT-singlejar-1.4.0.jar" \
    qsim.http.App &
cd -
for cell in "0.8 4" "0.95 8"; do
  set -- $cell
  python qsim_fj.py --rho "$1" --r "$2" --n-jobs 20000000 --url http://localhost:8080  # podman
  python qsim_fj.py --rho "$1" --r "$2" --n-jobs 20000000 --url http://localhost:8081  # local
done
```

## Results

**Agreement.** All 16 cells agree; **no cell disagrees**. Largest relative gap 0.342%
(`rho=0.95, r=4`), median 0.091%, largest gap/combined-half-width 0.56 — no cell
approaches its CI boundary. The `r=1` closed-form gate passed at all four `rho`. Full
table in `run_full.log` and `table1_qsim_comparison.tex`.

**Transient bias.** No shift was detectable in any of the 16 cells when run length was
doubled from 20 M to 40 M jobs. Two limits on how strongly that can be read:

- The comparison uses JMT's **intra-run** half-width, which this file classifies as a
  diagnostic and which runs up to **2.9x wider** than the cross-seed half-width actually
  reported (at `rho=0.95, r=1`: 0.561 vs 0.196). The largest observed shift is **1.12%**.
- `T(2N)` extends the *same* sample path as `T(N)` at the same seed, so if bias decays
  like `1/N` this statistic sees only about **half** of it.

Supported claim: *no transient bias larger than roughly 1% at high load is detectable by a
single-seed doubling against a wide spectral interval* — not "no transient bias". The
directly relevant evidence is the run-length trend measured against the exact form at a
single baseline: `rho = 0.95, r = 1` moves from **−6.95%** at 200 K jobs to **+0.69%** at
20 M.

## Artifact provenance

Everything here is committed, but the pieces were produced at different points and one of
them predates corrections made afterwards. Precisely:

- **`qsim_results.json`** (the cache) and the **means, CIs and verdicts** derived from it
  are the run itself and are authoritative.
- **`table1_qsim_results.json`**, **`table1_qsim_comparison.tex`** and the **figure** are
  regenerated from that cache by the current code, so they carry the corrected CI-nature
  table.
- **Two fields the code now records are absent from this run's data**: `completed` per run
  and `fresh_runs` per cell were added while fixing review findings, *after* the grid ran,
  so the cache does not contain them and nothing downstream can reconstruct them. A future
  run will have both.
- **`run_full.log` and `run_timing.log` are verbatim records of the original run and are
  deliberately not edited.** They therefore contain three sentences the current code can no
  longer emit, all retracted:
  - `"none observed (samplesDiscarded = 0)"` — false at 20 M; see the warm-up row above.
  - `"TIMING (interleaved per cell …)"` — the qsim arm came from cache, so that pass was
    **not** interleaved.
  - `"All cells: the un-removed transient is below the noise floor."` — overstated; see the
    bound in [Results](#results).

  Read the logs for the numbers, which are unaffected, and this README for the claims.

## Outputs

| File | Contents |
|---|---|
| `qsim_results.json` | Per-run cache (mean, samples, timings, per-run CI) |
| `table1_qsim_results.json` | Per-cell comparison, timing, transient check |
| `table1_qsim_comparison.tex` | Agreement table + the CI-nature table above |
| `table1_qsim_agreement.{png,pdf}` | Gap / combined half-width, with a ±1 agreement band |

## Files

| File | Responsibility |
|---|---|
| `qsim_fj.py` | One fork-join run: model JSON, POST, error mapping, `FJRun`. Also a CLI for single-run probes. |
| `validate_table1_qsim.py` | Grid walk, replication, aggregation, comparison, reporting, timing. |
| `expected_model.json` | Fixture pinning the emitted model JSON. |
| `test_qsim_fj.py` | 16 client tests (fake transport, no network). |
| `test_validate.py` | 26 tests for aggregation, reference loading, closed form, cache, the r=1 gate, timing normalisation, the verdict guards, and the undecidable-vs-failed, run-length-warning, and seed-independence regressions. |
| `podman_overhead_probe.json` | The 2-cell podman-vs-local-JVM timing measurement. |

Stdlib only (`urllib`, `json`, `math`, `os`) plus matplotlib for the figure. `qopt` is
never imported at runtime — only inside the guarded `cross_check_qopt`.
