# Validating Table 1 with qsim-service

**Date:** 2026-07-31
**Status:** implemented, runs complete, all 16 cells agree (2026-08-01)

## 1. Purpose

Table 1 of the Squillante & Tantawi *Performance2026* submission compares three
interpolation approximations for the heterogeneous two-queue fork-join sojourn time
against discrete-event simulation, over `rho in {0.4, 0.8, 0.9, 0.95}` and
`r in {2, 4, 8}`. Those simulation numbers come from this repo's own simulator
(`forkjoin.simulate`, a Lindley recursion), driven by
`experiments/table1_repro/reproduce_table1.py`.

This experiment re-measures the same 16 cells with an **independent** discrete-event
simulator — [qsim-service](https://github.com/modeling-analysis/qsim-service), an HTTP/JSON
wrapper around the headless JMT engine — and asks a single question:

> Do the two independent simulators agree, cell by cell, within their confidence
> intervals?

It is a validation experiment. It does not replace Table 1's numbers and does not
re-derive the approximation error columns.

## 2. Success criterion

For each cell, agreement is **CI overlap**:

```
|T_qsim - T_simple| <= hw_qsim + hw_simple
```

The report also prints `gap / (hw_qsim + hw_simple)` so near-misses are visible rather
than collapsing into a bare yes/no. Every disagreeing cell is printed to the console,
recorded in the results JSON, and listed in a summary block. No cell is dropped.

## 3. Scope

- 12 heterogeneous cells: `rho in {0.4, 0.8, 0.9, 0.95}` x `r in {2, 4, 8}`.
- 4 homogeneous cells at `r = 1`, which serve as a hard correctness gate (§7).
- `mu1 = 1`, `mu2 = r`, `lambda = rho * mu1` — the same parameterization as
  `reproduce_table1.py`, so `rho` is bottleneck utilization and `r = mu2 / mu1`.

Out of scope: re-deriving the UL / LH / LHe error columns against `T_qsim`; any change
to `forkjoin/`; any change to the paper.

## 4. Layout

```
experiments/table1_qsim/
  README.md                    # protocol differences, CI-nature table, how to run
  qsim_fj.py                   # model builder + HTTP client for ONE fork-join run
  validate_table1_qsim.py      # grid driver, replication, comparison, reporting
  expected_model.json          # fixture pinning the emitted model JSON
  qsim_results.json            # cache, one entry per individual run
  table1_qsim_comparison.tex   # side-by-side agreement table + CI-nature table
  table1_qsim_agreement.png/pdf
```

`forkjoin/` is untouched. `pyproject.toml` gains no dependency: the client is stdlib
`urllib`, and the figure uses the already-declared matplotlib.

### 4.1 Module split

Two modules, each with one job:

- **`qsim_fj.py`** — given `(lambda, mu1, mu2, seed, stopping)`, return one `FJRun`
  record: `mean`, `ci`, `samples_analyzed`, `samples_discarded`, `success`, `completed`,
  `wall_clock_seconds` (engine-side), `round_trip_seconds` (client-side). It owns the
  model JSON, the POST, and HTTP-status-to-error mapping. It knows nothing about Table 1.
  Runnable directly (`python qsim_fj.py --rho 0.8 --r 4`) to print a single run, which is
  what makes the timing probe cheap.
- **`validate_table1_qsim.py`** — walk the grid, replicate, aggregate,
  compare against `table1_repro`'s results, emit console report / `.tex` / figure. It
  knows nothing about HTTP.

Rejected alternatives: one self-contained script (mixes the HTTP and reporting concerns
into a file too large to edit reliably, and neither half is exercisable alone);
concurrency across seeds (qsim-service's design notes JMT carries static/global state and
claims cleanliness only for *sequential* requests in a warm JVM — if the grid proves too
slow, the honest parallelism is several service instances on different ports, as a
follow-up).

### 4.2 Dependency stance

The experiment talks to qsim-service over HTTP directly and does **not** import `qopt`.
The fork-join repo therefore gains no dependency on an unreleased sibling project, and
the experiment cannot break from a `qopt` refactor.

The cost is that the model JSON is transcribed from `qopt`'s emitter. Two mitigations:

1. A comment in `qsim_fj.py` citing `qopt/station.py:ForkJoinStation.sim_node` and
   `qopt/network.py:Network.to_model_dict` as the source.
2. `--cross-check-qopt`, which asserts byte-equality against `qopt`'s emitted model when
   `qopt` happens to be importable. This is what keeps the transcription from silently
   drifting.

## 5. The model

`qsim_fj.py` emits exactly this (`rho = 0.8`, `r = 4` shown):

```json
{"model": {
   "name": "fj-rho0.80-r4",
   "classes": [{"name": "jobs", "type": "open"}],
   "nodes": [
     {"name": "src", "type": "source",
      "arrivals": {"jobs": {"distribution": {"type": "exponential", "rate": 0.8}}}},
     {"name": "fj", "type": "fork-join", "join": "all",
      "branches": [
        {"service": {"jobs": {"distribution": {"type": "exponential", "rate": 1.0}}}},
        {"service": {"jobs": {"distribution": {"type": "exponential", "rate": 4.0}}}}]},
     {"name": "snk", "type": "sink"}],
   "routing": {"jobs": [{"from": "src", "to": "fj", "probability": 1.0},
                        {"from": "fj", "to": "snk", "probability": 1.0}]}},
 "seed": 0,
 "stopping": {"disableStatisticStop": true, "maxSamples": 20000000,
              "alpha": 0.05, "maxWallClockSeconds": 1800},
 "measures": ["response-time", "system-response-time"]}
```

`response-time` on a `fork-join` node is the whole fork-to-join sojourn per
qsim-service's contract — the quantity Table 1 reports.

**Only two measures are requested, deliberately.** `utilization` and `queue-length` on a
fork-join node return *join-station* numbers with `success: true` and no warning
(qsim-service#8). Omitting `measures` is worse still: the service substitutes defaults
that include exactly those two.

## 6. Protocol

### 6.1 The CI-nature problem

The two simulators' confidence intervals are **not the same kind of object**. Confirmed
by decompiling `jmt.engine.dataAnalysis.NewDynamicDataAnalyzer` (fields `batches`,
`batchMean`, `weightBatchMean`, `batchLen`, `polyOrder`, `transientLen`, `crossesNum`,
`hMeans`, `heuristicPassed`, plus a `PolyFit` class):

| | simple simulator | JMT / qsim per-run `lower`/`upper` |
|---|---|---|
| Unit of observation | 5 whole-run means | batch means **within one run** |
| Method | Student-t on 5 values | spectral analysis: batch means + polynomial fit to the log-spectrum (Heidelberger-Welch style) |
| Autocorrelation handled by | independent replications | spectral estimate of the variance of the mean |
| Warm-up | 500 K jobs, caller-specified | `transientLen`, detected per run by a crossings heuristic |
| Count | 1 per cell | 1 per seed |

Consequence for the design: the **headline** qsim CI is a Student-t interval over
per-seed means — the same construction as the simple simulator's — and JMT's intra-run CI
is retained only as a per-run diagnostic. `reproduce_table1.py` already documents the
same distinction: its `--paper-exact` mode uses a single-run normal-approximation CI and
notes it "ignores autocorrelation; too tight near rho=1".

This table is a **deliverable**: it appears in `README.md`, in the console report footer,
and as a LaTeX `tabular` block inside `table1_qsim_comparison.tex`, so it travels with
the numbers it qualifies.

### 6.2 Measured service behaviour (probe, 2026-07-31)

A single-cell probe at `rho = 0.8, r = 4` against
`localhost/qsim-service:0.1.0` (podman, temurin 17-jre) settled the protocol
empirically. Both requests used the §5 model; only `stopping` differed.

| | convergence run | fixed-length run |
|---|---|---|
| requested | `minSamples 1e6, precision 0.005, maxSamples 1e8` | `disableStatisticStop true, maxSamples 2e6` |
| `samplesAnalyzed` | **10,880** | **2,000,000** (exact) |
| mean | 4.7082 (-6.6% vs `table1`'s 5.0424) | 5.0585 (+0.3%) |
| achieved relative half-width | 7.9% | 1.4% |
| reported `success` | `true` (though the 0.005 target was missed 16x) | `true` |
| `samplesDiscarded` | 0 | 0 |
| `wallClockSeconds` | 0.26 | 10.85 |

Findings, in decreasing order of confidence:

1. **`minSamples` is silently dropped.** Code-verified: it reaches `model/Stopping.java`
   and `http/SimulationService.java`, but `translate/JsimgWriter.java` emits only
   `maxSamples`, `maxEvents`, `maxSimulated` and `disableStatisticStop`. The API accepts
   the field and it never reaches JMT. Filed upstream.
2. **`success` and `precision` in the response are not convergence evidence.** The
   convergence run stopped at 10,880 samples with a 7.9% half-width while reporting
   `success: true` and echoing `precision: 0.005`. `alpha`/`precision` *are* emitted
   (`JsimgWriter` 557-566), so this is not simply an unwired field, but the reported flags
   cannot be relied on either way.
3. **The fixed-length path is exact and is the only trustworthy run-length control.**
   `disableStatisticStop = true` with `maxSamples = N` returned exactly `N` samples.
4. **No transient removal occurs.** `samplesDiscarded` was 0 in both runs.
   **CORRECTED 2026-08-01 — this generalisation was false.** It holds at 2 M jobs but not
   at the 20 M protocol length, where JMT's crossings heuristic fires: 75 of the 80 grid
   runs discarded a non-zero transient, median 6,490 samples and max 58,255 (0.03% and
   0.29% of a run) against the reference's caller-fixed 500 K (2.44%). The real difference
   is magnitude and control, not presence versus absence. See §6.5 as amended.
5. **Throughput ~184 K jobs/s** (2 M jobs in 10.85 s engine-side), against ~2.2 M jobs/s
   for `forkjoin.simulate` — roughly 12x slower.

Two lower-confidence observations, raised upstream as questions rather than asserted:
`alpha` is written as `1 - alpha` (`0.95` for a requested `0.05`), and `variance`/`stdDev`
come back `null`.

### 6.2a Addendum (2026-07-31, later the same day): both bugs fixed upstream

qsim-service was updated while this experiment was being implemented. Findings 1 and 2
above are now **fixed**, and one of my own conclusions was **wrong**:

- **`minSamples` now reaches JMT** — qsim-service#10, fixed in `395ecfc`.
- **The `alpha` inversion was a real bug** — filed upstream as #12, fixed in `df45cd1`.
  I tested that hypothesis and dismissed it, which was an error. My probe compared
  half-widths at `alpha` 0.05 vs 0.95 on a *fixed-length* run, where the stopping logic
  is disabled, so it could not exercise the semantics; and I discarded the 1.2x ratio I
  measured because I was expecting the ~30x of a 95%-to-5% confidence flip. The real
  mechanism, per the fix's own note: JMT feeds the attribute to `TStudent.ICDF`, which
  expects a *significance* level; handed `1 - alpha` it returns the negated one-sided
  quantile, `confInt` goes negative, and `NewDynamicDataAnalyzer.HWtest` compares without
  `abs()`, so **the stopping rule passed unconditionally and every interval came back
  ~17% too narrow**. My measured 1.2x ratio *was* that signature (16.4% too narrow).
  This, not the missing floor alone, is why the convergence run stopped at 10,880 samples.

Re-probed on `qsim-service:0.2.0` (`df45cd1`), same cell and seed as §6.2:

| | 0.1.0 | 0.2.0 |
|---|---|---|
| `samplesAnalyzed` | 2,000,000 | 2,000,000 (unchanged) |
| mean | 5.058501 | 5.058501 (**bit-identical**) |
| per-run CI half-width | 0.069270 | 0.082844 (**1.196x wider**) |
| `success` | `true` | `false` |

**Consequences for this design:**

1. **Means are unaffected.** The sample path does not depend on `alpha`, so every mean
   measured under 0.1.0 remains valid. Only the per-run CI changed — and that CI is a
   diagnostic here, except in the §6.5 transient check, which must therefore run on
   0.2.0.
2. **`completed: false` is now the normal case for a fixed-length run**, not an error.
   Verified: `rho = 0.9, r = 1` returns `completed: false` with the full 200,000 samples
   analyzed and `wallClockSeconds` 1.154 against an 1800 s watchdog. §8's row treating it
   as a hard failure was correct against 0.1.0 (where the field read `true` spuriously)
   and is **wrong** against 0.2.0 — it rejected 4 of 16 cells. Run length is now gated
   solely on `samplesAnalyzed == N`, which catches a genuine watchdog kill strictly
   better, and `completed` joins `success` as a recorded diagnostic.
3. **The protocol does not change.** §6.3's fixed-length choice no longer rests on the
   bugs — it rests on its own merits: equal-length replications keep the five per-seed
   means i.i.d., which is what makes the cross-seed Student-t interval valid, and 20 M
   matches the reference protocol's job count. Convergence stopping is now *available*
   and still declined, for that reason.

### 6.3 Protocol: fixed run length at the paper's own job count

Convergence stopping is unusable per finding 2, and finding 5 shows the paper's own job
count is affordable (~2.4 h for the full grid), so the protocol is simply:

- **All 5 seeds (`0..4`) at `disableStatisticStop = true, maxSamples = 20_000_000`.**

This is simpler than the calibrate-then-fix scheme it replaces, matches the paper's job
count exactly, and dissolves the unequal-run-length problem that motivated calibration:
equal-length runs keep the 5 per-seed means i.i.d., so the cross-seed Student-t interval
is properly constructed. 5 runs per cell, not 6.

`N` is a CLI flag (`--n-jobs`, default 20 M) and part of the cache key.

### 6.4 Aggregation

Grand mean of the 5 per-seed means, half-width `t_{0.975,4} * s / sqrt(5)` with `s` the
ddof=1 sample std. The `T_975` quantile table is copied from `reproduce_table1.py` with a
comment pointing at the original — not imported, since that file is a CLI, not a library.

### 6.5 What cannot be matched: warm-up

qsim exposes **no warm-up knob**. It does, however, remove a transient — JMT detects one
per run — so the original framing of this section ("no transient removal happens at all",
inferred from a 2 M-job probe) was wrong. **Amended with measured values:** across the 80
grid runs at 20 M jobs, 75 discarded a non-zero transient, median 6,490 samples, max
58,255 — i.e. **0.03%** typical and **0.29%** worst case, against the reference's
caller-fixed 500 K = **2.44%**. The difference is magnitude and controllability, not
presence versus absence, and the direction of any residual bias is therefore not
predictable a priori. (The originally predicted downward bias is real only at short run
lengths: at 200 K jobs, where JMT discards nothing, `rho = 0.95, r = 1` came in 7.0% low.)

Since the knob does not exist, **quantify the residual bias instead of documenting it
away**: one seed per cell is additionally run at `N = 40 M`. This check is what
`--transient-check` runs.

**Two limits on how strongly its negative result can be read** — both must be stated
wherever the result is:

1. It compares against JMT's **intra-run** half-width, which §6.1 classifies as a
   diagnostic and which measures up to **2.9x wider** than the cross-seed half-width this
   experiment reports (`rho=0.95, r=1`: 0.561 vs 0.196). Shifts as large as 1.12% score
   "within CI".
2. `T(2N)` extends the **same sample path** as `T(N)` at the same seed, so if bias decays
   like `1/N` the statistic observes only about **half** of it.

So the supported claim is "no shift larger than roughly 1% at high load is detectable by
this test", not "no transient bias". The stronger evidence for the latter is the
200 K → 20 M behaviour above: −7.0% shrinking to +0.21% at `rho = 0.95, r = 1`.

Retired from the earlier draft: surfacing `samplesDiscarded` in the README as evidence of
*run length*. It does carry information about transient removal, just not about how long a
run was — `samplesAnalyzed` is the run-length field.

## 7. Correctness gates

Run in this order; each gates the next.

1. **Health preflight** — one `GET /health` before any cell, so a misconfigured URL fails
   immediately rather than on cell 1.
2. **Determinism** — the same seed twice must produce identical means. One extra run.
3. **`r = 1` closed form** — for `r = 1` the exact result is known (Nelson & Tantawi 1988):
   `T = (12 - rho) / (8 (mu - lambda))`. The four homogeneous cells are therefore an
   **analytic oracle for qsim itself**, not a sanity check. If qsim cannot reproduce
   `(12 - rho)/8 * T_1` inside its CI, nothing downstream is meaningful: the script says
   so and stops.
4. **`system-response-time` oracle** — this network has one station, so system response
   time must equal the fork-join `response-time`. The guard flags a relative difference
   above 2%; its job is catching *wrong semantics* (a join-residence number instead of a
   region number), not measuring precision. The probe (§6.2) returned the two measures
   **identical to the last digit** in both runs, so the loose threshold is slack, not a
   concession — and the fork-join region semantics the README documents are confirmed.

Only then are the 12 heterogeneous cells run.

## 8. Error handling

| Condition | Response |
|---|---|
| Service unreachable | Abort at preflight with both documented start commands. |
| HTTP 400 / 422 | Our model JSON is wrong, so every cell fails identically: dump the offending request body, abort. |
| HTTP 500 | Engine failure on one cell: record it, continue remaining cells, list all failures in the summary, exit non-zero. A two-hour grid must not die on one flaky cell. |
| `samplesAnalyzed != N` | The run did not execute the requested length, so it is not comparable to its 4 siblings: discard it, record the shortfall, and fail the cell. This is the **only** trustworthy run-length check, since `success` is not one (§6.2 finding 2). |
| `completed: false` | **Superseded by §6.2a — do not implement this row.** Recorded as a diagnostic, never acted on: on 0.2.0 it means "CI targets unmet", which is the expected state under `disableStatisticStop`, and gating on it rejected 4 of 16 valid cells. A genuine watchdog kill is caught by the `samplesAnalyzed` row above. |
| `success` / `precision` on any measure | **Recorded as a diagnostic, never acted on.** Finding 2 of §6.2: a run that missed its target 16x still reported `success: true`. Neither value gates anything; `qopt`'s `measures.extract` warns on `success: false`, and this experiment deliberately does not. |

**Exit codes distinguish a finding from a failure.** Cells that *disagree* are a
scientific result, not an error: the script reports them prominently and exits **0**. A
non-zero exit means the run itself is untrustworthy — a failed correctness gate (§7), an
HTTP 500 on any cell, or an unreachable service.

**Service lifecycle is not managed by the script.** qsim-service is a GPL Java process
with its own lifecycle; starting or killing it as a side effect of running an experiment
would be surprising. The script only ever reads `QSIM_URL` (default
`http://localhost:8080`) and fails fast with instructions.

**Cache and resumability.** `qsim_results.json` is written after **every individual run**,
not per cell, via tmp-file + `os.replace`. A `Ctrl-C` two hours in loses at most one run,
and re-invoking resumes. Cache key: `(rho, r, protocol params, seed)`.

## 9. Runtime

qsim-service runs the authoritative grid in a **podman container** on
`eclipse-temurin:17-jre` — the JVM the project targets and tests. Rationale: the local
JVM is Temurin 26, the built jar is Java-17 bytecode
(`major version: 61`, `maven.compiler.release=17`) so it would *load*, but
qsim-service's own design spec still lists "JMT 1.4.0 JVM compatibility (target 17)" as
an open item. For an artifact whose purpose is trustworthy numbers, running on the tested
JVM removes a variable that would otherwise need a caveat, and makes the run reproducible
elsewhere. Requires a one-time `podman build` (~2-3 min; Maven runs inside).

Podman on macOS runs a Linux VM, so container timing is not natively comparable to
`forkjoin.simulate` on the host. Mitigation: after the grid, re-run **2 cells** under the
local `-cp` Java process purely to quantify VM overhead, and report the timing ratio with
that caveat attached.

The JVM risk is partly self-policing: the determinism check and the `r = 1` closed-form
gate (§7) would catch a broken JMT on either JVM.

**Resolved by measurement (2026-08-01).** The 2-cell probe was run, and the JVM question
this section hedges against is settled: podman/temurin-17-jre and the local JDK 26 return
**bit-identical** numbers — mean and CI agreeing to 17 significant digits on a shared
2 M-job request, and identical means on both 20 M probe cells. JMT 1.4.0 on JDK 26 is not
a source of numerical difference, so the choice of runtime here was a precaution rather
than a necessity.

The VM's timing overhead, which was the reason for the probe, measures **+2.2%** at 20 M
jobs (+3.5% at `rho=0.8, r=4`; +0.9% at `rho=0.95, r=8`), so the reported qsim-vs-Lindley
ratio -- 8.61x work-normalised -- becomes 8.42x net of virtualization. (Do NOT quote the
raw seconds-ratio column, median 8.39x; §10 explains why it is not the honest figure.) At 2 M jobs the same overhead
measures 7.0%, because fixed startup cost weighs more at short run lengths — timings must
be read at the measurement scale.

**Measured cost (probe, §6.2):** ~184 K jobs/s engine-side, so 20 M jobs is ~109 s per
seed, ~9 min per cell, and **~2.4 h for all 16 cells** at 5 seeds. The `--transient-check`
runs (§6.5, one 40 M run per cell) add ~1.2 h. Both affordable; no reduction in `N`,
seeds, or grid size is needed.

Runs executed against `qsim-service:0.2.0` (`df45cd1`) -- the build that carries both
upstream fixes of §6.2a. Earlier builds must not be used: their per-run intervals are
~17% too narrow, and `--transient-check` gates on exactly that interval.

## 10. Timing comparison

Three numbers per qsim run, all cached: engine-side `wallClockSeconds`, client-side
round-trip seconds (their difference exposes HTTP/JSON overhead), and
`samplesAnalyzed` / `samplesDiscarded`, giving seconds per million samples.

The baseline is **re-timed fresh** under `--time-baseline`. The existing
`experiments/table1_repro/run_full.log` is *not* a usable baseline: most cells took ~46 s
for 5 seeds x 20 M jobs, but three cells logged 958 s, 1013 s and 489 s for identical
work — machine contention, not model difficulty.

Two rules make the numbers mean anything:

- **Never concurrent.** Concurrency is what produced the 958 s outlier.
- **Interleaved per cell** — simple, then qsim, cell by cell — so slow machine drift hits
  both arms equally instead of penalizing whichever ran second.
  **NOT ACHIEVED in the committed run:** `--time-baseline` ran as a second invocation
  over a warm cache, so the qsim arm was replayed from runs measured ~2.4 h earlier in
  another process. The code now detects and reports this; forcing it needs a cold cache.

**CORRECTED:** this is NOT equal work. `forkjoin.simulate` computes `total = warmup +
n_jobs`, so its 500 K warm-up is simulated IN ADDITION — 20.5 M against qsim's 20 M,
2.44% apart. Throughput is normalised by each arm's true total; the raw seconds ratio
is not the honest figure. Original text: "this **is** equal work… a genuine like-for-like
throughput comparison".
seed — so the timing ratio is a genuine like-for-like throughput comparison, which the
earlier convergence-stopping draft could not have claimed. The measured starting point is
~184 K jobs/s for qsim against ~2.2 M jobs/s for `forkjoin.simulate` (~12x), to be
confirmed across the grid.

Two caveats stay attached to the ratio: qsim's figure includes JMT's per-event machinery
and the podman VM (quantified by the 2-cell local-JVM probe), and neither simulator
removes the same transient (§6.5). Achieved CI half-width is printed beside elapsed time
in both columns so cost and accuracy are never read apart.

## 11. Testing

- `expected_model.json` pins the emitted model JSON.
- `--cross-check-qopt` asserts byte-equality with `qopt`'s emitter when importable.
- `--quick` runs the whole pipeline at `--n-jobs 200_000` (still fixed-length) for a
  seconds-long smoke test. It does **not** switch to convergence stopping, which §6.2
  finding 2 rules out entirely.
- Everything else is the correctness gate in §7 — these are real assertions against a
  closed form, not self-consistency checks.

## 12. Deliverables

1. `experiments/table1_qsim/` as laid out in §4.
2. A 16-row agreement table (`rho`, `r`, `T_simple +/- CI`, `T_qsim +/- CI`, gap,
   gap/combined-half-width, agree) in console, LaTeX, and JSON form.
2a. `table1_qsim_agreement.png/pdf`: one bar per cell, grouped by `rho`, height =
   `gap / (hw_qsim + hw_simple)`, with a shaded band at `+/-1` marking the agreement
   threshold. A bar inside the band is an agreeing cell, so the whole grid's verdict is
   readable at a glance and near-misses are visible.
3. The CI-nature table of §6.1 in the README, the console footer, and the LaTeX output.
4. A timing table per §10 (NOT equal work: 20.5 M vs 20 M — see §10), with the podman-VM caveat and the
   2-cell local probe.
5. The transient-bias result of §6.5: 20 M vs 40 M means per cell, stated as either
   "below the noise floor" or a measured bias with its sign.
6. An explicit list of disagreeing cells, or an explicit statement that there are none.
7. The §6.2 findings carried into `README.md`, so a reader knows the run deliberately
   avoids qsim's convergence stopping and why.
