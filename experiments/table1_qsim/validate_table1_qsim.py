#!/usr/bin/env python3
"""Validate Table 1's fork-join cells against qsim-service (headless JMT).

Table 1 of the Performance2026 submission reports simulated sojourn times from this
repo's own Lindley-recursion simulator. This script re-measures the same 16 cells
with an independent discrete-event simulator and asks one question: do the two
agree, cell by cell, within their confidence intervals?

Design: docs/superpowers/specs/2026-07-31-qsim-table1-validation-design.md

The two CIs are NOT the same kind of object. The simple simulator's is a Student-t
interval over 5 whole-run means; JMT's per-run interval is a batch-means/spectral
estimate within ONE run. To compare like with like, this script builds a
cross-seed Student-t interval from qsim's per-seed means using the same
construction, and treats JMT's intra-run interval as a diagnostic only (spec 6.1).
"""

import argparse
import json
import math
import os
import sys
import time

import qsim_fj

MU1 = 1.0
RHO_VALUES = [0.4, 0.8, 0.9, 0.95]
R_VALUES = [2, 4, 8]
R_CHECK = 1                      # homogeneous: an exact-closed-form gate (spec 7.3)
SEEDS = [0, 1, 2, 3, 4]
DEFAULT_N_JOBS = 20_000_000      # the paper's own job count (spec 6.3)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REFERENCE_FILE = os.path.join(SCRIPT_DIR, "..", "table1_repro", "table1_results.json")
REFERENCE_N_JOBS = 20_000_000
"""Run length baked into the reference cache's keys, hence into every comparison.

Named rather than inlined because it is the thing a short run silently violates:
whatever `--n-jobs` says, the reference arm was measured at THIS length, so any
other length compares against a CI it cannot match. See protocol_mismatch_warning.
"""
REFERENCE_KEY = f"{{rho}},{{r}}|{REFERENCE_N_JOBS}|500000|0,1,2,3,4"

# Student-t 0.975 quantiles by degrees of freedom. Copied verbatim from
# experiments/table1_repro/reproduce_table1.py so both arms of the comparison use
# an identical interval construction (that file is a CLI, not an importable lib).
T_975 = {1: 12.706205, 2: 4.302653, 3: 3.182446, 4: 2.776445, 5: 2.570582,
         6: 2.446912, 7: 2.364624, 8: 2.306004, 9: 2.262157}


def grand_mean_and_ci(per_seed_means):
    """Grand mean and 95% t-CI half-width across independent replications."""
    k = len(per_seed_means)
    if k == 0:
        raise ValueError("need at least one per-seed mean")
    gm = sum(per_seed_means) / k
    if k == 1:
        return gm, float("nan")
    if k - 1 not in T_975:
        raise ValueError(
            f"no t quantile for {k - 1} degrees of freedom ({k} replications); T_975 "
            f"covers up to {max(T_975) + 1}. Add the quantile or use fewer seeds."
        )
    s = math.sqrt(sum((m - gm) ** 2 for m in per_seed_means) / (k - 1))
    return gm, T_975[k - 1] * s / math.sqrt(k)


def homogeneous_exact(rho, mu=MU1):
    """Exact mean response time for the HOMOGENEOUS 2-queue fork-join (r = 1).

    Nelson & Tantawi (1988): T = (12 - rho) / (8 (mu - lam)), with lam = rho * mu.
    This is what makes the r=1 cells an analytic oracle for qsim itself (spec 7.3).
    """
    return (12.0 - rho) / (8.0 * (mu - rho * mu))


def load_reference(path=None):
    """Reference means/CIs from table1_repro's cache, keyed by (rho, r).

    The CI is REBUILT here from the cached per-seed means rather than read from a
    stored half-width, so the simple simulator's interval and qsim's are produced
    by one code path and cannot differ by construction.
    """
    path = REFERENCE_FILE if path is None else path
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"reference cache not found: {path}\n"
            f"Run experiments/table1_repro/reproduce_table1.py first."
        )
    with open(path) as f:
        cache = json.load(f)
    reference = {}
    for rho in RHO_VALUES:
        for r in list(R_VALUES) + [R_CHECK]:
            key = REFERENCE_KEY.format(rho=rho, r=r)
            entry = cache.get(key)
            if entry is None:
                continue
            reference[(rho, r)] = grand_mean_and_ci(entry["per_seed_means"])
    return reference


def agreement(t_a, hw_a, t_b, hw_b):
    """(gap, gap/combined-half-width, agree) for two independent estimates.

    `agree` is None when either half-width is unavailable: agreement is then
    undecidable, which must not be reported as success.
    """
    gap = t_b - t_a
    if math.isnan(hw_a) or math.isnan(hw_b):
        return gap, float("nan"), None
    combined = hw_a + hw_b
    return gap, abs(gap) / combined, abs(gap) <= combined


def run_key(rho, r, n_jobs, seed):
    """Cache key for one run. Includes n_jobs so runs of different length never mix."""
    return f"{rho},{r}|n={n_jobs}|seed={seed}"


class RunCache:
    """Per-run JSON cache, written after every single run.

    A full grid takes hours, so the cache is flushed per run rather than per cell:
    an interrupted session loses at most one run and re-invoking resumes. Writes go
    to a temp file and are renamed, so the JSON on disk is never half-written
    (spec 8).
    """

    def __init__(self, path):
        self.path = path
        self._data = {}
        if os.path.exists(path):
            with open(path) as f:
                self._data = json.load(f)

    def get(self, key):
        return self._data.get(key)

    def put(self, key, record):
        self._data[key] = record
        tmp = f"{self.path}.tmp"
        with open(tmp, "w") as f:
            json.dump(self._data, f, indent=2)
        os.replace(tmp, self.path)


CI_NATURE_ROWS = [
    ("Unit of observation", "5 whole-run means", "batch means within ONE run"),
    ("Method", "Student-t on 5 values",
     "spectral: batch means + polyfit to log-spectrum"),
    ("Autocorrelation", "independent replications",
     "spectral variance-of-the-mean estimate"),
    ("Warm-up", "500 K jobs (2.44%), caller-set",
     "auto-detected per run; at 20 M: median 6.5 K (0.03%)"),
    ("Count", "1 per cell", "1 per seed"),
]
"""Why qsim's per-run CI is not comparable to the simple simulator's (spec 6.1).

Reported everywhere the numbers are, so a reader cannot pick up a half-width
without also picking up what produced it.

The warm-up row corrects an earlier claim that qsim performs NO transient removal.
That was generalised from a 2 M-job probe, where `samplesDiscarded` was indeed 0.
At the 20 M protocol length JMT's crossings heuristic does fire: 75 of the 80 grid
runs discarded a non-zero transient, median 6,490 samples and max 58,255 (0.03%
and 0.29% of a run). So the difference from the reference is one of MAGNITUDE and
CONTROL -- 0.03% auto-detected per run versus 2.44% fixed by the caller -- not
presence versus absence.
"""

GATE_RELATIVE_FLOOR = 0.01
"""Relative tolerance floor added to the r=1 gate's CI-based tolerance.

Measured justification: BOTH simulators deviate from the exact closed form at high
load by more than sampling noise alone would suggest. The reference -- which
discards a fixed 500 K warm-up -- sits +0.01% above exact at rho=0.4 and +0.47% at
rho=0.95, the latter consuming 39% of its own CI. qsim, measured after the fact,
sits +0.69% above exact at rho=0.95, which a CI-only tolerance scores at 0.968 --
a 3% margin from failing the gate and skipping all 12 heterogeneous cells.

So the floor is load-bearing, and the reason is finite-run behaviour at rho -> 1
generally, NOT (as an earlier version of this docstring claimed) that qsim removes
no warm-up. It does remove one: auto-detected, median 0.03% of a run against the
reference's fixed 2.44%. Note also that the observed deviation is UPWARD, the
opposite sign to the downward bias an un-removed transient would produce -- another
reason not to attribute the floor's necessity to warm-up handling.

The gate's job is catching a BROKEN fork-join implementation -- which would be
wrong by tens of percent -- not measuring precision, exactly as with the 2%
semantics oracle in spec 7.4. The measured gap is printed either way, so a real
drift stays visible even when the gate passes.
"""


def protocol_mismatch_warning(n_jobs):
    """Warning text when a run's length departs from the reference's, else None.

    Keyed on the LENGTH, not on the --quick flag. `--n-jobs 200000` takes exactly the
    same trap as --quick -- a 200 K qsim CI compared against the 20 M reference's
    roughly 10x tighter one -- but used to reach "ALL CELLS AGREE" at exit 0 with no
    warning at all, because the guard tested the flag instead of the condition it
    stood for. At 200 K a rho=0.95 cell that is 7% wrong still scores "agree".
    """
    if n_jobs == REFERENCE_N_JOBS:
        return None
    return (f"WARNING: comparing {n_jobs:,}-job qsim runs against the "
            f"{REFERENCE_N_JOBS:,}-job reference, whose CI is far tighter (at 200 K a "
            f"rho=0.95 cell that is 7% wrong still scores 'agree'). This run's verdict "
            f"is a smoke test of the pipeline, NEVER evidence that cells agree.")


def check_determinism(url, n_jobs=200_000, rho=0.8, r=4.0):
    """Same seed twice must give the same mean (spec 7.2)."""
    kwargs = dict(seed=99, n_jobs=n_jobs, url=url)
    first = qsim_fj.run_one(rho, MU1, r * MU1, **kwargs)
    second = qsim_fj.run_one(rho, MU1, r * MU1, **kwargs)
    if first.mean != second.mean:
        raise AssertionError(
            f"determinism gate FAILED: same seed gave {first.mean!r} then "
            f"{second.mean!r}. Results are not reproducible; stopping."
        )
    return f"ok: seed 99 reproducible at {first.mean:.6f} ({n_jobs:,} jobs)"


def simulate_cell(rho, r, n_jobs, seeds, cache, url):
    """All seeds for one cell, cached per run. Returns a row dict (or raises)."""
    lam, mu2 = rho * MU1, r * MU1
    means, runs = [], []
    fresh = 0
    for seed in seeds:
        key = run_key(rho, r, n_jobs, seed)
        record = cache.get(key)
        if record is None:
            fresh += 1
            print(f"    seed {seed} ...", end="", flush=True)
            started = time.time()
            run = qsim_fj.run_one(lam, MU1, mu2, seed=seed, n_jobs=n_jobs, url=url,
                                  name=f"fj-rho{rho:g}-r{r:g}-s{seed}")
            record = {
                "mean": run.mean, "intra_run_ci": run.ci,
                "samples_analyzed": run.samples_analyzed,
                "samples_discarded": run.samples_discarded,
                "success": run.success,
                "wall_clock_seconds": run.wall_clock_seconds,
                "round_trip_seconds": run.round_trip_seconds,
                "system_mean": run.system_mean,
                # `completed` is a diagnostic, not a gate (spec 8 as superseded by
                # 6.2a), but it is the field that distinguishes pre- and post-fix
                # service behaviour, so it belongs in the record rather than nowhere.
                "completed": run.completed,
            }
            cache.put(key, record)
            print(f" {run.mean:.4f} in {time.time() - started:.0f}s", flush=True)
        means.append(record["mean"])
        runs.append(record)

    t_qsim, hw_qsim = grand_mean_and_ci(means)
    engine_seconds = sum(x["wall_clock_seconds"] for x in runs)
    return {
        "rho": rho, "r": r, "lam": lam, "mu1": MU1, "mu2": mu2,
        "n_jobs": n_jobs, "seeds": list(seeds),
        "per_seed_means": means, "t_qsim": t_qsim, "hw_qsim": hw_qsim,
        "engine_seconds": engine_seconds,
        "round_trip_seconds": sum(x["round_trip_seconds"] for x in runs),
        "jobs_per_second": (n_jobs * len(seeds) / engine_seconds
                            if engine_seconds else None),
        "intra_run_half_widths": [
            None if x["intra_run_ci"] is None
            else (x["intra_run_ci"][1] - x["intra_run_ci"][0]) / 2 for x in runs],
        "samples_discarded": [x["samples_discarded"] for x in runs],
        "completed": [x.get("completed") for x in runs],
        # How many of this cell's runs executed now vs came from cache. The timing
        # comparison needs this: a cached qsim arm was measured in a different
        # process at a different time, which voids the interleaving guarantee.
        "fresh_runs": fresh,
    }


def gate_homogeneous(rows):
    """The r=1 cells must reproduce the exact closed form (spec 7.3).

    A gate that inspects zero cells is not a passing gate. If any r=1 cell is
    missing -- because its runs raised and were recorded as failures -- this
    reports a failure rather than vacuously approving the heterogeneous cells,
    which spec 7.3 says must not run when the oracle is unsatisfied.

    Three distinct outcomes, all of which leave the oracle UNSATISFIED and so all of
    which block the heterogeneous cells, but which must not be reported as the same
    thing: the cell is missing, the cell has no usable CI (undecidable), or the cell
    genuinely misses the closed form (failed). Only the last is evidence about qsim.
    """
    notes, failures = [], []
    checked = {row["rho"] for row in rows if row["r"] == R_CHECK}
    for rho in RHO_VALUES:
        if rho not in checked:
            failures.append(
                f"r=1 rho={rho}: cell missing, so the closed-form oracle could not "
                f"be evaluated (see the run failures above)"
            )
    for row in rows:
        if row["r"] != R_CHECK:
            continue
        exact = homogeneous_exact(row["rho"])
        # hw_a is the relative floor, not 0: see GATE_RELATIVE_FLOOR.
        gap, ratio, ok = agreement(exact, GATE_RELATIVE_FLOOR * exact,
                                   row["t_qsim"], row["hw_qsim"])
        # `ok` is TRISTATE. None means this cell has no usable CI -- a single seed
        # makes hw_qsim NaN -- which is undecidable, NOT a missed closed form. Bare
        # truthiness here treated the two identically, so `--seeds 0` reported an
        # EXACT match (gap +0.0000, ratio NaN) as "qsim does not reproduce the exact
        # homogeneous result ... nan tolerance units apart" and skipped all 12
        # heterogeneous cells. Both still block those cells; only one blames qsim.
        verdict = {True: "ok", False: "FAILED", None: "UNDECIDABLE"}[ok]
        tol = " n/a" if ok is None else f"{ratio:>4.2f}"
        notes.append(f"  rho={row['rho']:<5} exact={exact:9.4f} "
                     f"qsim={row['t_qsim']:9.4f} +/-{row['hw_qsim']:>7.4f}  "
                     f"gap={gap:+.4f} ({tol} tolerance units)  {verdict}")
        if ok is False:
            failures.append(f"r=1 rho={row['rho']}: qsim {row['t_qsim']:.4f} vs exact "
                            f"{exact:.4f}, {ratio:.2f} tolerance units apart")
        elif ok is None:
            failures.append(
                f"r=1 rho={row['rho']}: no usable CI, so the closed-form oracle is "
                f"UNDECIDABLE rather than failed -- this says nothing about qsim. "
                f"Observed gap {gap:+.4f} against exact {exact:.4f}. A single seed "
                f"yields no interval; re-run this cell with at least two seeds."
            )
    return notes, failures


def print_report(rows, reference, failures, gate_notes):
    print("\n" + "=" * 100)
    print("qsim-service vs. table1_repro simulator  (fixed-length runs, "
          "cross-seed 95% t-CI both arms)")
    print("=" * 100)
    for line in gate_notes:
        print(line)
    print("-" * 100)
    print(f"{'rho':>5} {'r':>2} | {'T_simple':>9} {'+-CI':>7} | "
          f"{'T_qsim':>9} {'+-CI':>7} | {'gap':>9} {'gap/hw':>7} {'agree':>6} | "
          f"{'jobs/s':>9}")
    print("-" * 100)
    disagreements, undecidable, uncompared = [], [], []
    for row in rows:
        ref = reference.get((row["rho"], row["r"]))
        if ref is None:
            print(f"{row['rho']:>5} {row['r']:>2} | {'-- no reference cell --':>60}")
            uncompared.append((row["rho"], row["r"]))
            continue
        t_ref, hw_ref = ref
        gap, ratio, agree = agreement(t_ref, hw_ref, row["t_qsim"], row["hw_qsim"])
        mark = {True: "yes", False: "NO", None: "?"}[agree]
        if agree is False:
            disagreements.append((row, t_ref, hw_ref, gap, ratio))
        elif agree is None:
            undecidable.append((row["rho"], row["r"]))
        print(f"{row['rho']:>5} {row['r']:>2} | {t_ref:>9.4f} {hw_ref:>7.4f} | "
              f"{row['t_qsim']:>9.4f} {row['hw_qsim']:>7.4f} | "
              f"{gap:>+9.4f} {ratio:>7.2f} {mark:>6} | "
              f"{row['jobs_per_second']:>9,.0f}")
    print("=" * 100)

    if disagreements:
        print(f"\nDISAGREEING CELLS ({len(disagreements)}):")
        for row, t_ref, hw_ref, gap, ratio in disagreements:
            print(f"  rho={row['rho']} r={row['r']}: simple {t_ref:.4f}+/-{hw_ref:.4f} "
                  f"vs qsim {row['t_qsim']:.4f}+/-{row['hw_qsim']:.4f}  "
                  f"gap {gap:+.4f} = {ratio:.2f} combined half-widths")
        print("  These are findings, not errors: the run itself is sound.")
        verdict_complete = True          # a disagreement IS a verdict
    elif undecidable or uncompared or not rows or failures:
        # An empty `disagreements` list is NOT evidence of agreement: a cell with no
        # reference entry never reaches the comparison, and a cell whose CI is
        # unavailable yields agree=None. Claiming agreement over those would assert
        # a result the data cannot support -- with an empty reference dict every
        # cell prints "no reference cell" and the old code still said ALL AGREE.
        decided = len(rows) - len(undecidable) - len(uncompared)
        print(f"\nNO VERDICT: {decided} of {len(rows)} cells were actually compared "
              f"and decided" + (f", and the run has {len(failures)} failure(s)"
                                if failures else "") + ".")
        if failures and len(rows) < len(RHO_VALUES) * (1 + len(R_VALUES)):
            print(f"  Only {len(rows)} of {RHO_VALUES.__len__() * (1 + len(R_VALUES))} "
                  f"grid cells ran at all, so this is a partial run.")
        if uncompared:
            print(f"  {len(uncompared)} cell(s) had no reference entry: {uncompared}")
        if undecidable:
            print(f"  {len(undecidable)} cell(s) had no usable CI, so agreement is "
                  f"undecidable: {undecidable}")
        if decided:
            print("  Of the cells that WERE decided, none disagreed -- but that is "
                  "not a whole-grid result.")
        # No whole-grid verdict was reached, so this run must not be mistaken for a
        # success by its exit status either (main turns this into a non-zero exit).
        verdict_complete = False
    else:
        print(f"\nALL CELLS AGREE: all {len(rows)} cells compared and decided; every "
              f"qsim mean is within the combined 95% CI of the simple simulator's.")
        verdict_complete = True

    if failures:
        print(f"\nRUN FAILURES ({len(failures)}) -- results above are not trustworthy:")
        for f in failures:
            print(f"  - {f}")

    print("\nThe two CIs compared above are built identically (cross-seed t over "
          "per-seed means).")
    print("qsim's own PER-RUN interval is a different object entirely:")
    width = max(len(a) for a, _, _ in CI_NATURE_ROWS)
    print(f"  {'':<{width}}   {'simple simulator':<28} qsim per-run CI")
    for label, simple, qsim in CI_NATURE_ROWS:
        print(f"  {label:<{width}}   {simple:<28} {qsim}")
    return disagreements, verdict_complete


_REPO_ROOT = os.path.join(SCRIPT_DIR, "..", "..")

BASELINE_WARMUP = 500_000
"""Warm-up jobs the reference simulator discards, matching reproduce_table1.py.

Note this is ADDITIONAL work: forkjoin.simulate runs `warmup + n_jobs` jobs and
keeps the last n_jobs. qsim has no warm-up knob and runs exactly n_jobs (spec 6.5),
so the two arms do not simulate identical totals -- 2.44% apart at 20 M jobs. The
timing table reports each arm's own total so the gap is visible rather than
implied.
"""


def time_baseline_cell(rho, r, n_jobs, seeds):
    """Re-time forkjoin.simulate fresh for one cell (spec 10).

    `forkjoin` is imported lazily -- inside the function, not at module scope -- so
    the qsim-only path works without numpy installed. The sys.path insertion is
    guarded by a membership test, so calling this once per cell cannot grow sys.path.
    """
    if _REPO_ROOT not in sys.path:
        sys.path.insert(0, _REPO_ROOT)
    from forkjoin import simulate

    lam, mu2 = rho * MU1, r * MU1
    started = time.time()
    means = []
    for seed in seeds:
        result = simulate(lam, MU1, mu2, n_jobs=n_jobs, warmup=BASELINE_WARMUP,
                          seed=seed)
        means.append(result.mean_response_time)
    elapsed = time.time() - started
    mean, hw = grand_mean_and_ci(means)
    # forkjoin.simulate computes total = warmup + n_jobs, so the warmup jobs are
    # SIMULATED IN ADDITION to n_jobs, not carved out of it. Throughput must
    # therefore be divided by the total simulated, or the simple simulator's rate
    # is understated -- by 2.5% at 20 M jobs and 3.5x at --quick's 200 K.
    total_per_seed = BASELINE_WARMUP + n_jobs
    return {"rho": rho, "r": r, "seconds": elapsed, "t_simple": mean, "hw_simple": hw,
            "warmup": BASELINE_WARMUP, "total_jobs_per_seed": total_per_seed,
            "jobs_per_second": total_per_seed * len(seeds) / elapsed}


def run_timing_comparison(rows, n_jobs, seeds, cache, url):
    """Interleaved, never concurrent: simple then qsim, cell by cell (spec 10).

    Interleaving matters because slow machine drift then hits both arms equally
    instead of penalising whichever ran second.
    """
    print("\nTIMING (never concurrent; qsim "
          f"{n_jobs:,} jobs/seed vs simple {BASELINE_WARMUP + n_jobs:,})")
    print(f"  {'rho':>5} {'r':>2} | {'simple s':>9} {'qsim s':>9} {'ratio':>7} | "
          f"{'simple jobs/s':>13} {'qsim jobs/s':>12} | "
          f"{'hw_simple':>9} {'hw_qsim':>9}")
    results = []
    for row in rows:
        rho, r = row["rho"], row["r"]
        base = time_baseline_cell(rho, r, n_jobs, seeds)          # simple first
        fresh = simulate_cell(rho, r, n_jobs, seeds, cache, url)  # then qsim
        ratio = fresh["engine_seconds"] / base["seconds"]
        entry = {**base, "qsim_seconds": fresh["engine_seconds"], "ratio": ratio,
                 "qsim_jobs_per_second": fresh["jobs_per_second"],
                 "hw_qsim": fresh["hw_qsim"],
                 "qsim_runs_fresh": fresh["fresh_runs"],
                 "qsim_runs_total": len(seeds),
                 "interleaved": fresh["fresh_runs"] == len(seeds)}
        results.append(entry)
        print(f"  {rho:>5} {r:>2} | {base['seconds']:>9.1f} "
              f"{fresh['engine_seconds']:>9.1f} {ratio:>7.1f}x | "
              f"{base['jobs_per_second']:>13,.0f} "
              f"{fresh['jobs_per_second']:>12,.0f} | "
              f"{base['hw_simple']:>9.4f} {fresh['hw_qsim']:>9.4f}")
    print(f"  NOT identical work: the simple simulator additionally runs "
          f"{BASELINE_WARMUP:,} warm-up jobs per seed "
          f"({BASELINE_WARMUP + n_jobs:,} total vs qsim's {n_jobs:,}), so read the "
          f"jobs/s columns -- which are work-normalised -- rather than the raw ratio.")
    stale = [(x["rho"], x["r"]) for x in results if not x["interleaved"]]
    if stale:
        print(f"  NOT INTERLEAVED for {len(stale)} of {len(results)} cells: their qsim "
              f"runs came from cache, i.e. were measured in an earlier process at an "
              f"earlier time, while the simple arm was re-timed just now. The "
              f"per-cell pairing that is supposed to cancel machine drift did NOT "
              f"happen for those cells. Delete qsim_results.json (or pass a fresh "
              f"--n-jobs) to force a genuinely interleaved measurement.")
    else:
        print(f"  Interleaved per cell: both arms of all {len(results)} cells were "
              f"measured fresh, alternating, in this process.")
    print("  Caveats: qsim's seconds are engine-side inside a podman VM (see the "
          "2-cell local-JVM probe in the README), and neither simulator removes the "
          "same transient (spec 6.5).")
    return results


TRANSIENT_SEED = 0
"""Single seed used for the N-vs-2N comparison."""


def run_transient_check(rows, n_jobs, cache, url):
    """Compare N against 2N at one seed: is residual transient bias visible?

    qsim exposes no warm-up PARAMETER. It does remove a transient -- JMT detects one
    per run, median 0.03% of the run against the reference's caller-fixed 2.44% (see
    CI_NATURE_ROWS) -- but the amount is neither settable nor equal to the
    reference's, so any residual bias has to be measured rather than assumed away.

    A shift smaller than the N-run's own intra-run CI means no bias is DETECTABLE at
    this sensitivity. That is a bound, not a demonstration of zero: see the
    SENSITIVITY note printed with the result (spec 6.5).
    """
    print(f"\nTRANSIENT-BIAS CHECK ({n_jobs:,} vs {2 * n_jobs:,} jobs, "
          f"seed {TRANSIENT_SEED})")
    print(f"  {'rho':>5} {'r':>2} | {'T(N)':>9} {'T(2N)':>9} {'diff':>9} "
          f"{'diff %':>8} {'within CI':>10}")
    results = []
    for row in rows:
        rho, r = row["rho"], row["r"]
        short = simulate_cell(rho, r, n_jobs, [TRANSIENT_SEED], cache, url)
        long = simulate_cell(rho, r, 2 * n_jobs, [TRANSIENT_SEED], cache, url)
        t_n, t_2n = short["t_qsim"], long["t_qsim"]
        diff = t_2n - t_n
        hw = short["intra_run_half_widths"][0]
        within = None if hw is None else abs(diff) <= hw
        mark = {True: "yes", False: "NO", None: "?"}[within]
        results.append({"rho": rho, "r": r, "t_n": t_n, "t_2n": t_2n, "diff": diff,
                        "intra_run_hw": hw, "within_ci": within})
        print(f"  {rho:>5} {r:>2} | {t_n:>9.4f} {t_2n:>9.4f} {diff:>+9.4f} "
              f"{100 * diff / t_n:>+7.2f}% {mark:>10}")
    outside = [x for x in results if x["within_ci"] is False]
    unknown = [(x["rho"], x["r"]) for x in results if x["within_ci"] is None]
    worst = max((abs(x["diff"]) / x["t_n"] * 100 for x in results), default=0.0)
    if outside:
        signs = {"low" if x["diff"] > 0 else "high" for x in outside}
        print(f"  {len(outside)} cell(s) shifted beyond the noise floor; the N-job "
              f"mean is biased {'/'.join(sorted(signs))}. Reportable bias, not an error.")
    elif unknown or len(results) < len(rows):
        # Same trap as the agreement verdict: absent and undecidable cells must not
        # be absorbed into a whole-grid negative result.
        print(f"  NO WHOLE-GRID VERDICT: {len(results) - len(unknown)} of {len(rows)} "
              f"cells were decided.")
        if unknown:
            print(f"    {len(unknown)} had no intra-run CI, so undecidable: {unknown}")
    else:
        print(f"  All {len(results)} cells: no shift detectable at this sensitivity.")
    # State the bound rather than implying "no bias". Two real limits:
    print(f"  SENSITIVITY: this compares against JMT's INTRA-RUN half-width, which "
          f"spec 6.1 treats as a diagnostic only and which runs up to 2.9x wider than "
          f"the cross-seed half-width the experiment actually reports. Largest observed "
          f"shift {worst:.2f}%. Supported claim: 'no shift larger than roughly that is "
          f"detectable here' -- NOT 'no transient bias'.")
    print(f"  Also: T(2N) at seed {TRANSIENT_SEED} extends the same sample path as T(N), so "
          f"if bias decays like 1/N this statistic sees only about half of it.")
    return results


def _tex_escape(text):
    """Escape the two LaTeX-special characters that appear in CI_NATURE_ROWS."""
    return text.replace("%", r"\%").replace("_", r"\_")


def _tex_ci_nature_table():
    """The spec 6.1 table as LaTeX, so it travels with the numbers it qualifies."""
    lines = [r"\begin{tabular}{@{}lll@{}}", r"\toprule",
             r"& Simple simulator & qsim per-run CI \\", r"\midrule"]
    for label, simple, qsim in CI_NATURE_ROWS:
        lines.append(f"{_tex_escape(label)} & {_tex_escape(simple)} & "
                     f"{_tex_escape(qsim)} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return lines


def write_outputs(rows, reference, make_figure=True, timing=None, transient=None,
                  url=None):
    """Write the JSON record, the LaTeX tables, and the agreement figure."""
    by_cell = {(row["rho"], row["r"]): row for row in rows}
    records = []
    for row in rows:
        ref = reference.get((row["rho"], row["r"]))
        record = dict(row)
        if ref is not None:
            t_ref, hw_ref = ref
            gap, ratio, agree = agreement(t_ref, hw_ref, row["t_qsim"], row["hw_qsim"])
            record.update({"t_simple": t_ref, "hw_simple": hw_ref, "gap": gap,
                           "gap_over_combined_hw": ratio, "agree": agree})
        records.append(record)
    out_json = os.path.join(SCRIPT_DIR, "table1_qsim_results.json")
    with open(out_json, "w") as f:
        json.dump({"rows": records,
                   "ci_nature": [list(r) for r in CI_NATURE_ROWS],
                   "timing": timing, "transient_check": transient,
                   "qsim_build": os.environ.get("QSIM_BUILD"),
                   # The URL actually used, not the env default: --url overrides it.
                   "qsim_url": url or qsim_fj.DEFAULT_URL}, f, indent=2)
    print(f"\nResults written to {out_json}")

    lines = [
        r"% Auto-generated by experiments/table1_qsim/validate_table1_qsim.py",
        r"% Both CIs are cross-seed 95% Student-t intervals over per-seed means.",
        r"\begin{table}[ht]", r"\centering",
        (r"\caption{Table~1's simulated sojourn times re-measured with an "
         r"independent discrete-event simulator (qsim-service / headless JMT). "
         r"Both columns are cross-seed 95\% Student-t intervals over five per-seed "
         r"means at a fixed run length. `Agree' is CI overlap: "
         r"$|T_{qsim}-T_{simple}| \le hw_{qsim}+hw_{simple}$.}"),
        r"\label{tab:qsim-validation}", r"\small{",
        r"\begin{tabular}{@{}ccccrrc@{}}", r"\toprule",
        (r"$\rho$ & $r$ & $T_\mathrm{simple}\pm$CI & $T_\mathrm{qsim}\pm$CI"
         r" & Gap & Gap/hw & Agree \\"),
        r"\midrule",
    ]
    for i, rho in enumerate(RHO_VALUES):
        if i > 0:
            lines.append(r"\midrule")
        for r in [R_CHECK] + list(R_VALUES):
            row = by_cell.get((rho, r))
            ref = reference.get((rho, r))
            if row is None or ref is None:
                continue
            t_ref, hw_ref = ref
            gap, ratio, agree = agreement(t_ref, hw_ref, row["t_qsim"], row["hw_qsim"])
            mark = {True: r"\checkmark", False: r"\textbf{no}", None: "?"}[agree]
            lines.append(
                f"${rho}$ & ${r}$ & ${t_ref:.4f} \\pm {hw_ref:.4f}$ & "
                f"${row['t_qsim']:.4f} \\pm {row['hw_qsim']:.4f}$ & "
                f"${gap:+.4f}$ & ${ratio:.2f}$ & {mark} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}%", r"}", r"\end{table}", "",
              r"% The two simulators' CIs differ in nature; see spec 6.1:", ""]
    lines += _tex_ci_nature_table()
    tex = os.path.join(SCRIPT_DIR, "table1_qsim_comparison.tex")
    with open(tex, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"LaTeX written to {tex}")

    if make_figure:
        _make_figure(by_cell, reference)


def _make_figure(by_cell, reference):
    """Bars of gap / combined half-width, with a shaded +/-1 agreement band."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping figure.")
        return

    all_r = [R_CHECK] + list(R_VALUES)
    fig, axes = plt.subplots(1, len(RHO_VALUES), figsize=(3.6 * len(RHO_VALUES), 3.4),
                             sharey=True)
    for ax, rho in zip(axes, RHO_VALUES):
        xs, heights, colors, labels = [], [], [], []
        for i, r in enumerate(all_r):
            row, ref = by_cell.get((rho, r)), reference.get((rho, r))
            if row is None or ref is None:
                continue
            gap, ratio, agree = agreement(ref[0], ref[1], row["t_qsim"], row["hw_qsim"])
            signed = ratio if gap >= 0 else -ratio
            xs.append(i)
            # An undecidable cell (agree None, ratio NaN) is drawn as a zero-height
            # GRAY bar labelled "n/a", never firebrick: colouring it like a genuine
            # disagreement would invert the reader's conclusion, the same trap the
            # bar labels below already guard against for near-zero bars.
            heights.append(0.0 if agree is None else signed)
            colors.append({True: "steelblue", False: "firebrick",
                           None: "gray"}[agree])
            labels.append("n/a" if agree is None else f"{abs(signed):.2f}")
        ax.axhspan(-1, 1, color="lightgray", alpha=0.6, zorder=0,
                   label="agreement band")
        ax.bar(xs, heights, 0.6, color=colors, zorder=2)
        # Label every bar. Without this a near-zero bar (excellent agreement, e.g.
        # 0.01) is visually indistinguishable from an absent cell, which inverts the
        # reader's conclusion.
        for x, h, label in zip(xs, heights, labels):
            ax.annotate(label, (x, h),
                        textcoords="offset points", xytext=(0, 3 if h >= 0 else -11),
                        ha="center", fontsize=7, zorder=4)
        ax.axhline(0, color="k", linewidth=0.7, zorder=3)
        ax.set_title(rf"$\rho = {rho}$")
        ax.set_xticks(range(len(all_r)))
        ax.set_xticklabels([f"$r={r}$" for r in all_r])
        ax.grid(True, axis="y", linestyle="--", linewidth=0.5, color="lightgray")
    axes[0].set_ylabel("gap / combined 95% CI half-width")
    axes[0].legend(fontsize=8, loc="best")
    fig.suptitle("qsim-service vs. table1_repro simulator: inside the band = agree")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = os.path.join(SCRIPT_DIR, f"table1_qsim_agreement.{ext}")
        fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Figure written to {SCRIPT_DIR}/table1_qsim_agreement.{{png,pdf}}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--n-jobs", type=int, default=None)
    parser.add_argument("--quick", action="store_true",
                        help="200K jobs per seed (smoke test, still fixed-length)")
    parser.add_argument("--seeds", type=int, nargs="+", default=SEEDS)
    parser.add_argument("--url", default=qsim_fj.DEFAULT_URL)
    parser.add_argument("--no-fig", action="store_true")
    parser.add_argument("--skip-gates", action="store_true",
                        help="skip health/determinism/r=1 gates (debugging only)")
    parser.add_argument("--time-baseline", action="store_true",
                        help="re-time forkjoin.simulate and compare cost (spec 10)")
    parser.add_argument("--transient-check", action="store_true",
                        help="compare N vs 2N to quantify un-removed warm-up (spec 6.5)")
    args = parser.parse_args(argv)

    if args.n_jobs is not None and args.n_jobs <= 0:
        parser.error(f"--n-jobs must be > 0, got {args.n_jobs}")
    if len(set(args.seeds)) != len(args.seeds):
        # Checked BEFORE the count guard: with `--seeds 0 0 ... 0` the count message
        # would blame the arity when the real defect is the repetition.
        dupes = sorted({s for s in args.seeds if args.seeds.count(s) > 1})
        parser.error(
            f"--seeds: repeated seed(s) {dupes}. Each seed is one INDEPENDENT "
            f"replication -- the cross-seed t-interval is valid only because the "
            f"per-seed means are i.i.d. A repeat re-reads the same run from cache, so "
            f"it adds no information while shrinking the sample variance: "
            f"`--seeds 0 0 0` yields s=0 and reports a ZERO-WIDTH CI as certainty."
        )
    if len(args.seeds) > max(T_975) + 1:
        # Fail here rather than after hours of simulation: grand_mean_and_ci has no
        # t quantile past this many replications, and it is called only at the end.
        parser.error(f"--seeds: {len(args.seeds)} seeds given, but the t table covers "
                     f"at most {max(T_975) + 1}; add quantiles to T_975 to go higher")
    n_jobs = args.n_jobs or (200_000 if args.quick else DEFAULT_N_JOBS)
    mismatch = protocol_mismatch_warning(n_jobs)
    if mismatch:
        print(mismatch)
    build = os.environ.get("QSIM_BUILD")
    # Per-run CIs differ 1.196x between service builds 0.1.0 and 0.2.0, and
    # --transient-check gates on exactly that CI, so a cache mixing builds would be
    # undetectable. Recording the build makes it checkable rather than inferable.
    print(f"qsim build: {build or 'UNRECORDED -- set QSIM_BUILD to stamp results'}")
    cache = RunCache(os.path.join(SCRIPT_DIR, "qsim_results.json"))
    reference = load_reference()
    print(f"Protocol: {len(args.seeds)} seeds {args.seeds} x {n_jobs:,} jobs "
          f"(fixed length), url={args.url}")

    gate_notes, failures = [], []
    if not args.skip_gates:
        qsim_fj.health(args.url)          # raises with start instructions
        gate_notes.append("health: ok")
        gate_notes.append(f"determinism: {check_determinism(args.url)}")

    rows = []
    print("\nHomogeneous r=1 cells (exact-closed-form gate):")
    for rho in RHO_VALUES:
        try:
            rows.append(simulate_cell(rho, R_CHECK, n_jobs, args.seeds, cache, args.url))
        except qsim_fj.QsimError as exc:
            failures.append(f"rho={rho} r={R_CHECK}: {exc}")
    hom_notes, hom_failures = gate_homogeneous(rows)
    gate_notes.append("r=1 closed-form gate (Nelson & Tantawi):")
    gate_notes.extend(hom_notes)
    if hom_failures and not args.skip_gates:
        print_report(rows, reference, failures + hom_failures, gate_notes)
        # Do not assert a cause here: the gate is equally unsatisfied by a missing
        # cell and by an undecidable one, neither of which is evidence that qsim is
        # wrong. The per-cell verdicts above say which of the three it was.
        print("\nGATE NOT SATISFIED: the r=1 closed-form oracle did not pass (see the "
              "per-cell verdicts above), so the heterogeneous cells were not run.")
        return 1
    failures.extend(hom_failures)

    print("\nHeterogeneous cells:")
    for rho in RHO_VALUES:
        for r in R_VALUES:
            print(f"  rho={rho} r={r}")
            try:
                rows.append(simulate_cell(rho, r, n_jobs, args.seeds, cache, args.url))
            except qsim_fj.QsimError as exc:
                failures.append(f"rho={rho} r={r}: {exc}")

    _, verdict_complete = print_report(rows, reference, failures, gate_notes)
    timing = (run_timing_comparison(rows, n_jobs, args.seeds, cache, args.url)
              if args.time_baseline else None)
    transient = (run_transient_check(rows, n_jobs, cache, args.url)
                 if args.transient_check else None)
    write_outputs(rows, reference, make_figure=not args.no_fig,
                  timing=timing, transient=transient, url=args.url)
    # Disagreement is a finding (exit 0). A run that reached no whole-grid verdict is
    # NOT a success, though: an empty or mis-keyed reference would otherwise be
    # indistinguishable from agreement by exit status alone.
    if not verdict_complete:
        print("\nExiting non-zero: no whole-grid verdict was reached (see NO VERDICT "
              "above). This is not a disagreement -- it is an unusable run.")
    return 0 if (verdict_complete and not failures) else 1


if __name__ == "__main__":
    raise SystemExit(main())
