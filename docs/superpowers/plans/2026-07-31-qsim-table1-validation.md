# qsim-service Table 1 Validation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Re-measure Table 1's 16 fork-join cells with qsim-service (headless JMT over
HTTP) and report, cell by cell, whether the two independent simulators agree within their
confidence intervals.

**Architecture:** Two stdlib-only modules under `experiments/table1_qsim/`. `qsim_fj.py`
turns `(lambda, mu1, mu2, seed, n_jobs)` into one validated `FJRun` record and knows
nothing about Table 1. `validate_table1_qsim.py` walks the grid, replicates over 5 seeds,
aggregates a cross-seed Student-t CI, compares against `table1_repro`'s cached per-seed
means, and emits console / LaTeX / figure output. No HTTP knowledge in the driver, no
Table 1 knowledge in the client.

**Tech Stack:** Python 3.10+ stdlib (`urllib`, `json`, `dataclasses`, `argparse`),
matplotlib (already a project dependency), pytest 9.0.3 for tests, qsim-service 0.1.0 in
podman.

> **Note (2026-08-01): this plan was executed; two of its claims were later falsified.**
> Both are corrected in the spec and README rather than here, since the plan is kept as
> the historical record of what was intended:
> (a) *"qsim performs no transient removal / `samplesDiscarded` is always 0"* — false at
> the 20 M protocol length, where 75 of 80 runs discarded one (median 6,490 samples);
> (b) *"this **is** equal work"* for the timing comparison — it is not, because
> `forkjoin.simulate` runs `warmup + n_jobs`, i.e. 20.5 M against qsim's 20 M.
> A third item was amended during execution: the `r=1` gate gained a 1% relative
> tolerance floor.

## Global Constraints

- **Spec is authoritative:** `docs/superpowers/specs/2026-07-31-qsim-table1-validation-design.md`. Section references below (§N) point at it.
- **No new dependencies.** `pyproject.toml` is not modified. `qopt` is never imported at runtime; only inside `--cross-check-qopt`, guarded by `try/except ImportError`.
- **Fixed-length runs only.** Every request sends `{"disableStatisticStop": true, "maxSamples": N}`. Never send `minSamples` (dropped by the service — qsim-service#10) and never rely on convergence stopping (§6.2 finding 2).
- **`measures` is always exactly `["response-time", "system-response-time"]`.** Never omit it; never add `utilization` or `queue-length` (join-station numbers at a fork-join node, qsim-service#8).
- **Run length is verified by `samplesAnalyzed == N`, never by `success`.** `success` and `precision` are recorded as diagnostics and gate nothing (§6.2 finding 2, §8).
- **Parameterization:** `MU1 = 1.0`, `mu2 = r * MU1`, `lam = rho * MU1`.
- **Grid:** `RHO_VALUES = [0.4, 0.8, 0.9, 0.95]`, `R_VALUES = [2, 4, 8]`, plus `R_CHECK = 1`.
- **Protocol:** `SEEDS = [0, 1, 2, 3, 4]`, `DEFAULT_N_JOBS = 20_000_000` (§6.3).
- **Service URL** comes from `QSIM_URL`, default `http://localhost:8080`. The code never starts or stops the service (§8).
- **Exit codes:** disagreeing cells exit **0** (a finding). Non-zero only for a failed gate, an HTTP 500, or an unreachable service (§8).
- **Never run the two simulators concurrently**, and interleave them per cell when timing (§10).

## File Structure

| File | Responsibility |
|---|---|
| `experiments/table1_qsim/qsim_fj.py` | One fork-join run: model JSON, POST, error mapping, `FJRun`. Plus a `__main__` CLI for single-run probes. |
| `experiments/table1_qsim/validate_table1_qsim.py` | Grid walk, replication, aggregation, comparison, reporting, timing. |
| `experiments/table1_qsim/expected_model.json` | Fixture pinning the emitted model JSON. |
| `experiments/table1_qsim/test_qsim_fj.py` | Unit tests for the client (fake transport, no network). |
| `experiments/table1_qsim/test_validate.py` | Unit tests for aggregation, reference loading, closed form, cache. |
| `experiments/table1_qsim/README.md` | Protocol differences, CI-nature table, how to run. |
| `experiments/table1_qsim/qsim_results.json` | Run cache (generated). |
| `experiments/table1_qsim/table1_qsim_comparison.tex` | Agreement + CI-nature tables (generated). |
| `experiments/table1_qsim/table1_qsim_agreement.{png,pdf}` | Agreement figure (generated). |

Tests live beside the modules because the repo has no `tests/` directory and is
script-oriented (`verify_formulas.py` at root). With pytest's default `prepend` import
mode, `import qsim_fj` resolves from a test file in the same directory.

---

### Task 1: Model builder and fixture

**Files:**
- Create: `experiments/table1_qsim/qsim_fj.py`
- Create: `experiments/table1_qsim/expected_model.json`
- Test: `experiments/table1_qsim/test_qsim_fj.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `build_model(lam, mu1, mu2, *, name=None, job_class="jobs") -> dict`,
  `build_request(lam, mu1, mu2, *, seed, n_jobs, alpha=0.05, max_wall_clock=1800, name=None) -> dict`,
  `canonical(obj) -> str`, `MEASURES: tuple[str, str]`.

- [ ] **Step 1: Write the failing test**

```python
# experiments/table1_qsim/test_qsim_fj.py
import json
import os

import pytest

import qsim_fj

HERE = os.path.dirname(os.path.abspath(__file__))


def test_build_model_matches_fixture():
    model = qsim_fj.build_model(0.8, 1.0, 4.0, name="fj-rho0.80-r4")
    with open(os.path.join(HERE, "expected_model.json")) as f:
        expected = json.load(f)
    assert qsim_fj.canonical(model) == qsim_fj.canonical(expected)


def test_build_model_topology():
    model = qsim_fj.build_model(0.8, 1.0, 4.0)
    names = [n["name"] for n in model["nodes"]]
    assert names == ["src", "fj", "snk"]
    fj = model["nodes"][1]
    assert fj["type"] == "fork-join"
    assert fj["join"] == "all"
    rates = [b["service"]["jobs"]["distribution"]["rate"] for b in fj["branches"]]
    assert rates == [1.0, 4.0]
    assert model["nodes"][0]["arrivals"]["jobs"]["distribution"]["rate"] == 0.8


def test_build_request_is_fixed_length_only():
    req = qsim_fj.build_request(0.8, 1.0, 4.0, seed=3, n_jobs=20_000_000)
    stopping = req["stopping"]
    assert stopping["disableStatisticStop"] is True
    assert stopping["maxSamples"] == 20_000_000
    # qsim-service#10: minSamples is silently dropped by the service, so never send it.
    assert "minSamples" not in stopping
    # Convergence stopping is unusable (spec 6.2 finding 2), so no precision target.
    assert "precision" not in stopping
    assert req["seed"] == 3
    assert req["measures"] == ["response-time", "system-response-time"]


def test_build_request_rejects_bad_inputs():
    with pytest.raises(ValueError):
        qsim_fj.build_request(1.5, 1.0, 4.0, seed=0, n_jobs=1000)   # lam >= mu1: unstable
    with pytest.raises(ValueError):
        qsim_fj.build_request(0.8, 1.0, 4.0, seed=0, n_jobs=0)      # n_jobs must be > 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd experiments/table1_qsim && python -m pytest test_qsim_fj.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'qsim_fj'`

- [ ] **Step 3: Write minimal implementation**

```python
# experiments/table1_qsim/qsim_fj.py
#!/usr/bin/env python3
"""One heterogeneous 2-queue fork-join run against qsim-service.

Stdlib only, by design: this experiment must not add a dependency on the
unreleased `qopt` sibling project (spec 4.2). The model JSON below is transcribed
from qopt's own emitters --

    qopt/station.py :: ForkJoinStation.sim_node
    qopt/network.py :: Network.to_model_dict

-- and `python qsim_fj.py --cross-check-qopt` asserts the two still agree, so the
transcription cannot drift silently.

Protocol note (spec 6.2, qsim-service#10): only FIXED-LENGTH runs are trustworthy.
`minSamples` is accepted by the API and never reaches JMT, and a run that stops
early still reports `success: true`, so run length is verified by comparing
`samplesAnalyzed` against the requested N and by nothing else.
"""

import json

MEASURES = ("response-time", "system-response-time")
"""Closed, always-explicit measure list.

Omitting `measures` makes the service substitute defaults that include
`utilization` and `queue-length`, which at a fork-join node are *join-station*
numbers returned with success:true and no warning (qsim-service#8).
"""

SOURCE, FJ, SINK = "src", "fj", "snk"
JOB_CLASS = "jobs"


def canonical(obj):
    """Key-order-independent JSON rendering, for comparing two model dicts."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"))


def build_model(lam, mu1, mu2, *, name=None, job_class=JOB_CLASS):
    """The qsim `model` block for src -> fork-join(mu1, mu2) -> snk."""
    def exponential(rate):
        return {"distribution": {"type": "exponential", "rate": rate}}

    return {
        "name": name or f"fj-lam{lam:g}-mu{mu1:g}-{mu2:g}",
        "classes": [{"name": job_class, "type": "open"}],
        "nodes": [
            {"name": SOURCE, "type": "source",
             "arrivals": {job_class: exponential(lam)}},
            {"name": FJ, "type": "fork-join", "join": "all",
             "branches": [{"service": {job_class: exponential(mu1)}},
                          {"service": {job_class: exponential(mu2)}}]},
            {"name": SINK, "type": "sink"},
        ],
        "routing": {job_class: [
            {"from": SOURCE, "to": FJ, "probability": 1.0},
            {"from": FJ, "to": SINK, "probability": 1.0},
        ]},
    }


def build_request(lam, mu1, mu2, *, seed, n_jobs, alpha=0.05,
                  max_wall_clock=1800, name=None):
    """Full POST /simulate body for one fixed-length run of exactly `n_jobs` jobs."""
    if not (0 < lam < min(mu1, mu2)):
        raise ValueError(
            f"unstable or invalid: need 0 < lam < min(mu1, mu2), got lam={lam}, "
            f"mu1={mu1}, mu2={mu2}"
        )
    if n_jobs <= 0:
        raise ValueError(f"n_jobs must be > 0, got {n_jobs}")
    return {
        "model": build_model(lam, mu1, mu2, name=name),
        "seed": seed,
        # No minSamples (dropped by the service) and no precision target: fixed
        # length is the only run-length control that works. See module docstring.
        "stopping": {"alpha": alpha,
                     "maxSamples": n_jobs,
                     "disableStatisticStop": True,
                     "maxWallClockSeconds": max_wall_clock},
        "measures": list(MEASURES),
    }
```

- [ ] **Step 4: Generate the fixture, then run tests to verify they pass**

Generate `expected_model.json` from the implementation, then eyeball it against spec §5
before trusting it:

```bash
cd experiments/table1_qsim
python -c "
import json, qsim_fj
m = qsim_fj.build_model(0.8, 1.0, 4.0, name='fj-rho0.80-r4')
json.dump(m, open('expected_model.json','w'), indent=2, sort_keys=True)
"
python -c "
import json
m = json.load(open('expected_model.json'))
fj = [n for n in m['nodes'] if n['name']=='fj'][0]
assert fj['type']=='fork-join' and fj['join']=='all'
assert [b['service']['jobs']['distribution']['rate'] for b in fj['branches']]==[1.0,4.0]
print('fixture sane')
"
python -m pytest test_qsim_fj.py -v
```

Expected: `fixture sane`, then 4 passed.

- [ ] **Step 5: Commit**

```bash
git add experiments/table1_qsim/qsim_fj.py \
        experiments/table1_qsim/expected_model.json \
        experiments/table1_qsim/test_qsim_fj.py
git commit -m "Add fork-join model builder for qsim validation

Fixed-length requests only: no minSamples (qsim-service#10) and no
precision target, since convergence stopping is unusable. Measure list
is closed to avoid join-station numbers (qsim-service#8)."
```

---

### Task 2: HTTP client, error mapping, and per-run gates

**Files:**
- Modify: `experiments/table1_qsim/qsim_fj.py`
- Test: `experiments/table1_qsim/test_qsim_fj.py`

**Interfaces:**
- Consumes: `build_request`, `MEASURES` from Task 1.
- Produces: `FJRun` dataclass with fields `mean, ci, samples_analyzed, samples_discarded, success, completed, wall_clock_seconds, round_trip_seconds, system_mean`; exceptions `QsimError`, `QsimTransportError`, `QsimRequestError`, `QsimEngineError`, `QsimRunLengthError`, `QsimSemanticsError`; functions `health(url, timeout=10, transport=None) -> dict` and `run_one(lam, mu1, mu2, *, seed, n_jobs, url=DEFAULT_URL, transport=None, ...) -> FJRun`; constant `DEFAULT_URL`; `ORACLE_TOLERANCE = 0.02`.

- [ ] **Step 1: Write the failing test**

Append to `test_qsim_fj.py`:

```python
def _response(mean=5.0585, n=20_000_000, system_mean=None, completed=True,
              success=True, discarded=0):
    """A well-formed /simulate response body, shaped like the real service's."""
    if system_mean is None:
        system_mean = mean

    def measure(station, mtype, value):
        return {"station": station, "type": mtype, "class": "jobs",
                "mean": value, "lower": value - 0.07, "upper": value + 0.07,
                "alpha": 0.05, "precision": 0.005, "success": success,
                "samplesAnalyzed": n, "samplesDiscarded": discarded,
                "variance": None, "stdDev": None}

    return {"modelName": "t", "wallClockSeconds": 108.5, "completed": completed,
            "measures": [measure("fj", "response-time", mean),
                         measure("", "system-response-time", system_mean)]}


def _transport(status=200, body=None, capture=None):
    """Fake transport: records the decoded request, returns a canned response."""
    def fake(url, data, timeout):
        if capture is not None:
            capture.append((url, json.loads(data) if data else None))
        payload = _response() if body is None else body
        return status, json.dumps(payload).encode()
    return fake


def test_run_one_extracts_fj_measure():
    run = qsim_fj.run_one(0.8, 1.0, 4.0, seed=0, n_jobs=20_000_000,
                          transport=_transport())
    assert run.mean == pytest.approx(5.0585)
    assert run.samples_analyzed == 20_000_000
    assert run.samples_discarded == 0
    assert run.completed is True
    assert run.wall_clock_seconds == pytest.approx(108.5)
    assert run.round_trip_seconds >= 0.0
    assert run.ci == pytest.approx((5.0585 - 0.07, 5.0585 + 0.07))


def test_run_one_posts_to_simulate_endpoint():
    capture = []
    qsim_fj.run_one(0.8, 1.0, 4.0, seed=7, n_jobs=1000,
                    url="http://h:9", transport=_transport(capture=capture))
    url, body = capture[0]
    assert url == "http://h:9/simulate"
    assert body["seed"] == 7
    assert body["stopping"]["maxSamples"] == 1000


def test_run_one_rejects_short_run():
    """samplesAnalyzed != N is the ONLY trustworthy run-length check (spec 6.2)."""
    short = _response(n=10_880, success=True)      # the qsim-service#10 symptom
    with pytest.raises(qsim_fj.QsimRunLengthError) as exc:
        qsim_fj.run_one(0.8, 1.0, 4.0, seed=0, n_jobs=1_000_000,
                        transport=_transport(body=short))
    assert "10880" in str(exc.value).replace(",", "")


def test_run_one_rejects_incomplete_run():
    with pytest.raises(qsim_fj.QsimEngineError):
        qsim_fj.run_one(0.8, 1.0, 4.0, seed=0, n_jobs=20_000_000,
                        transport=_transport(body=_response(completed=False)))


def test_run_one_enforces_system_response_time_oracle():
    """One station, so system response time must equal the fork-join's (spec 7.4)."""
    bad = _response(mean=5.0, system_mean=6.0)     # 20% apart: wrong semantics
    with pytest.raises(qsim_fj.QsimSemanticsError):
        qsim_fj.run_one(0.8, 1.0, 4.0, seed=0, n_jobs=20_000_000,
                        transport=_transport(body=bad))


def test_run_one_tolerates_success_false():
    """disableStatisticStop makes success:false expected; it must not fail the run."""
    run = qsim_fj.run_one(0.8, 1.0, 4.0, seed=0, n_jobs=20_000_000,
                          transport=_transport(body=_response(success=False)))
    assert run.success is False
    assert run.mean == pytest.approx(5.0585)


@pytest.mark.parametrize("status,expected", [
    (400, "QsimRequestError"), (422, "QsimRequestError"),
    (500, "QsimEngineError"), (503, "QsimEngineError"),
])
def test_run_one_maps_http_status(status, expected):
    body = {"error": "boom", "details": ["d1"]}
    with pytest.raises(getattr(qsim_fj, expected)):
        qsim_fj.run_one(0.8, 1.0, 4.0, seed=0, n_jobs=1000,
                        transport=_transport(status=status, body=body))


def test_health_rejects_non_200():
    with pytest.raises(qsim_fj.QsimTransportError):
        qsim_fj.health("http://h:9", transport=_transport(status=500, body={}))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd experiments/table1_qsim && python -m pytest test_qsim_fj.py -v`
Expected: FAIL — `AttributeError: module 'qsim_fj' has no attribute 'run_one'`

- [ ] **Step 3: Write minimal implementation**

Add to `qsim_fj.py` (imports go at the top with the existing `import json`):

```python
import os
import time
import urllib.error
import urllib.request
from dataclasses import dataclass

DEFAULT_URL = os.environ.get("QSIM_URL", "http://localhost:8080")
ORACLE_TOLERANCE = 0.02
"""Max relative gap between the fj and system response times (spec 7.4).

Loose on purpose: each measure has its own analyzer, so this guards against wrong
*semantics* (a join-residence number instead of a fork-join-region number), not
against imprecision. In practice the probe returned the two identical.
"""

_REQUEST_STATUSES = (400, 405, 413, 422)


class QsimError(RuntimeError):
    """Base for every qsim failure."""


class QsimTransportError(QsimError):
    """Could not reach the service at all."""


class QsimRequestError(QsimError):
    """4xx: our JSON is wrong, so every cell will fail identically."""


class QsimEngineError(QsimError):
    """5xx, an unreadable body, or completed:false."""


class QsimRunLengthError(QsimError):
    """samplesAnalyzed != requested N, so the run is not comparable to its siblings."""


class QsimSemanticsError(QsimError):
    """system-response-time disagreed with the fork-join response-time."""


@dataclass(frozen=True)
class FJRun:
    """One fixed-length fork-join run.

    `ci` is JMT's *intra-run* interval (batch-means/spectral). It is a diagnostic
    only -- the reported CI for a cell is the cross-seed Student-t interval built
    by the driver, which is the same construction the simple simulator uses
    (spec 6.1). Likewise `success` gates nothing (spec 6.2 finding 2).
    """
    mean: float
    ci: tuple
    samples_analyzed: int
    samples_discarded: int
    success: bool
    completed: bool
    wall_clock_seconds: float
    round_trip_seconds: float
    system_mean: float


def _urllib_transport(url, data, timeout):
    """POST when `data` is bytes, GET when None. Returns (status, body_bytes)."""
    request = urllib.request.Request(
        url, data=data, method="GET" if data is None else "POST",
        headers={} if data is None else {"Content-Type": "application/json"},
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return response.status, response.read()
    except urllib.error.HTTPError as exc:
        return exc.code, exc.read()
    except (urllib.error.URLError, OSError) as exc:
        raise QsimTransportError(f"{url}: {exc}") from exc


def _decode(raw):
    try:
        return json.loads(raw)
    except (ValueError, TypeError) as exc:
        raise QsimEngineError(f"unreadable response body: {raw[:200]!r}") from exc


def _detail(raw):
    payload = None
    try:
        payload = json.loads(raw)
    except (ValueError, TypeError):
        return repr(raw[:200])
    if isinstance(payload, dict):
        details = payload.get("details") or []
        joined = "; ".join(details) if isinstance(details, list) else str(details)
        return f"{payload.get('error')}: {joined}" if joined else str(payload.get("error"))
    return repr(raw[:200])


def health(url=DEFAULT_URL, timeout=10, transport=None):
    """One GET /health, so a bad URL fails before the first cell (spec 7.1)."""
    transport = transport or _urllib_transport
    status, raw = transport(f"{url.rstrip('/')}/health", None, timeout)
    if status != 200:
        raise QsimTransportError(
            f"{url}/health returned HTTP {status}: {raw[:200]!r}\n"
            f"Start the service with one of:\n"
            f"  podman run -d --rm --name qsim -p 8080:8080 qsim-service:0.1.0\n"
            f"  java -cp 'target/qsim-service.jar:target/dependency/*:"
            f"lib/JMT-singlejar-1.4.0.jar' qsim.http.App"
        )
    return _decode(raw)


def _measure(response, station, mtype):
    for m in response.get("measures", []):
        if m.get("station") == station and m.get("type") == mtype:
            return m
    return None


def run_one(lam, mu1, mu2, *, seed, n_jobs, url=DEFAULT_URL, transport=None,
            alpha=0.05, max_wall_clock=1800, name=None, timeout=None):
    """Run exactly `n_jobs` jobs and return a validated FJRun."""
    transport = transport or _urllib_transport
    request = build_request(lam, mu1, mu2, seed=seed, n_jobs=n_jobs, alpha=alpha,
                            max_wall_clock=max_wall_clock, name=name)
    body = json.dumps(request).encode("utf-8")
    started = time.monotonic()
    status, raw = transport(f"{url.rstrip('/')}/simulate", body,
                            max_wall_clock + 60 if timeout is None else timeout)
    round_trip = time.monotonic() - started

    if status != 200:
        detail = _detail(raw)
        if status in _REQUEST_STATUSES:
            raise QsimRequestError(
                f"HTTP {status} from /simulate: {detail}\n"
                f"request was: {json.dumps(request)[:1000]}"
            )
        if 500 <= status < 600:
            raise QsimEngineError(f"HTTP {status} from /simulate: {detail}")
        raise QsimTransportError(f"unexpected HTTP {status} from /simulate: {detail}")

    response = _decode(raw)
    if not response.get("completed", True):
        raise QsimEngineError(
            f"run reported completed:false (watchdog fired); a fixed-length run has "
            f"no usable partial result"
        )

    fj = _measure(response, FJ, "response-time")
    if fj is None or fj.get("mean") is None:
        raise QsimEngineError(f"no 'response-time' for station {FJ!r} in response")

    analyzed = fj.get("samplesAnalyzed")
    if analyzed != n_jobs:
        raise QsimRunLengthError(
            f"requested {n_jobs} jobs but samplesAnalyzed={analyzed}; the run is not "
            f"comparable to its siblings. Note `success` is not a run-length signal "
            f"(qsim-service#10)."
        )

    system = _measure(response, "", "system-response-time")
    system_mean = None if system is None else system.get("mean")
    if system_mean is None:
        raise QsimEngineError(
            "no 'system-response-time' in response; the single-station oracle "
            "(spec 7.4) cannot run"
        )
    gap = abs(system_mean - fj["mean"]) / fj["mean"]
    if gap > ORACLE_TOLERANCE:
        raise QsimSemanticsError(
            f"system-response-time {system_mean:.6f} differs from fork-join "
            f"response-time {fj['mean']:.6f} by {100 * gap:.2f}% (> "
            f"{100 * ORACLE_TOLERANCE:g}%). This network has one station, so they "
            f"must agree: the fork-join measure may not carry region semantics."
        )

    lower, upper = fj.get("lower"), fj.get("upper")
    return FJRun(
        mean=fj["mean"],
        ci=None if lower is None or upper is None else (lower, upper),
        samples_analyzed=analyzed,
        samples_discarded=fj.get("samplesDiscarded"),
        success=bool(fj.get("success", True)),
        completed=True,
        wall_clock_seconds=response.get("wallClockSeconds"),
        round_trip_seconds=round_trip,
        system_mean=system_mean,
    )
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd experiments/table1_qsim && python -m pytest test_qsim_fj.py -v`
Expected: 14 passed (4 from Task 1 + 10 here).

- [ ] **Step 5: Commit**

```bash
git add experiments/table1_qsim/qsim_fj.py experiments/table1_qsim/test_qsim_fj.py
git commit -m "Add qsim HTTP client with per-run validation gates

Run length is checked against samplesAnalyzed, never against success
(qsim-service#10 showed a 10,880-sample run reporting success:true).
Enforces the single-station system-response-time oracle, and treats
success:false as expected under disableStatisticStop."
```

---

### Task 3: Single-run CLI, qopt cross-check, and live smoke test

**Files:**
- Modify: `experiments/table1_qsim/qsim_fj.py`

**Interfaces:**
- Consumes: `run_one`, `build_model`, `canonical`, `health` from Tasks 1-2.
- Produces: `python qsim_fj.py --rho R --r R [--n-jobs N] [--seed S] [--url U] [--cross-check-qopt]`.

- [ ] **Step 1: Write the CLI**

Append to `qsim_fj.py`:

```python
def cross_check_qopt(lam, mu1, mu2, name):
    """Assert our transcribed model still matches qopt's emitter (spec 4.2).

    Compared as canonical JSON rather than raw bytes: key order is not part of the
    contract, values are. Returns a status string; raises AssertionError on drift.
    """
    try:
        from qopt import ForkJoinStation, Network, Route
    except ImportError as exc:
        return f"skipped: qopt not importable ({exc})"

    network = Network(
        [ForkJoinStation(mu=mu1, r=mu2 / mu1, c1=1.0, c2=1.0, name=FJ)],
        [Route(Network.SOURCE, FJ, 1.0), Route(FJ, Network.SINK, 1.0)],
        arrival_rate=lam, name=name,
    )
    theirs = network.to_model_dict([1.0])
    ours = build_model(lam, mu1, mu2, name=name)
    if canonical(ours) != canonical(theirs):
        raise AssertionError(
            "model JSON has drifted from qopt's emitter.\n"
            f"ours:   {canonical(ours)}\n"
            f"theirs: {canonical(theirs)}"
        )
    return "ok: identical to qopt's emitter"


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--rho", type=float, required=True, help="lam / mu1")
    parser.add_argument("--r", type=float, required=True, help="mu2 / mu1")
    parser.add_argument("--n-jobs", type=int, default=2_000_000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--url", default=DEFAULT_URL)
    parser.add_argument("--cross-check-qopt", action="store_true")
    args = parser.parse_args(argv)

    mu1, mu2, lam = 1.0, args.r, args.rho
    name = f"fj-rho{args.rho:.2f}-r{args.r:g}"

    if args.cross_check_qopt:
        print(f"qopt cross-check: {cross_check_qopt(lam, mu1, mu2, name)}")

    health(args.url)
    run = run_one(lam, mu1, mu2, seed=args.seed, n_jobs=args.n_jobs,
                  url=args.url, name=name)
    hw = None if run.ci is None else (run.ci[1] - run.ci[0]) / 2
    print(f"rho={args.rho} r={args.r:g} seed={args.seed} n_jobs={args.n_jobs:,}")
    print(f"  mean                = {run.mean:.6f}")
    print(f"  intra-run CI        = {run.ci}"
          + (f"  half-width {hw:.6f} ({100 * hw / run.mean:.2f}%)" if hw else "")
          + "   [diagnostic only]")
    print(f"  samplesAnalyzed     = {run.samples_analyzed:,}"
          f"   discarded = {run.samples_discarded}")
    print(f"  success             = {run.success}   [gates nothing, spec 6.2]")
    print(f"  wallClockSeconds    = {run.wall_clock_seconds:.3f}"
          f"   round-trip = {run.round_trip_seconds:.3f}")
    print(f"  jobs/s (engine)     = {run.samples_analyzed / run.wall_clock_seconds:,.0f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Run the live smoke test**

The service must be up (`podman ps` shows `qsim`). This is a real network call:

```bash
cd experiments/table1_qsim
python qsim_fj.py --rho 0.8 --r 4 --n-jobs 2000000 --cross-check-qopt
```

Expected: `samplesAnalyzed = 2,000,000` exactly, `mean` near `5.058`, `jobs/s` around
180,000. The qopt line prints either `ok: identical to qopt's emitter` or
`skipped: qopt not importable` — both acceptable; a drift AssertionError is not.

- [ ] **Step 3: Verify the unit tests still pass**

Run: `cd experiments/table1_qsim && python -m pytest test_qsim_fj.py -v`
Expected: 14 passed.

- [ ] **Step 4: Commit**

```bash
git add experiments/table1_qsim/qsim_fj.py
git commit -m "Add single-run CLI and qopt cross-check to qsim_fj

Compares as canonical JSON (key order is not part of the contract) and
degrades to a skip when qopt is not importable, keeping the experiment
dependency-free."
```

---

### Task 4: Aggregation, reference loading, and the homogeneous closed form

**Files:**
- Create: `experiments/table1_qsim/validate_table1_qsim.py`
- Test: `experiments/table1_qsim/test_validate.py`

**Interfaces:**
- Consumes: nothing from earlier tasks (pure functions).
- Produces: `MU1`, `RHO_VALUES`, `R_VALUES`, `R_CHECK`, `SEEDS`, `DEFAULT_N_JOBS`, `T_975`,
  `grand_mean_and_ci(means) -> (mean, half_width)`,
  `homogeneous_exact(rho, mu=MU1) -> float`,
  `load_reference(path=None) -> dict[(rho, r), (mean, half_width)]`,
  `agreement(t_a, hw_a, t_b, hw_b) -> (gap, ratio, agree)`.

- [ ] **Step 1: Write the failing test**

```python
# experiments/table1_qsim/test_validate.py
import json
import math

import pytest

import validate_table1_qsim as v


def test_grand_mean_and_ci_matches_textbook_t_interval():
    means = [2.0, 2.1, 1.9, 2.05, 1.95]
    mean, hw = v.grand_mean_and_ci(means)
    assert mean == pytest.approx(2.0)
    # s = sample std (ddof=1); hw = t_{0.975,4} * s / sqrt(5)
    s = math.sqrt(sum((m - 2.0) ** 2 for m in means) / 4)
    assert hw == pytest.approx(2.776445 * s / math.sqrt(5))


def test_grand_mean_and_ci_single_seed_has_no_interval():
    mean, hw = v.grand_mean_and_ci([3.0])
    assert mean == 3.0
    assert math.isnan(hw)


def test_homogeneous_exact_is_nelson_tantawi():
    # T = (12 - rho) / (8 (mu - lam)), mu = 1 so lam = rho
    assert v.homogeneous_exact(0.4) == pytest.approx((12 - 0.4) / (8 * 0.6))
    assert v.homogeneous_exact(0.4) == pytest.approx(2.4166666, abs=1e-6)
    assert v.homogeneous_exact(0.95) == pytest.approx((12 - 0.95) / (8 * 0.05))


def test_load_reference_reads_paper_protocol_cells():
    ref = v.load_reference()
    # The r=1 rho=0.4 cell must agree with the closed form to ~1e-4.
    mean, hw = ref[(0.4, 1)]
    assert mean == pytest.approx(v.homogeneous_exact(0.4), abs=1e-3)
    # All 16 cells present.
    for rho in v.RHO_VALUES:
        for r in list(v.R_VALUES) + [v.R_CHECK]:
            assert (rho, r) in ref, f"missing reference cell rho={rho} r={r}"
    # And the CI is a real cross-seed interval, not a placeholder.
    assert hw > 0


def test_load_reference_rebuilds_ci_from_per_seed_means(tmp_path):
    """The reference CI must be recomputed here, so both arms use one construction."""
    cache = tmp_path / "ref.json"
    means = [5.0, 5.1, 4.9, 5.05, 4.95]
    cache.write_text(json.dumps({
        "0.8,4|20000000|500000|0,1,2,3,4": {"per_seed_means": means}}))
    ref = v.load_reference(str(cache))
    expected_mean, expected_hw = v.grand_mean_and_ci(means)
    assert ref[(0.8, 4)] == (pytest.approx(expected_mean), pytest.approx(expected_hw))


def test_load_reference_missing_file_is_an_error(tmp_path):
    with pytest.raises(FileNotFoundError):
        v.load_reference(str(tmp_path / "nope.json"))


@pytest.mark.parametrize("t_b,expect_agree", [
    (5.00, True),    # inside the combined half-width
    (5.30, False),   # well outside
])
def test_agreement(t_b, expect_agree):
    gap, ratio, agree = v.agreement(5.00, 0.02, t_b, 0.03)
    assert agree is expect_agree
    assert gap == pytest.approx(t_b - 5.00)
    assert ratio == pytest.approx(abs(t_b - 5.00) / 0.05)


def test_agreement_handles_nan_half_width():
    """A single-seed arm has no CI; agreement is then undecidable, not True."""
    gap, ratio, agree = v.agreement(5.0, float("nan"), 5.4, 0.01)
    assert agree is None
    assert math.isnan(ratio)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd experiments/table1_qsim && python -m pytest test_validate.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'validate_table1_qsim'`

- [ ] **Step 3: Write minimal implementation**

```python
# experiments/table1_qsim/validate_table1_qsim.py
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

import json
import math
import os

MU1 = 1.0
RHO_VALUES = [0.4, 0.8, 0.9, 0.95]
R_VALUES = [2, 4, 8]
R_CHECK = 1                      # homogeneous: an exact-closed-form gate (spec 7.3)
SEEDS = [0, 1, 2, 3, 4]
DEFAULT_N_JOBS = 20_000_000      # the paper's own job count (spec 6.3)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REFERENCE_FILE = os.path.join(SCRIPT_DIR, "..", "table1_repro", "table1_results.json")
REFERENCE_KEY = "{rho},{r}|20000000|500000|0,1,2,3,4"

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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd experiments/table1_qsim && python -m pytest test_validate.py -v`
Expected: 9 passed.

- [ ] **Step 5: Commit**

```bash
git add experiments/table1_qsim/validate_table1_qsim.py \
        experiments/table1_qsim/test_validate.py
git commit -m "Add aggregation, reference loading, and r=1 closed-form oracle

The reference CI is rebuilt from table1_repro's cached per-seed means so
both arms of the comparison use one interval construction. Adds the
Nelson-Tantawi exact form that gates qsim itself at r=1."
```

---

### Task 5: Resumable run cache

**Files:**
- Modify: `experiments/table1_qsim/validate_table1_qsim.py`
- Test: `experiments/table1_qsim/test_validate.py`

**Interfaces:**
- Consumes: Task 4's module-level constants.
- Produces: `run_key(rho, r, n_jobs, seed) -> str`, `RunCache` class with
  `__init__(path)`, `get(key) -> dict | None`, `put(key, record) -> None`, `path` attribute.

- [ ] **Step 1: Write the failing test**

Append to `test_validate.py`:

```python
def test_run_key_includes_protocol_parameters():
    """A tighter/longer run must not collide with a looser one (spec 6.5)."""
    a = v.run_key(0.8, 4, 20_000_000, 3)
    b = v.run_key(0.8, 4, 40_000_000, 3)
    c = v.run_key(0.8, 4, 20_000_000, 4)
    assert a != b and a != c


def test_cache_round_trip(tmp_path):
    path = tmp_path / "c.json"
    cache = v.RunCache(str(path))
    assert cache.get("k") is None
    cache.put("k", {"mean": 5.0})
    assert cache.get("k") == {"mean": 5.0}
    # Persisted immediately, so an interrupted grid loses at most one run.
    assert json.loads(path.read_text())["k"] == {"mean": 5.0}


def test_cache_reloads_from_disk(tmp_path):
    path = tmp_path / "c.json"
    v.RunCache(str(path)).put("k", {"mean": 1.0})
    assert v.RunCache(str(path)).get("k") == {"mean": 1.0}


def test_cache_writes_atomically(tmp_path):
    """No .tmp litter left behind, and the file is always valid JSON."""
    path = tmp_path / "c.json"
    cache = v.RunCache(str(path))
    for i in range(5):
        cache.put(f"k{i}", {"mean": float(i)})
    assert json.loads(path.read_text())["k4"] == {"mean": 4.0}
    assert list(tmp_path.iterdir()) == [path]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd experiments/table1_qsim && python -m pytest test_validate.py -k cache -v`
Expected: FAIL — `AttributeError: module 'validate_table1_qsim' has no attribute 'RunCache'`

- [ ] **Step 3: Write minimal implementation**

Add to `validate_table1_qsim.py`:

```python
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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd experiments/table1_qsim && python -m pytest test_validate.py -v`
Expected: 13 passed.

- [ ] **Step 5: Commit**

```bash
git add experiments/table1_qsim/validate_table1_qsim.py \
        experiments/table1_qsim/test_validate.py
git commit -m "Add resumable per-run cache with atomic writes

Flushed after every run, not every cell: a multi-hour grid must survive
interruption losing at most one run."
```

---

### Task 6: Gates, grid walk, and console report

**Files:**
- Modify: `experiments/table1_qsim/validate_table1_qsim.py`

**Interfaces:**
- Consumes: `qsim_fj.health`, `qsim_fj.run_one`, `qsim_fj.QsimError` (Tasks 2-3);
  `grand_mean_and_ci`, `homogeneous_exact`, `load_reference`, `agreement`, `RunCache`,
  `run_key` (Tasks 4-5).
- Produces: `CI_NATURE_ROWS`, `check_determinism(url, ...) -> str`,
  `simulate_cell(rho, r, n_jobs, seeds, cache, url) -> dict`,
  `gate_homogeneous(rows) -> list[str]`, `print_report(rows, failures, gate_notes)`,
  `main(argv=None) -> int`.

- [ ] **Step 1: Write the implementation**

Add to `validate_table1_qsim.py` (put `import argparse`, `import sys`, `import time`, and
`import qsim_fj` at the top with the existing imports):

```python
CI_NATURE_ROWS = [
    ("Unit of observation", "5 whole-run means", "batch means within ONE run"),
    ("Method", "Student-t on 5 values",
     "spectral: batch means + polyfit to log-spectrum"),
    ("Autocorrelation", "independent replications",
     "spectral variance-of-the-mean estimate"),
    ("Warm-up", "500 K jobs, caller-specified",
     "none observed (samplesDiscarded = 0)"),   # RETRACTED: false at 20 M; see spec 6.5
    ("Count", "1 per cell", "1 per seed"),
]
"""Why qsim's per-run CI is not comparable to the simple simulator's (spec 6.1).

Reported everywhere the numbers are, so a reader cannot pick up a half-width
without also picking up what produced it.
"""


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
    for seed in seeds:
        key = run_key(rho, r, n_jobs, seed)
        record = cache.get(key)
        if record is None:
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
    }


GATE_RELATIVE_FLOOR = 0.01
"""Relative tolerance floor added to the r=1 gate's CI-based tolerance.

Measured justification: the REFERENCE simulator -- which discards 500 K warm-up
jobs -- itself deviates from the exact closed form by +0.01% at rho=0.4 and +0.47%
at rho=0.95, the latter consuming 39% of its own CI. qsim discards no warm-up at
all (spec 6.5), so its finite-run bias at high load is likely larger. Without a
floor, this gate would fail on that bias and skip all 12 heterogeneous cells.

The gate's job is catching a BROKEN fork-join implementation -- which would be
wrong by tens of percent -- not measuring precision, exactly as with the 2%
semantics oracle in spec 7.4. The measured gap is printed either way, so a real
drift stays visible even when the gate passes.
"""


def gate_homogeneous(rows):
    """The r=1 cells must reproduce the exact closed form (spec 7.3)."""
    notes, failures = [], []
    for row in rows:
        if row["r"] != R_CHECK:
            continue
        exact = homogeneous_exact(row["rho"])
        # hw_a is the relative floor, not 0: see GATE_RELATIVE_FLOOR.
        gap, ratio, ok = agreement(exact, GATE_RELATIVE_FLOOR * exact,
                                   row["t_qsim"], row["hw_qsim"])
        verdict = "ok" if ok else "FAILED"
        notes.append(f"  rho={row['rho']:<5} exact={exact:9.4f} "
                     f"qsim={row['t_qsim']:9.4f} +/-{row['hw_qsim']:.4f}  "
                     f"gap={gap:+.4f} ({ratio:.2f} half-widths)  {verdict}")
        if not ok:
            failures.append(f"r=1 rho={row['rho']}: qsim {row['t_qsim']:.4f} vs exact "
                            f"{exact:.4f}, {ratio:.2f} half-widths apart")
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
    disagreements = []
    for row in rows:
        ref = reference.get((row["rho"], row["r"]))
        if ref is None:
            print(f"{row['rho']:>5} {row['r']:>2} | {'-- no reference cell --':>60}")
            continue
        t_ref, hw_ref = ref
        gap, ratio, agree = agreement(t_ref, hw_ref, row["t_qsim"], row["hw_qsim"])
        mark = {True: "yes", False: "NO", None: "?"}[agree]
        if agree is False:
            disagreements.append((row, t_ref, hw_ref, gap, ratio))
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
    else:
        print("\nALL CELLS AGREE: every qsim mean is within the combined 95% CI of "
              "the simple simulator's.")

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
    return disagreements


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
    args = parser.parse_args(argv)

    n_jobs = args.n_jobs or (200_000 if args.quick else DEFAULT_N_JOBS)
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
        print("\nGATE FAILED: qsim does not reproduce the exact homogeneous result, "
              "so the heterogeneous cells were not run.")
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

    disagreements = print_report(rows, reference, failures, gate_notes)
    write_outputs(rows, reference, make_figure=not args.no_fig)
    # Disagreement is a finding (exit 0); only an untrustworthy run is a failure.
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Add a placeholder for the Task 7 output writer**

`main` calls `write_outputs`, which Task 7 implements. Add this stub now so `--quick`
runs end-to-end, and replace it in Task 7:

```python
def write_outputs(rows, reference, make_figure=True):
    """Replaced in Task 7 by the LaTeX + figure writer."""
    print("\n[write_outputs not implemented yet]")
```

- [ ] **Step 3: Run the smoke test end-to-end**

Service must be up. 16 cells x 5 seeds x 200 K jobs, roughly a minute:

```bash
cd experiments/table1_qsim
python validate_table1_qsim.py --quick --no-fig; echo "exit=$?"
```

Expected: the health and determinism gates pass; the `r=1` gate prints four rows (at
200 K jobs the CI is wide, so it should pass easily); the 16-row table prints. Cells may
legitimately disagree at 200 K jobs — that is fine and must still `exit=0`.

- [ ] **Step 4: Verify unit tests still pass**

Run: `cd experiments/table1_qsim && python -m pytest -v`
Expected: 27 passed (14 + 13).

- [ ] **Step 5: Commit**

```bash
git add experiments/table1_qsim/validate_table1_qsim.py
git commit -m "Add gates, grid walk, and console agreement report

Gate order per spec 7: health, determinism, then the r=1 closed-form
oracle, which stops the run before the heterogeneous cells if qsim
cannot reproduce the exact homogeneous result. Disagreeing cells exit 0
(a finding); only a failed gate or engine error exits non-zero."
```

---

### Task 7: LaTeX and figure output

**Files:**
- Modify: `experiments/table1_qsim/validate_table1_qsim.py`

**Interfaces:**
- Consumes: row dicts from `simulate_cell`, `reference` from `load_reference`,
  `agreement`, `CI_NATURE_ROWS`.
- Produces: `write_outputs(rows, reference, make_figure=True) -> None`, replacing the
  Task 6 stub. Writes `table1_qsim_comparison.tex`, `table1_qsim_results.json`, and
  `table1_qsim_agreement.{png,pdf}`.

- [ ] **Step 1: Replace the stub with the real writer**

```python
def _tex_ci_nature_table():
    """The spec 6.1 table as LaTeX, so it travels with the numbers it qualifies."""
    lines = [r"\begin{tabular}{@{}lll@{}}", r"\toprule",
             r"& Simple simulator & qsim per-run CI \\", r"\midrule"]
    for label, simple, qsim in CI_NATURE_ROWS:
        esc = lambda s: s.replace("%", r"\%").replace("_", r"\_")
        lines.append(f"{esc(label)} & {esc(simple)} & {esc(qsim)} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return lines


def write_outputs(rows, reference, make_figure=True):
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
                   "ci_nature": [list(r) for r in CI_NATURE_ROWS]}, f, indent=2)
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
        xs, heights, colors = [], [], []
        for i, r in enumerate(all_r):
            row, ref = by_cell.get((rho, r)), reference.get((rho, r))
            if row is None or ref is None:
                continue
            gap, ratio, agree = agreement(ref[0], ref[1], row["t_qsim"], row["hw_qsim"])
            signed = ratio if gap >= 0 else -ratio
            xs.append(i)
            heights.append(signed)
            colors.append("steelblue" if agree else "firebrick")
        ax.axhspan(-1, 1, color="lightgray", alpha=0.6, zorder=0,
                   label="agreement band")
        ax.bar(xs, heights, 0.6, color=colors, zorder=2)
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
```

- [ ] **Step 2: Verify the outputs are produced and well-formed**

```bash
cd experiments/table1_qsim
python validate_table1_qsim.py --quick; echo "exit=$?"
python -c "
import json; d=json.load(open('table1_qsim_results.json'))
print('rows:', len(d['rows']))
print('cells with a verdict:', sum(1 for r in d['rows'] if 'agree' in r))
assert len(d['ci_nature']) == 5
"
grep -c "checkmark\|textbf{no}" table1_qsim_comparison.tex
ls -la table1_qsim_agreement.png table1_qsim_agreement.pdf
```

Expected: `rows: 16`, 16 cells with a verdict, the `.tex` containing 16 verdict marks,
and both figure files non-empty.

- [ ] **Step 3: Verify unit tests still pass**

Run: `cd experiments/table1_qsim && python -m pytest -v`
Expected: 27 passed.

- [ ] **Step 4: Commit**

```bash
git add experiments/table1_qsim/validate_table1_qsim.py
git commit -m "Add LaTeX and figure output for the agreement comparison

The figure plots gap/combined-half-width against a shaded +/-1 band, so
the grid's verdict is readable at a glance and near-misses stay visible.
The CI-nature table is emitted as LaTeX alongside the comparison table."
```

---

### Task 8: Timing baseline and transient-bias check

**Files:**
- Modify: `experiments/table1_qsim/validate_table1_qsim.py`

**Interfaces:**
- Consumes: `simulate_cell`, `RunCache`, `run_key`, `grand_mean_and_ci`, `agreement`.
- Produces: `time_baseline_cell(rho, r, n_jobs, seeds) -> dict`,
  `run_timing_comparison(rows, n_jobs, seeds, cache, url) -> list[dict]`,
  `run_transient_check(rows, n_jobs, cache, url) -> list[dict]`; CLI flags
  `--time-baseline`, `--transient-check`.

- [ ] **Step 1: Add the timing baseline**

The baseline is re-timed fresh: `table1_repro/run_full.log` is unusable because three
cells logged 958 s / 1013 s / 489 s against ~46 s for identical work (machine contention).

```python
_REPO_ROOT = os.path.join(SCRIPT_DIR, "..", "..")


def time_baseline_cell(rho, r, n_jobs, seeds):
    """Re-time forkjoin.simulate fresh for one cell (spec 10).

    `forkjoin` is imported lazily -- inside the function, not at module scope -- so
    the qsim-only path works without numpy installed. The sys.path entry is added
    once at module level rather than per call, which would grow sys.path per cell.
    """
    if _REPO_ROOT not in sys.path:
        sys.path.insert(0, _REPO_ROOT)
    from forkjoin import simulate

    lam, mu2 = rho * MU1, r * MU1
    started = time.time()
    means = []
    for seed in seeds:
        result = simulate(lam, MU1, mu2, n_jobs=n_jobs, warmup=500_000, seed=seed)
        means.append(result.mean_response_time)
    elapsed = time.time() - started
    mean, hw = grand_mean_and_ci(means)
    return {"rho": rho, "r": r, "seconds": elapsed, "t_simple": mean, "hw_simple": hw,
            "jobs_per_second": n_jobs * len(seeds) / elapsed}


def run_timing_comparison(rows, n_jobs, seeds, cache, url):
    """Interleaved, never concurrent: simple then qsim, cell by cell (spec 10).

    Interleaving matters because slow machine drift then hits both arms equally
    instead of penalising whichever ran second.
    """
    print("\nTIMING (interleaved per cell, never concurrent; equal work at "
          f"{n_jobs:,} jobs/seed)")
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
                 "hw_qsim": fresh["hw_qsim"]}
        results.append(entry)
        print(f"  {rho:>5} {r:>2} | {base['seconds']:>9.1f} "
              f"{fresh['engine_seconds']:>9.1f} {ratio:>7.1f}x | "
              f"{base['jobs_per_second']:>13,.0f} "
              f"{fresh['jobs_per_second']:>12,.0f} | "
              f"{base['hw_simple']:>9.4f} {fresh['hw_qsim']:>9.4f}")
    print("  Caveats: qsim's seconds are engine-side inside a podman VM (see the "
          "2-cell local-JVM probe in the README), and neither simulator removes the "
          "same transient (spec 6.5).")
    return results
```

- [ ] **Step 2: Add the transient-bias check**

qsim has no warm-up knob and `samplesDiscarded` was 0 in every probe, so the bias is
quantified rather than documented away (§6.5).

```python
def run_transient_check(rows, n_jobs, cache, url):
    """Compare N against 2N at one seed: is the un-removed transient visible?

    qsim exposes no warm-up parameter and discards nothing, whereas the paper
    discards 500 K jobs. If the N and 2N means agree inside the N-run's own
    intra-run CI, the transient is below the noise floor (spec 6.5).
    """
    print(f"\nTRANSIENT-BIAS CHECK ({n_jobs:,} vs {2 * n_jobs:,} jobs, seed 0)")
    print(f"  {'rho':>5} {'r':>2} | {'T(N)':>9} {'T(2N)':>9} {'diff':>9} "
          f"{'diff %':>8} {'within CI':>10}")
    results = []
    for row in rows:
        rho, r = row["rho"], row["r"]
        short = simulate_cell(rho, r, n_jobs, [0], cache, url)
        long = simulate_cell(rho, r, 2 * n_jobs, [0], cache, url)
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
    if outside:
        signs = {"low" if x["diff"] > 0 else "high" for x in outside}
        print(f"  {len(outside)} cell(s) shifted beyond the noise floor; the N-job "
              f"mean is biased {'/'.join(sorted(signs))}. Reportable bias, not an error.")
    else:
        print("  All cells: the un-removed transient is below the noise floor.")
    return results
```

- [ ] **Step 3: Wire both into `main`**

Add the two flags to the parser, and call them after `print_report` but before
`write_outputs`:

```python
    parser.add_argument("--time-baseline", action="store_true",
                        help="re-time forkjoin.simulate and compare cost (spec 10)")
    parser.add_argument("--transient-check", action="store_true",
                        help="compare N vs 2N to quantify un-removed warm-up (spec 6.5)")
```

```python
    timing = (run_timing_comparison(rows, n_jobs, args.seeds, cache, args.url)
              if args.time_baseline else None)
    transient = (run_transient_check(rows, n_jobs, cache, args.url)
                 if args.transient_check else None)
    write_outputs(rows, reference, make_figure=not args.no_fig,
                  timing=timing, transient=transient)
```

And extend `write_outputs` to persist them:

```python
def write_outputs(rows, reference, make_figure=True, timing=None, transient=None):
```

then include them in the JSON payload:

```python
        json.dump({"rows": records,
                   "ci_nature": [list(r) for r in CI_NATURE_ROWS],
                   "timing": timing, "transient_check": transient}, f, indent=2)
```

- [ ] **Step 4: Smoke-test both flags**

```bash
cd experiments/table1_qsim
python validate_table1_qsim.py --quick --time-baseline --transient-check --no-fig
echo "exit=$?"
python -c "
import json; d=json.load(open('table1_qsim_results.json'))
assert d['timing'] and len(d['timing'])==16, d['timing']
assert d['transient_check'] and len(d['transient_check'])==16
print('timing and transient sections present')
"
```

Expected: both tables print, `exit=0`, and both JSON sections have 16 entries.

- [ ] **Step 5: Commit**

```bash
git add experiments/table1_qsim/validate_table1_qsim.py
git commit -m "Add interleaved timing comparison and transient-bias check

Timing re-times forkjoin.simulate fresh (run_full.log has 20x contention
outliers) and interleaves the two simulators per cell, never running them
concurrently. Since qsim has no warm-up knob and discards nothing, the
N-vs-2N check quantifies the bias instead of documenting it away."
```

---

### Task 9: README, full run, and committed artifacts

**Files:**
- Create: `experiments/table1_qsim/README.md`
- Generated: `qsim_results.json`, `table1_qsim_results.json`,
  `table1_qsim_comparison.tex`, `table1_qsim_agreement.{png,pdf}`

**Interfaces:**
- Consumes: everything above.
- Produces: the committed validation artifacts.

- [ ] **Step 1: Write the README**

It must contain: what this validates, the parameterization, the protocol *and why it
differs* from the paper's, the §6.1 CI-nature table, the qsim-service findings (§6.2) with
a link to qsim-service#10, how to start the service both ways, the commands, and how to
read the outputs.

```bash
cd experiments/table1_qsim
cat > README.md <<'MARKDOWN'
# Validating Table 1 with qsim-service

Re-measures the 16 simulated cells of **Table 1** (`\label{tab:approx-results}`) with an
**independent** discrete-event simulator — [qsim-service](https://github.com/atantawi/qsim-service),
an HTTP/JSON wrapper around the headless JMT engine — and reports whether the two
simulators agree, cell by cell, within their confidence intervals.

This does **not** replace Table 1's numbers and does not re-derive the approximation
error columns. Design: `docs/superpowers/specs/2026-07-31-qsim-table1-validation-design.md`.

## Parameterization

`mu_1 = 1`, `mu_2 = r`, `lambda = rho * mu_1` — identical to `../table1_repro/`, so
`rho` is bottleneck utilization and `r = mu_2 / mu_1`.

Grid: `rho in {0.4, 0.8, 0.9, 0.95}`, `r in {1, 2, 4, 8}` (`r = 1` is the closed-form gate).

## Protocol, and why it is not the paper's

Per seed: a **fixed-length** run of 20 M jobs (`disableStatisticStop: true`,
`maxSamples: 20000000`) — the paper's own job count. Five seeds (`0..4`), aggregated by
the independent-replications method: grand mean plus `t_{0.975,4} * s / sqrt(5)`.

Two deliberate departures, both forced by measured service behaviour:

1. **No convergence stopping.** Requesting `minSamples: 1e6` with `precision: 0.005`
   produced a run of **10,880 samples** that still reported `success: true` while its
   actual relative half-width was 7.9%. `minSamples` is accepted by the API and never
   reaches JMT — filed as [qsim-service#10](https://github.com/atantawi/qsim-service/issues/10).
   Run length here is verified by `samplesAnalyzed == N` and by nothing else.
2. **No warm-up.** qsim exposes no warm-up knob and `samplesDiscarded` is always 0,
   whereas the paper discards 500 K jobs. Rather than assume this away,
   `--transient-check` compares N against 2N; see the output for the measured verdict.

## The two CIs are not the same kind of object

Both columns in the comparison table are cross-seed Student-t intervals, built by one
code path, so they are comparable. But qsim *also* reports a per-run interval, and that
one is a different animal — recorded as a diagnostic only:

|  | Simple simulator | qsim per-run CI |
|---|---|---|
| Unit of observation | 5 whole-run means | batch means **within one run** |
| Method | Student-t on 5 values | spectral: batch means + polynomial fit to the log-spectrum |
| Autocorrelation | independent replications | spectral variance-of-the-mean estimate |
| Warm-up | 500 K jobs, caller-specified | none observed (`samplesDiscarded = 0`) |
| Count | 1 per cell | 1 per seed |

## Running

Start the service (it is never started or stopped by these scripts):

```bash
# podman/docker -- the authoritative runtime (temurin 17-jre, the tested JVM)
cd ~/Projects/quantum/qsim-service
podman build -t qsim-service:0.1.0 .
podman run -d --rm --name qsim -p 8080:8080 qsim-service:0.1.0

# or a local JVM (note: JMT on JDK > 17 is untested by that project)
mvn package
java -cp "target/qsim-service.jar:target/dependency/*:lib/JMT-singlejar-1.4.0.jar" \
     qsim.http.App
```

Then:

```bash
python qsim_fj.py --rho 0.8 --r 4 --cross-check-qopt   # one run, ~11 s
python validate_table1_qsim.py --quick                 # smoke test, ~1 min
python validate_table1_qsim.py                         # full grid, ~2.4 h
python validate_table1_qsim.py --time-baseline --transient-check   # + ~2.5 h
python -m pytest                                       # unit tests, no network
```

`QSIM_URL` overrides the default `http://localhost:8080`. Every run is cached in
`qsim_results.json` immediately, so an interrupted grid resumes.

**Exit codes:** disagreeing cells exit **0** — that is the experiment's finding, not its
failure. Non-zero means the run itself is untrustworthy: a failed gate, an engine error,
or an unreachable service.

## Outputs

| File | Contents |
|---|---|
| `qsim_results.json` | Per-run cache (mean, samples, timings, per-run CI) |
| `table1_qsim_results.json` | Per-cell comparison, timing, transient check |
| `table1_qsim_comparison.tex` | Agreement table + the CI-nature table above |
| `table1_qsim_agreement.{png,pdf}` | Gap / combined half-width, with a ±1 agreement band |

## Correctness gates

Run in order, each gating the next: `GET /health`; determinism (same seed twice must
give the same mean); and the **`r = 1` closed form** — for equal rates the exact result
is known (Nelson & Tantawi 1988, `T = (12 - rho) / (8 (mu - lambda))`), making those four
cells an analytic oracle for qsim itself. If they fail, the heterogeneous cells are not
run. Additionally, every run asserts `system-response-time == response-time` (this
network has one station, so they must agree; the probe returned them identical to the
last digit, confirming fork-join *region* semantics).
MARKDOWN
```

- [ ] **Step 2: Run the full grid**

Long-running (~2.4 h). Log it, since the log is the record of what actually happened:

```bash
cd experiments/table1_qsim
python validate_table1_qsim.py 2>&1 | tee run_full.log; echo "exit=$?"
```

Expected: gates pass, 16 rows, an explicit agreement verdict. Read the disagreement
block and the gate output before proceeding — if the `r=1` gate failed, stop and
investigate rather than committing artifacts.

- [ ] **Step 3: Run the timing and transient passes**

```bash
cd experiments/table1_qsim
python validate_table1_qsim.py --time-baseline --transient-check \
    2>&1 | tee run_timing.log; echo "exit=$?"
```

The qsim cells are cached from Step 2, so this pays only for the baseline re-timing and
the 2N runs.

- [ ] **Step 4: Quantify the podman VM timing overhead (2 cells)**

Spec §9: the grid's timing is measured inside a Linux VM, so it is not natively
comparable. Start a local-JVM instance on a different port and re-run two cells:

```bash
cd ~/Projects/quantum/qsim-service
QSIM_PORT=8081 java -cp "target/qsim-service.jar:target/dependency/*:lib/JMT-singlejar-1.4.0.jar" \
    qsim.http.App &
cd -
for cell in "0.8 4" "0.95 8"; do
  set -- $cell
  echo "=== rho=$1 r=$2 ==="
  python qsim_fj.py --rho "$1" --r "$2" --n-jobs 2000000 --url http://localhost:8081
done
kill %1
```

Record both `jobs/s` figures in the README under a "Timing caveat" heading, stating the
podman-to-local ratio explicitly.

- [ ] **Step 5: Commit artifacts and the README**

```bash
cd /Users/tantawi/Projects/fork-join
git add experiments/table1_qsim/README.md \
        experiments/table1_qsim/qsim_results.json \
        experiments/table1_qsim/table1_qsim_results.json \
        experiments/table1_qsim/table1_qsim_comparison.tex \
        experiments/table1_qsim/table1_qsim_agreement.png \
        experiments/table1_qsim/table1_qsim_agreement.pdf \
        experiments/table1_qsim/run_full.log \
        experiments/table1_qsim/run_timing.log
git commit -m "Add qsim-service validation results for Table 1

Full grid: 16 cells x 5 seeds x 20M fixed-length jobs, cross-seed t-CIs
on both arms. Includes the interleaved timing comparison, the N-vs-2N
transient-bias check, and the podman-vs-local-JVM timing caveat."
```

- [ ] **Step 6: Open the issue and the PR**

```bash
cd /Users/tantawi/Projects/fork-join
gh issue create --title "Validate Table 1's simulation numbers with an independent simulator" \
  --body "Table 1's simulated sojourn times come from this repo's Lindley-recursion
simulator. This tracks re-measuring the same 16 cells with qsim-service (headless JMT
over HTTP) and reporting cell-by-cell CI overlap.

Design: \`docs/superpowers/specs/2026-07-31-qsim-table1-validation-design.md\`

Constrained by https://github.com/atantawi/qsim-service/issues/10: \`minSamples\` never
reaches JMT and \`success\` does not report convergence, so only fixed-length runs are
used and run length is verified via \`samplesAnalyzed\`."

git push -u origin experiments/qsim-table1-validation
gh pr create --fill --body "$(cat <<'BODY'
Validates Table 1's simulated sojourn times against an independent discrete-event
simulator (qsim-service / headless JMT), per
`docs/superpowers/specs/2026-07-31-qsim-table1-validation-design.md`.

Contains the design spec, this plan, `experiments/table1_qsim/`, and the result
artifacts. Nothing in `forkjoin/` changes and `pyproject.toml` gains no dependency —
the client is stdlib `urllib`.

See `experiments/table1_qsim/README.md` for the agreement verdict, the protocol
departures forced by atantawi/qsim-service#10, and why the two simulators' per-run
CIs are not comparable objects.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
BODY
)"
```

---

## Self-Review

**Spec coverage.** §1-3 → Task 4 (grid constants, reference, agreement) + Task 9 README.
§4 layout → all tasks; §4.1 split → Tasks 1-3 vs 4-8; §4.2 dependency stance → Task 1
(stdlib) + Task 3 (`cross_check_qopt`). §5 model + closed measure list → Task 1. §6.1
CI-nature table → Task 6 (`CI_NATURE_ROWS`), Task 7 (LaTeX), Task 9 (README) — all three
required outlets. §6.2 findings → Task 1 (no `minSamples`/`precision`), Task 2
(`QsimRunLengthError`, `success` ignored), Task 9 README. §6.3 fixed-N protocol → Tasks 1,
6. §6.4 aggregation → Task 4. §6.5 transient → Task 8. §7 gates 1-3 → Task 6
(`health`, `check_determinism`, `gate_homogeneous`); gate 4 → Task 2 (`QsimSemanticsError`).
§8 error table → Task 2 (status mapping, run-length, `completed`) + Task 6 (per-cell
`try/except` continuing the grid); exit codes → Task 6 `main`; cache/resumability →
Task 5. §9 runtime → Task 9 Steps 2-4. §10 timing → Task 8. §11 testing → Tasks 1-5
(`expected_model.json`, `--cross-check-qopt`, `--quick`). §12 deliverables 1-7 → Task 9.

**Placeholder scan.** One intentional stub: Task 6 Step 2's `write_outputs`, replaced in
Task 7 Step 1 — flagged in both places. No TBD/TODO, no "add error handling", no "similar
to Task N"; every code step carries runnable code.

**Type consistency.** `FJRun` fields (Task 2) are consumed by name in `simulate_cell`
(Task 6): `mean, ci, samples_analyzed, samples_discarded, success, wall_clock_seconds,
round_trip_seconds, system_mean` — all defined. `agreement()` returns
`(gap, ratio, agree)` with `agree in {True, False, None}`, and every consumer (Task 6
`print_report`, Task 7 `.tex` and figure, Task 8 transient) handles `None`. Row dict keys
written by `simulate_cell` (`t_qsim, hw_qsim, engine_seconds, jobs_per_second,
intra_run_half_widths`) are exactly those read in Tasks 6-8. `grand_mean_and_ci` returns
NaN for k=1, which `agreement` maps to `agree=None` — consistent with Task 8's
single-seed transient runs.

**Fixes applied during review.** (1) `write_outputs` gained `timing=`/`transient=`
keyword arguments in Task 8 Step 3, since Task 7 defined it with three parameters — the
signature change is now explicit rather than implied. (2) Task 4's
`test_agreement_handles_nan_half_width` was added after noticing `simulate_cell` with one
seed yields a NaN half-width, which an earlier draft of `agreement` would have reported
as agreement. (3) Task 6's `gate_homogeneous` passes `hw_a=0.0` for the exact value rather
than NaN, so the gate produces a real verdict instead of `None`. (4) Amended during
execution: that half-width is now `GATE_RELATIVE_FLOOR * exact` rather than `0.0`. Task 4
revealed that the reference simulator itself deviates from the exact form by +0.47% at
rho=0.95 -- 39% of its own CI -- so a CI-only tolerance risked failing the gate on
finite-run bias and skipping all 12 heterogeneous cells. Human partner approved the 1%
relative floor.
