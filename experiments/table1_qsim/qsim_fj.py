#!/usr/bin/env python3
"""One heterogeneous 2-queue fork-join run against qsim-service.

Stdlib only, by design: this experiment must not add a dependency on the
unreleased `qopt` sibling project (spec 4.2). The model JSON below is transcribed
from qopt's own emitters --

    qopt/station.py :: ForkJoinStation.sim_node
    qopt/network.py :: Network.to_model_dict

-- and `python qsim_fj.py --cross-check-qopt` asserts the two still agree, so the
transcription cannot drift silently.

Protocol note: only FIXED-LENGTH runs are used, and run length is verified by
comparing `samplesAnalyzed` against the requested N -- never by `success` or
`completed`, both of which report "CI targets unmet" under disableStatisticStop
and are recorded as diagnostics only.

Two upstream bugs originally forced this and are now FIXED (qsim-service 0.2.0):
`minSamples` never reached JMT (#10, fixed in 395ecfc) and the measure alpha was
written as its complement, inverting every interval and making the stopping rule
pass unconditionally (#12, fixed in df45cd1). Fixed-length runs remain the choice
on their own merits: equal-length replications keep the five per-seed means
i.i.d., which is what makes the cross-seed Student-t interval valid, and 20 M
matches the reference protocol's job count.
"""

import json
import os
import time
import urllib.error
import urllib.request
from dataclasses import dataclass

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
    """5xx, an unreadable body, or a missing required measure."""


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
            f"  podman run -d --name qsim -p 8080:8080 qsim-service:0.2.0\n"
            f"    (needs df45cd1 or later: earlier builds invert the measure alpha)\n"
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
    # `completed` is a DIAGNOSTIC, not a gate. Measured against qsim-service 0.2.0
    # (df45cd1): a rho=0.9 r=1 run returned completed:false with
    # samplesAnalyzed == the full 200,000 requested and wallClockSeconds 1.154
    # against an 1800 s watchdog. Under disableStatisticStop the per-measure CI
    # targets are deliberately not pursued, so `completed:false` means "targets
    # unmet" -- the same non-event as `success:false` -- not "watchdog fired".
    #
    # A genuine watchdog kill is still caught, and caught better, by the
    # samplesAnalyzed check below: a truncated run cannot have run the full N.
    # (On 0.1.0 this field read True spuriously, because the alpha inversion of
    # qsim-service#12 made the stopping rule pass unconditionally.)

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
        completed=bool(response.get("completed", True)),
        wall_clock_seconds=response.get("wallClockSeconds"),
        round_trip_seconds=round_trip,
        system_mean=system_mean,
    )


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
