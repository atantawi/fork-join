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
                    url="http://h:9",
                    transport=_transport(capture=capture, body=_response(n=1000)))
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


def test_run_one_accepts_completed_false_when_full_length():
    """completed:false is a DIAGNOSTIC, not a failure, under disableStatisticStop.

    Measured on qsim-service 0.2.0: a rho=0.9 r=1 run returned completed:false with
    the full 200,000 samples analyzed and wallClockSeconds 1.154 against an 1800 s
    watchdog. It means "CI targets unmet", the same non-event as success:false.
    Failing on it would reject every high-load fixed-length run.
    """
    run = qsim_fj.run_one(0.8, 1.0, 4.0, seed=0, n_jobs=20_000_000,
                          transport=_transport(body=_response(completed=False)))
    assert run.completed is False        # recorded faithfully, not coerced to True
    assert run.mean == pytest.approx(5.0585)


def test_run_one_still_rejects_a_truncated_run_reporting_completed_false():
    """A genuine watchdog kill is caught by the sample count, which is the real gate."""
    truncated = _response(n=7_500, completed=False)
    with pytest.raises(qsim_fj.QsimRunLengthError):
        qsim_fj.run_one(0.8, 1.0, 4.0, seed=0, n_jobs=20_000_000,
                        transport=_transport(body=truncated))


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
