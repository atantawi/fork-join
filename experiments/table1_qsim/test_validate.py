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


# --- regression tests for the two defects found during in-session review ---

def test_gate_homogeneous_fails_when_a_cell_is_missing():
    """A gate that inspects zero cells is not a passing gate (spec 7.3).

    Regression test for a real defect: when an r=1 cell's runs raised, the cell was
    absent from `rows`, gate_homogeneous iterated nothing, reported no failures, and
    the 12 heterogeneous cells ran anyway despite an unsatisfied oracle.
    """
    notes, failures = v.gate_homogeneous([])
    assert len(failures) == len(v.RHO_VALUES)        # one per missing rho, not zero
    assert notes == []

    good = [{"rho": rho, "r": 1, "t_qsim": v.homogeneous_exact(rho), "hw_qsim": 0.01}
            for rho in v.RHO_VALUES]
    notes, failures = v.gate_homogeneous(good)
    assert failures == []
    assert len(notes) == len(v.RHO_VALUES)

    notes, failures = v.gate_homogeneous(good[:-1])  # one cell short
    assert len(failures) == 1
    assert str(v.RHO_VALUES[-1]) in failures[0]


def test_gate_homogeneous_uses_the_relative_floor():
    """The floor must widen the tolerance, not be ignored (measured: 0.968 -> 0.40)."""
    rho = 0.95
    exact = v.homogeneous_exact(rho)
    # A deviation just outside the CI alone but inside CI + 1% of exact.
    hw = 0.01 * exact * 0.5
    row = {"rho": rho, "r": 1, "t_qsim": exact + hw * 1.5, "hw_qsim": hw}
    notes, failures = v.gate_homogeneous([row] + [
        {"rho": r, "r": 1, "t_qsim": v.homogeneous_exact(r), "hw_qsim": 0.01}
        for r in v.RHO_VALUES if r != rho])
    assert failures == [], "the 1% floor should admit this deviation"
    # Without any floor the same deviation is 1.5 half-widths out, i.e. a failure.
    _, _, agree_ci_only = v.agreement(exact, 0.0, row["t_qsim"], hw)
    assert agree_ci_only is False


def test_time_baseline_normalises_by_warmup_plus_n_jobs():
    """forkjoin.simulate runs warmup + n_jobs, so throughput must use the true total.

    Regression test: the original formula divided by n_jobs alone, understating the
    reference simulator's rate by 2.5% at 20M jobs and 3.5x at --quick's 200K.
    """
    n_jobs = 1_000
    result = v.time_baseline_cell(0.8, 4, n_jobs, [0])
    total = v.BASELINE_WARMUP + n_jobs
    assert result["total_jobs_per_seed"] == total
    assert result["warmup"] == v.BASELINE_WARMUP
    # jobs_per_second must be built from the total, not from n_jobs.
    assert result["jobs_per_second"] == pytest.approx(total / result["seconds"])
    assert result["jobs_per_second"] > n_jobs / result["seconds"]


def _row(rho, r, t_qsim, hw):
    return {"rho": rho, "r": r, "t_qsim": t_qsim, "hw_qsim": hw,
            "jobs_per_second": 1.0}


def test_verdict_requires_every_cell_compared_and_decided(capsys):
    """"ALL CELLS AGREE" must not be printed over data that cannot support it.

    Regression test for a real defect: `disagreements` stayed empty when cells were
    never compared (no reference entry) or were undecidable (NaN half-width), so the
    script claimed whole-grid agreement at exit 0. With an empty reference dict it
    said so over ZERO compared cells.
    """
    rows = [_row(0.8, 4, 5.0, 0.01)]

    # (a) no reference entry for the cell -> uncompared, so no verdict.
    _, complete = v.print_report(rows, {}, [], [])
    assert complete is False
    assert "NO VERDICT" in capsys.readouterr().out

    # (b) undecidable: single-seed NaN half-width, and a wildly wrong value.
    _, complete = v.print_report([_row(0.8, 4, 99.0, float("nan"))],
                                 {(0.8, 4): (5.0, 0.01)}, [], [])
    assert complete is False
    out = capsys.readouterr().out
    assert "NO VERDICT" in out and "ALL CELLS AGREE" not in out

    # (c) a genuine agreement still gets the positive verdict.
    _, complete = v.print_report(rows, {(0.8, 4): (5.0, 0.01)}, [], [])
    assert complete is True
    assert "ALL CELLS AGREE" in capsys.readouterr().out

    # (d) a real disagreement is a verdict too, and stays exit-0 material.
    _, complete = v.print_report(rows, {(0.8, 4): (9.0, 0.01)}, [], [])
    assert complete is True
    assert "DISAGREEING CELLS" in capsys.readouterr().out


def test_verdict_withheld_when_the_run_has_failures(capsys):
    """A partial run must not print a positive whole-grid claim (gate-abort path)."""
    _, complete = v.print_report([_row(0.8, 4, 5.0, 0.01)], {(0.8, 4): (5.0, 0.01)},
                                 ["r=1 rho=0.95: gate failed"], [])
    assert complete is False
    out = capsys.readouterr().out
    assert "NO VERDICT" in out and "ALL CELLS AGREE" not in out


# --- regression tests for the two defects found in PR review ---

def test_gate_homogeneous_separates_undecidable_from_a_failed_oracle():
    """A cell with no CI is UNDECIDABLE, not a missed closed form.

    Regression test for a real defect: `agreement` returns tristate `ok`, but the gate
    used bare truthiness, so `None` (NaN half-width, i.e. a single seed) took the same
    branch as `False`. `--seeds 0` therefore reported an EXACT match -- gap +0.0000 --
    as "qsim does not reproduce the exact homogeneous result ... nan tolerance units
    apart", and skipped all 12 heterogeneous cells on that false diagnosis.
    """
    exact_rows = [{"rho": rho, "r": 1, "t_qsim": v.homogeneous_exact(rho),
                   "hw_qsim": float("nan")} for rho in v.RHO_VALUES]
    notes, failures = v.gate_homogeneous(exact_rows)

    # Still blocks the heterogeneous cells: an undecidable oracle is not a satisfied
    # one, and spec 7.3 requires satisfaction. Fail-closed is deliberate.
    assert len(failures) == len(v.RHO_VALUES)
    # But it must not claim qsim missed the closed form, nor print "nan" as a number.
    for note, failure in zip(notes, failures):
        assert "UNDECIDABLE" in note
        assert "FAILED" not in note
        assert "nan tolerance units" not in note
        assert "UNDECIDABLE rather than failed" in failure
        assert "tolerance units apart" not in failure

    # A genuine miss is still reported as a failure, with a real ratio. Give every rho
    # so the missing-cell failures (emitted first) do not crowd the list.
    bad = [{"rho": rho, "r": 1, "hw_qsim": 0.001,
            "t_qsim": v.homogeneous_exact(rho) * (1.5 if rho == 0.4 else 1.0)}
           for rho in v.RHO_VALUES]
    notes, failures = v.gate_homogeneous(bad)
    assert "FAILED" in notes[0] and "UNDECIDABLE" not in notes[0]
    assert len(failures) == 1 and "tolerance units apart" in failures[0]


@pytest.mark.parametrize("n_jobs,expect_warning", [
    (v.REFERENCE_N_JOBS, False),   # the protocol length: nothing to warn about
    (200_000, True),               # --quick's length, reached via --n-jobs too
    (2_000_000, True),             # any other short run
    (40_000_000, True),            # longer is also a mismatch, not a free pass
])
def test_short_runs_warn_on_length_not_on_the_quick_flag(n_jobs, expect_warning):
    """The trap is the run length, so the guard must key on it (not on --quick).

    Regression test: the warning was gated on `args.quick`, so `--n-jobs 200000` took
    the identical 200K-vs-20M comparison and printed "ALL CELLS AGREE" at exit 0 with
    no warning at all. REFERENCE_KEY pins the reference arm to 20M jobs.
    """
    warning = v.protocol_mismatch_warning(n_jobs)
    assert (warning is not None) is expect_warning
    if expect_warning:
        assert f"{n_jobs:,}" in warning and f"{v.REFERENCE_N_JOBS:,}" in warning
        assert "NEVER evidence" in warning


def test_reference_key_agrees_with_the_named_run_length():
    """REFERENCE_N_JOBS must be the length actually encoded in the cache key."""
    assert f"|{v.REFERENCE_N_JOBS}|" in v.REFERENCE_KEY
    assert v.load_reference()[(0.8, 4)][1] > 0      # the real cache still resolves


def test_duplicate_seeds_are_rejected_before_any_simulation(capsys):
    """A repeated seed is not a replication: it shrinks the CI without adding data.

    `--seeds 0 0 0` re-reads one cached run three times, so the per-seed means collapse
    to a single value, s = 0, and the cell reports a ZERO-WIDTH interval as certainty --
    the same species of defect as claiming agreement over undecidable cells. Rejected
    at argument-parse time, so no simulation is wasted first.
    """
    # The failure mode this guard exists to prevent, shown on the aggregation directly.
    _, hw = v.grand_mean_and_ci([7.0, 7.0, 7.0])
    assert hw == 0.0

    with pytest.raises(SystemExit) as exc:
        v.main(["--seeds", "0", "0", "0"])
    assert exc.value.code == 2
    assert "repeated seed(s) [0]" in capsys.readouterr().err

    with pytest.raises(SystemExit):
        v.main(["--seeds", "0", "0", "1", "1"])
    assert "repeated seed(s) [0, 1]" in capsys.readouterr().err

    # Ordering: with 11 repeats the message must blame the repetition, not the arity.
    with pytest.raises(SystemExit):
        v.main(["--seeds"] + ["0"] * 11)
    err = capsys.readouterr().err
    assert "repeated seed" in err and "t table covers" not in err


def test_grand_mean_and_ci_rejects_more_replications_than_the_t_table_covers():
    """Past 10 seeds there is no quantile; that must be a clear error, not KeyError."""
    with pytest.raises(ValueError, match="degrees of freedom"):
        v.grand_mean_and_ci([1.0] * (max(v.T_975) + 2))
    # The largest supported count still works.
    mean, hw = v.grand_mean_and_ci([1.0, 2.0] * ((max(v.T_975) + 1) // 2))
    assert hw > 0
