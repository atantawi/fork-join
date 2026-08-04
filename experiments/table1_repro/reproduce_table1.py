#!/usr/bin/env python3
"""Reproduce Table 1 of the Performance2026 paper (Squillante & Tantawi).

Paper: "...quantum-centric supercomputing..." (Performance2026 submission),
Table 1 (\\label{tab:approx-results}, p.17): the expected sojourn time of the
core parallel computational phase from the three interpolation approximations,
compared against discrete-event simulation, as a function of system load
rho = gamma/mu^* and heterogeneity ratio r = mu2/mu1 > 1.

Paper <-> repo notation
-----------------------
  E[T_FJ]      (eq. 6,  \\eqref{eq:thm:bounds-interpolation})  <->  mean_response_time            (T_UL)
  E[T_FJ^(0)]  (eq. 13, \\eqref{eq:...-approximation0})        <->  mean_response_time_lh         (T_LH,  beta=10)
  E[T_FJ^(1)]  (eq. 14, \\eqref{eq:...-approximation1})        <->  mean_response_time_lh_enhanced(T_LHe, beta=10)

Parameterization (verified against the paper's own numbers): with mu^* = mu1 the
slower/bottleneck server, set mu1 = 1, mu2 = r, and gamma = lambda = rho * mu1.
Then rho = gamma/mu^* is the bottleneck utilization and r = mu2/mu1.

Simulation protocol (paper Sec. results, line "20-million jobs ... 500-thousand
... five independent replicas ... t-distribution"):
  20 M jobs, 500 K warmup, independent seeds; T_sim is the grand mean over
  seeds and the 95% CI half-width is t_{0.975, k-1} * s / sqrt(k) with s the
  sample std (ddof=1) of the per-seed means. This is the "independent
  replications" method -- the correct basis for a CI given the strong
  autocorrelation of consecutive response times near saturation.

  The paper text says five replicas; the default here is k=10 to tighten the
  interval (the t-quantile drop plus the extra sqrt(2) narrows the half-width
  by ~42% at fixed replica spread). The paper text must then say ten, since
  the generated caption reports the actual k.

We drive the library simulate() once per seed (one continuous Lindley recursion,
warmup discarded) -- the same routine used everywhere else in the repo. We do
NOT reimplement the simulator.

Outputs (all in this directory):
  table1_results.json  -- full numeric results + per-seed means (cache-backed)
  table1.tex           -- LaTeX table in the paper's structure/notation
  table1_errors.png/pdf-- relative-error figure

Usage:
  python reproduce_table1.py             # paper protocol: 20M jobs, 10 seeds (~25 min cold)
  python reproduce_table1.py --quick     # 200K jobs, 10 seeds (smoke test, ~seconds)
"""

import argparse
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from forkjoin import (
    simulate,
    mean_response_time,             # E[T_FJ]      -- UL  (eq. 6)
    mean_response_time_lh,          # E[T_FJ^(0)]  -- LH  (eq. 13)
    mean_response_time_lh_enhanced,  # E[T_FJ^(1)] -- LHe (eq. 14)
)

# ---- Table 1 grid ----------------------------------------------------------
MU1 = 1.0                       # bottleneck rate mu^* = mu_min
RHO_VALUES = [0.4, 0.8, 0.9, 0.95]
R_VALUES = [2, 4, 8]            # r = 1 (homogeneous) is exact; shown as a check
R_CHECK = 1                     # extra homogeneous sanity row (not in main table)
BETA = 10.0

# ---- Simulation protocol (paper) -------------------------------------------
DEFAULT_N_JOBS = 20_000_000
DEFAULT_WARMUP = 500_000
# 10 replicas: the t-quantile drop (2.776 -> 2.262) plus the extra sqrt(2) in
# sqrt(k) narrows the CI half-width by ~42% at fixed replica-to-replica spread.
DEFAULT_SEEDS = list(range(10))

# Student-t 0.975 quantiles by degrees of freedom (avoids a scipy dependency).
T_975 = {1: 12.706205, 2: 4.302653, 3: 3.182446, 4: 2.776445, 5: 2.570582,
         6: 2.446912, 7: 2.364624, 8: 2.306004, 9: 2.262157}

# Published Table 1 values (Performance2026), keyed by (rho, r):
#   (T_sim, ci_halfwidth, T_FJ, T_FJ0, T_FJ1)
# Used only for a side-by-side sanity check in the console report.
PAPER = {
    (0.4, 2): (1.816, 0.001, 1.824, 1.834, 1.825),
    (0.4, 4): (1.704, 0.001, 1.704, 1.717, 1.705),
    (0.4, 8): (1.677, 0.001, 1.676, 1.681, 1.676),
    (0.8, 2): (5.109, 0.003, 5.101, 5.168, 5.151),
    (0.8, 4): (5.042, 0.003, 5.016, 5.050, 5.026),
    (0.8, 8): (5.032, 0.003, 5.003, 5.014, 5.005),
    (0.9, 2): (10.133, 0.006, 10.063, 10.170, 10.150),
    (0.9, 4): (10.094, 0.006, 10.009, 10.050, 10.023),
    (0.9, 8): (10.088, 0.006, 10.002, 10.014, 10.004),
    (0.95, 2): (20.260, 0.012, 20.036, 20.174, 20.153),
    (0.95, 4): (20.239, 0.012, 20.005, 20.050, 20.021),
    (0.95, 8): (20.236, 0.012, 20.001, 20.014, 20.003),
}

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_FILE = os.path.join(SCRIPT_DIR, "table1_results.json")
TEX_FILE = os.path.join(SCRIPT_DIR, "table1.tex")

# The published Table 1 was actually produced with this protocol (verified: it
# reproduces every T_sim and CI in PAPER to the third decimal) -- NOT the
# 20M/5-seed/t-CI protocol described in the paper text. --paper-exact selects it.
PAPER_EXACT_SEEDS = [42]
PAPER_EXACT_N_JOBS = 10_000_000


def paper_eq14_literal(gamma, mu1, mu2):
    """E[T_FJ^(1)] exactly as typeset in eq. 14 of the paper.

    NOTE: this first-order form does NOT reproduce the paper's own Table 1
    (it omits the quadratic / heavy-traffic anchoring that
    mean_response_time_lh_enhanced carries). Computed here only to document
    the discrepancy between eq. 14 as typeset and the table numbers.
    """
    lead = (1 / mu1 + mu1 / mu2**2 - 2 * mu1 / (mu1 + mu2) ** 2
            - 2 * mu1**2 * mu2 / (mu1 + mu2) ** 4)
    const = (mu1**2 + mu1 * mu2 + mu2**2) / (mu1 * mu2 * (mu1 + mu2))
    return gamma / (mu1 - gamma) * lead + const


def grand_mean_and_ci(per_seed_means):
    """Grand mean and 95% t-CI half-width across independent replications."""
    k = len(per_seed_means)
    gm = sum(per_seed_means) / k
    if k == 1:
        return gm, float("nan")
    var = sum((m - gm) ** 2 for m in per_seed_means) / (k - 1)
    s = math.sqrt(var)
    half = T_975[k - 1] * s / math.sqrt(k)
    return gm, half


def seed_key(rho, r, n_jobs, warmup, seed):
    """Cache key for one replication. Per-seed (not per-seed-list) so that
    growing the replica count reuses the replications already run: simulate()
    is deterministic in the seed, so a given (cell, n_jobs, warmup, seed) mean
    is the same number regardless of which run produced it."""
    return f"{rho},{r}|{n_jobs}|{warmup}|seed={seed}"


def migrate_cache(cache):
    """Expand legacy per-seed-list entries into per-seed entries, in place.

    Legacy keys look like "0.4,2|20000000|500000|0,1,2,3,4" and carry a
    per_seed_means list aligned with that seed list; only the first seed's
    normal-approximation half-width was recorded.
    """
    migrated = 0
    for key, entry in list(cache.items()):
        if "|seed=" in key or "per_seed_means" not in entry:
            continue
        cell, n_jobs, warmup, _ = key.split("|")
        rho, r = cell.split(",")
        for i, (s, m) in enumerate(zip(entry["seeds"], entry["per_seed_means"])):
            k = f"{rho},{r}|{n_jobs}|{warmup}|seed={s}"
            if k not in cache:
                cache[k] = {"mean": m,
                            "normal_ci_hw": entry["normal_ci_hw"] if i == 0 else None}
                migrated += 1
        del cache[key]
    if migrated:
        print(f"  [cache] migrated {migrated} legacy replication(s) to per-seed keys")
    return cache


def load_cache():
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE) as f:
            return migrate_cache(json.load(f))
    return {}


def simulate_cell(rho, r, n_jobs, warmup, seeds, cache):
    """Run/load all seeds for one cell. Returns (per_seed_means, normal_ci_hw).

    normal_ci_hw is the single-run normal-approximation half-width from the
    first seed (used only for --paper-exact, to mirror the published table).
    Replications are cached individually, so only seeds not already on disk run.
    """
    lam, mu2 = rho * MU1, r * MU1
    todo = [s for s in seeds if seed_key(rho, r, n_jobs, warmup, s) not in cache]
    if todo:
        print(f"  [run]   rho={rho} r={r} (lam={lam}, mu1={MU1}, mu2={mu2}): "
              f"{len(todo)} of {len(seeds)} seeds {todo} x {n_jobs:,} jobs ...",
              flush=True)
        t0 = time.time()
        for s in todo:
            res = simulate(lam, MU1, mu2, n_jobs=n_jobs, warmup=warmup, seed=s)
            cache[seed_key(rho, r, n_jobs, warmup, s)] = {
                "mean": res.mean_response_time,
                "normal_ci_hw": (res.ci_95[1] - res.ci_95[0]) / 2,
            }
            with open(CACHE_FILE, "w") as f:   # checkpoint after every replica
                json.dump(cache, f, indent=2)
        print(f"          done in {time.time() - t0:.0f}s", flush=True)
    else:
        print(f"  [cache] rho={rho} r={r}")

    entries = [cache[seed_key(rho, r, n_jobs, warmup, s)] for s in seeds]
    means = [e["mean"] for e in entries]
    normal_ci_hw = next((e["normal_ci_hw"] for e in entries
                         if e["normal_ci_hw"] is not None), float("nan"))
    print(f"          per-seed={['%.4f' % m for m in means]}", flush=True)
    return means, normal_ci_hw


def compute_row(rho, r, n_jobs, warmup, seeds, cache, normal_ci=False):
    lam, mu2 = rho * MU1, r * MU1
    means, normal_ci_hw = simulate_cell(rho, r, n_jobs, warmup, seeds, cache)
    if normal_ci:
        # Faithful reproduction of the published table: single-run mean and the
        # normal-approximation CI (ignores autocorrelation; too tight near rho=1).
        t_sim, ci = means[0], normal_ci_hw
    else:
        t_sim, ci = grand_mean_and_ci(means)

    t_fj = mean_response_time(lam, MU1, mu2)                         # UL
    t_fj0 = mean_response_time_lh(lam, MU1, mu2, beta=BETA)          # LH
    t_fj1 = mean_response_time_lh_enhanced(lam, MU1, mu2, beta=BETA)  # LHe
    t_fj1_lit = paper_eq14_literal(lam, MU1, mu2)                    # eq.14 as typeset

    def err(t):
        return (t - t_sim) / t_sim * 100.0

    return {
        "rho": rho, "r": r, "lam": lam, "mu1": MU1, "mu2": mu2,
        "t_sim": t_sim, "ci": ci, "n_seeds": len(seeds),
        "t_fj": t_fj, "err_fj": err(t_fj),
        "t_fj0": t_fj0, "err_fj0": err(t_fj0),
        "t_fj1": t_fj1, "err_fj1": err(t_fj1),
        "t_fj1_literal_eq14": t_fj1_lit, "err_fj1_literal_eq14": err(t_fj1_lit),
    }


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
def print_report(rows, n_jobs, seeds):
    df = len(seeds) - 1
    print("\n" + "=" * 118)
    print(f"REPRODUCED Table 1  ({len(seeds)} seeds x {n_jobs:,} jobs, t-CI df={df})   "
          f"vs  PAPER (Performance2026)")
    print("=" * 118)
    hdr = (f"{'rho':>4} {'r':>2} | {'Tsim(new)':>10} {'+-CI':>7} | "
           f"{'Tsim(pap)':>9} | {'T_FJ':>8} {'err%':>7} | "
           f"{'T_FJ0':>8} {'err%':>7} | {'T_FJ1':>8} {'err%':>7}")
    print(hdr)
    print("-" * 118)
    for r in rows:
        pap = PAPER.get((r["rho"], r["r"]))
        p_ts = f"{pap[0]:>9.3f}" if pap else f"{'--':>9}"
        note = ""
        if pap and not math.isnan(r["ci"]):
            if abs(r["t_sim"] - pap[0]) > max(r["ci"], pap[1]):
                note = "  <-- outside CI"
        print(f"{r['rho']:>4} {r['r']:>2} | {r['t_sim']:>10.4f} {r['ci']:>7.4f} | "
              f"{p_ts} | {r['t_fj']:>8.3f} {r['err_fj']:>+7.2f} | "
              f"{r['t_fj0']:>8.3f} {r['err_fj0']:>+7.2f} | "
              f"{r['t_fj1']:>8.3f} {r['err_fj1']:>+7.2f}{note}")
    print("=" * 118)
    print("T_FJ = UL (eq.6) | T_FJ0 = LH (eq.13) | T_FJ1 = LHe (eq.14, repo quadratic form).")
    print("Note: eq.14 as *typeset* differs from the table's T_FJ1 column; see err_fj1_literal_eq14 in JSON.")


def fmt_err_bold(val, is_best):
    sign = "+" if val >= 0 else ""
    body = f"\\mathbf{{{sign}{abs(val):.2f}}}" if is_best else f"{sign}{val:.2f}"
    # keep sign outside \mathbf for negatives to match paper style (+\mathbf{..} / -\mathbf{..})
    if is_best:
        s = "+" if val >= 0 else "-"
        body = f"{s}\\mathbf{{{abs(val):.2f}}}"
    else:
        body = f"{'+' if val >= 0 else ''}{val:.2f}"
    return f"${body}\\%$"


def generate_tex(rows, n_jobs, seeds, tex_file=TEX_FILE, normal_ci=False,
                 warmup=DEFAULT_WARMUP):
    df = len(seeds) - 1
    ci_desc = (rf"single-run normal-approximation CI (seed {seeds[0]}, "
               rf"{n_jobs:,} jobs)" if normal_ci else
               rf"{len(seeds)} independent seeds x {n_jobs:,} jobs, t-CI (df={df})")
    by = {(r["rho"], r["r"]): r for r in rows}

    # The multi-replica protocol changes what the CI column means (a Student-t
    # interval across independent replica means, rather than a within-run normal
    # approximation), so the caption has to say so explicitly.
    if normal_ci:
        protocol_sentence = ""
    else:
        protocol_sentence = (
            rf" Each simulation entry is the grand mean over {len(seeds)} independent "
            rf"replications of {n_jobs // 1_000_000} million jobs each, with a warm-up "
            rf"period of {warmup // 1000} thousand jobs discarded per replication; the "
            rf"reported $95\%$ confidence interval is the Student-$t$ interval across "
            rf"the {len(seeds)} replica means.")

    lines = [
        r"% Auto-generated by experiments/table1_repro/reproduce_table1.py",
        rf"% Protocol: {ci_desc}, 500K warmup, beta={BETA:g}.",
        (r"% NOTE: the $\ex[T_{FJ}^{(1)}]$ column is the quadratic enhanced-LH form "
         r"(docs/paper eq. 23),"),
        (r"% which retains the heavy-traffic anchoring term $c_2\rho^2$. Eq. 14 must be "
         r"typeset with that"),
        (r"% term for the equation and this column to agree -- see "
         r"eq23_vs_eq14_comparison.md."),
        r"\begin{table}[ht]",
        r"\centering",
        (r"\caption{The expected sojourn time for the core parallel computational "
         r"phase from our three interpolation approximations as a function of the "
         r"system load $\rho=\gamma/\mu^\ast$ and the heterogeneity ratio $r>1$ in "
         r"comparison against simulation, with the lowest error highlighted in bold. "
         r"All three approximations are exact when $r=1$." + protocol_sentence + r"}"),
        r"\label{tab:approx-results}",
        r"\small{",
        r"\begin{tabular}{@{}ccrcccccc@{}}",
        r"\toprule",
        (r"& & \multicolumn{1}{c}{Simulation}"
         r" & \multicolumn{2}{c}{$\ex[T_{FJ}]$ in~\eqref{eq:thm:bounds-interpolation}}"
         r" & \multicolumn{2}{c}{$\ex[T_{FJ}^{(0)}]$ in~\eqref{eq:thm:interpolation-approximation0}}"
         r" & \multicolumn{2}{c}{$\ex[T_{FJ}^{(1)}]$ in~\eqref{eq:thm:interpolation-approximation1}} \\"),
        r"\cmidrule(lr){3-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}\cmidrule(lr){8-9}",
        (r"$\rho$ & $r$ & $\ex[T_\mathrm{sim}] \pm 95\%\,\mathrm{CI}$"
         r" & $\ex[T]$ & Err & $\ex[T]$ & Err & $\ex[T]$ & Err \\"),
        r"\midrule",
    ]
    for i, rho in enumerate(RHO_VALUES):
        if i > 0:
            lines.append(r"\midrule")
        for j, r in enumerate(R_VALUES):
            d = by[(rho, r)]
            errs = {"fj": d["err_fj"], "fj0": d["err_fj0"], "fj1": d["err_fj1"]}
            best = min(errs, key=lambda k: abs(errs[k]))
            rho_cell = rf"\multirow{{{len(R_VALUES)}}}{{*}}{{${rho}$}}" if j == 0 else ""
            ci = d["ci"]
            sim = (f"${d['t_sim']:.3f} \\pm {ci:.3f}$" if not math.isnan(ci)
                   else f"${d['t_sim']:.3f}$")
            row = (f"{rho_cell} & ${r}$ & {sim}"
                   f" & {d['t_fj']:.3f} & {fmt_err_bold(d['err_fj'], best == 'fj')}"
                   f" & {d['t_fj0']:.3f} & {fmt_err_bold(d['err_fj0'], best == 'fj0')}"
                   f" & {d['t_fj1']:.3f} & {fmt_err_bold(d['err_fj1'], best == 'fj1')} \\\\")
            lines.append(row)
    lines += [r"\bottomrule", r"\end{tabular}%", r"}", r"\end{table}"]
    with open(tex_file, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nLaTeX table written to {tex_file}")


def make_figure(rows):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping figure.")
        return
    by = {(r["rho"], r["r"]): r for r in rows}
    fig, axes = plt.subplots(1, len(RHO_VALUES), figsize=(4 * len(RHO_VALUES), 3.6),
                             sharey=True)
    width = 0.25
    x = list(range(len(R_VALUES)))
    labels = [("$T_{FJ}$ (UL)", "err_fj", "steelblue"),
              ("$T_{FJ}^{(0)}$ (LH)", "err_fj0", "darkorange"),
              ("$T_{FJ}^{(1)}$ (LHe)", "err_fj1", "forestgreen")]
    for ax, rho in zip(axes, RHO_VALUES):
        for k, (lab, key, col) in enumerate(labels):
            vals = [by[(rho, r)][key] for r in R_VALUES]
            ax.bar([xi + (k - 1) * width for xi in x], vals, width,
                   label=lab, color=col)
        ax.axhline(0, color="k", linewidth=0.7)
        ax.set_title(rf"$\rho = {rho}$")
        ax.set_xticks(x)
        ax.set_xticklabels([f"$r={r}$" for r in R_VALUES])
        ax.grid(True, axis="y", linestyle="--", linewidth=0.5, color="lightgray")
    axes[0].set_ylabel("Relative error (\\%)")
    axes[0].legend(fontsize=8, loc="best")
    fig.suptitle("Approximation relative error vs. simulation (Table 1 grid)")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = os.path.join(SCRIPT_DIR, f"table1_errors.{ext}")
        fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Figure written to {os.path.join(SCRIPT_DIR, 'table1_errors.{png,pdf}')}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="200K jobs (smoke test)")
    ap.add_argument("--n-jobs", type=int, default=None)
    ap.add_argument("--seeds", type=int, nargs="+", default=DEFAULT_SEEDS)
    ap.add_argument("--paper-exact", action="store_true",
                    help="reproduce the PUBLISHED numbers exactly: single seed 42, "
                         "10M jobs, normal-approximation CI (as Table 1 was made)")
    ap.add_argument("--no-fig", action="store_true", help="skip the figure")
    args = ap.parse_args()

    if args.paper_exact:
        seeds = PAPER_EXACT_SEEDS
        n_jobs = args.n_jobs or PAPER_EXACT_N_JOBS
        warmup = DEFAULT_WARMUP
        normal_ci = True
        tex_file = os.path.join(SCRIPT_DIR, "table1_paperexact.tex")
    else:
        seeds = args.seeds
        n_jobs = args.n_jobs or (200_000 if args.quick else DEFAULT_N_JOBS)
        warmup = 10_000 if args.quick else DEFAULT_WARMUP
        normal_ci = False
        tex_file = TEX_FILE

    print(f"Protocol: {len(seeds)} seed(s) {seeds}, n_jobs={n_jobs:,}, "
          f"warmup={warmup:,}, beta={BETA:g}, "
          f"CI={'normal-approx' if normal_ci else 't-dist across seeds'}")
    cache = load_cache()

    rows = []
    # r=1 homogeneous sanity check (not part of the main table), then r in {2,4,8}.
    for rho in RHO_VALUES:
        rows.append(compute_row(rho, R_CHECK, n_jobs, warmup, seeds, cache, normal_ci))
    main_rows = []
    for rho in RHO_VALUES:
        for r in R_VALUES:
            main_rows.append(compute_row(rho, r, n_jobs, warmup, seeds, cache, normal_ci))

    print_report(main_rows, n_jobs, seeds)
    generate_tex(main_rows, n_jobs, seeds, tex_file=tex_file, normal_ci=normal_ci,
                 warmup=warmup)
    if not args.no_fig:
        make_figure(main_rows)


if __name__ == "__main__":
    main()
