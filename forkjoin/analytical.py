"""Closed-form formulas for 2-queue fork-join response time."""

import math


def _validate(lam, mu1, mu2):
    if lam <= 0 or mu1 <= 0 or mu2 <= 0:
        raise ValueError("All rates must be positive")
    if lam >= mu1 or lam >= mu2:
        raise ValueError(
            f"Stability violated: need lam < min(mu1, mu2), "
            f"got lam={lam}, mu1={mu1}, mu2={mu2}"
        )


def upper_bound_independent(lam, mu1, mu2):
    """Independent upper bound T_UB = E[max(R1,R2)] assuming independence."""
    _validate(lam, mu1, mu2)
    return 1 / (mu1 - lam) + 1 / (mu2 - lam) - 1 / (mu1 + mu2 - 2 * lam)


def lower_bound_bottleneck(lam, mu1, mu2):
    """Bottleneck lower bound T_bot = max(T1, T2)."""
    _validate(lam, mu1, mu2)
    return max(1 / (mu1 - lam), 1 / (mu2 - lam))


def upper_bound_split_merge(lam, mu1, mu2):
    """Split-merge upper bound via Pollaczek-Khinchine formula.

    Only valid when lam * E[S_max] < 1 (more restrictive than FJ stability).
    """
    _validate(lam, mu1, mu2)
    es = 1 / mu1 + 1 / mu2 - 1 / (mu1 + mu2)
    es2 = 2 / mu1**2 + 2 / mu2**2 - 2 / (mu1 + mu2) ** 2
    if lam * es >= 1:
        raise ValueError(
            f"Split-merge system unstable: lam*E[S_max] = {lam * es:.4f} >= 1"
        )
    return es + lam * es2 / (2 * (1 - lam * es))


def nelson_tantawi(lam, mu):
    """Exact homogeneous 2-queue result: T = (12 - rho) / (8 * (mu - lam))."""
    if lam <= 0 or mu <= 0:
        raise ValueError("All rates must be positive")
    if lam >= mu:
        raise ValueError(f"Stability violated: need lam < mu, got lam={lam}, mu={mu}")
    rho = lam / mu
    return (12 - rho) / (8 * (mu - lam))


def mean_response_time_lh(lam, mu1, mu2, beta=10.0):
    """Light-heavy traffic interpolation approximation (Reiman-Simon framework).

    Approximates T(rho) as a rational function:

        T_LH = (a0 + a1 * rho) / (mu_min * (1 - rho))

    with rho = lam / mu_min, matching conditions:
      - Light traffic (rho=0): T_LH = T0 = 1/mu1 + 1/mu2 - 1/(mu1+mu2)
      - Heavy traffic (rho->1): mu_min*(1-rho)*T_LH -> h(r)

    where the heavy-traffic factor is:

        h(r) = 1 + (3/8) * r^(-beta),   r = mu_max / mu_min >= 1

    This satisfies h(1) = 11/8 (recovers Nelson-Tantawi exactly for the
    homogeneous case) and h(r) -> 1 as r -> inf (bottleneck M/M/1 behavior
    for highly heterogeneous systems).

    The transition from h=11/8 to h=1 is sharp in practice: simulation data
    shows h ~= 1 already at r=1.5. Large beta (>=10) is needed to match
    observed behavior. The default beta=10 reflects this; it can be refined
    via systematic simulation calibration.

    Args:
        lam: Poisson arrival rate.
        mu1: Service rate of server 1.
        mu2: Service rate of server 2.
        beta: Shape parameter for the heavy-traffic factor h(r). Larger beta
              means faster decay from h=11/8 (homogeneous) to h=1 (heterogeneous).
              Default 10.0 (calibrated to simulation data).
    """
    _validate(lam, mu1, mu2)
    mu_min = min(mu1, mu2)
    mu_max = max(mu1, mu2)
    t0 = 1 / mu1 + 1 / mu2 - 1 / (mu1 + mu2)
    a0 = mu_min * t0
    r = mu_max / mu_min
    h = 1.0 + 0.375 * r ** (-beta)
    a1 = h - a0
    rho = lam / mu_min
    return (a0 + a1 * rho) / (mu_min - lam)


def mean_response_time(lam, mu1, mu2):
    """Approximate mean response time for heterogeneous 2-queue fork-join.

    Uses convex combination: T = (1 - alpha) * T_UB + alpha * T_bot
    where alpha = (rho1 + rho2) / 8.

    Exact for the homogeneous case (mu1 == mu2).
    """
    _validate(lam, mu1, mu2)
    rho1 = lam / mu1
    rho2 = lam / mu2
    alpha = (rho1 + rho2) / 8
    t_ub = upper_bound_independent(lam, mu1, mu2)
    t_bot = lower_bound_bottleneck(lam, mu1, mu2)
    return (1 - alpha) * t_ub + alpha * t_bot
