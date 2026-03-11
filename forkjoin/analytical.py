"""Closed-form formulas for 2-queue fork-join response time."""


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
