"""Discrete-event simulation of a 2-queue fork-join system."""

from dataclasses import dataclass

import numpy as np


@dataclass
class SimResult:
    mean_response_time: float
    std_response_time: float
    n_samples: int
    ci_95: tuple[float, float]


def simulate(lam, mu1, mu2, n_jobs=1_000_000, warmup=100_000, seed=None):
    """Simulate a 2-queue M/M/1 fork-join system.

    Jobs arrive as Poisson(lam). Each job forks into two tasks with
    Exp(mu1) and Exp(mu2) service times. Response time = max(departure1, departure2) - arrival.
    """
    if lam <= 0 or mu1 <= 0 or mu2 <= 0:
        raise ValueError("All rates must be positive")
    if lam >= mu1 or lam >= mu2:
        raise ValueError(
            f"Stability violated: need lam < min(mu1, mu2), "
            f"got lam={lam}, mu1={mu1}, mu2={mu2}"
        )

    rng = np.random.default_rng(seed)
    total = warmup + n_jobs

    # Generate all random variates at once
    interarrivals = rng.exponential(1 / lam, size=total)
    service1 = rng.exponential(1 / mu1, size=total)
    service2 = rng.exponential(1 / mu2, size=total)

    arrivals = np.cumsum(interarrivals)

    # Track last departure time per server
    depart1 = np.empty(total)
    depart2 = np.empty(total)

    depart1[0] = arrivals[0] + service1[0]
    depart2[0] = arrivals[0] + service2[0]

    for i in range(1, total):
        depart1[i] = max(depart1[i - 1], arrivals[i]) + service1[i]
        depart2[i] = max(depart2[i - 1], arrivals[i]) + service2[i]

    # Response time = max(depart1, depart2) - arrival
    response = np.maximum(depart1[warmup:], depart2[warmup:]) - arrivals[warmup:]

    mean_rt = float(np.mean(response))
    std_rt = float(np.std(response, ddof=1))
    n = len(response)
    margin = 1.96 * std_rt / np.sqrt(n)

    return SimResult(
        mean_response_time=mean_rt,
        std_response_time=std_rt,
        n_samples=n,
        ci_95=(mean_rt - margin, mean_rt + margin),
    )
