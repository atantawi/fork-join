# Traffic-Limits Approximation for Heterogeneous 2-Queue Fork-Join Systems

## 1. Overview

This document describes the design and implementation of the light-heavy traffic interpolation approximation (`mean_response_time_lh`) for the mean response time of a heterogeneous 2-queue fork-join system with Poisson arrivals and exponential service times.

**System:** Two parallel servers with rates $\mu_1 \geq \mu_2$ (i.e., server 1 is faster). Jobs arrive as a Poisson process with rate $\lambda$. Each job forks into two tasks; the job completes when both tasks finish. Stability requires $\lambda < \mu_2 = \mu_{\min}$.

---

## 2. Background and Motivation

### 2.1 Known Exact Results

- **Homogeneous case** ($\mu_1 = \mu_2 = \mu$): Nelson and Tantawi (1988) derived:
  $$T_2^{\text{hom}} = \frac{12 - \rho}{8(\mu - \lambda)}, \quad \rho = \lambda/\mu$$

- **Heterogeneous case**: No closed-form is known. Flatto and Hahn (1984) give an exact analysis via elliptic function parametrization, but it does not yield a simple expression.

### 2.2 Theoretical Bounds

Two simple bounds bracket the true response time:

- **Independent upper bound:** $T_{\text{UB}} = \frac{1}{\mu_1-\lambda} + \frac{1}{\mu_2-\lambda} - \frac{1}{\mu_1+\mu_2-2\lambda}$

- **Bottleneck lower bound:** $T_{\text{bot}} = \frac{1}{\mu_{\min} - \lambda}$

---

## 3. The Reiman-Simon Framework

Following Appendix A of the QC/HPC Resource Allocation paper (Squillante and Tantawi), we approximate $T(\lambda)$ with a rational function of the bottleneck utilization $\rho = \lambda/\mu_{\min}$:

$$\tilde{T}(\rho) = \frac{a_0 + a_1 \rho}{\mu_{\min}(1 - \rho)}$$

The two coefficients are determined by matching two conditions:

| Condition | Expression | Gives |
|-----------|-----------|-------|
| Light traffic: $T(\rho=0) = T_0$ | $T_0 = \frac{1}{\mu_1} + \frac{1}{\mu_2} - \frac{1}{\mu_1+\mu_2}$ | $a_0 = \mu_{\min} \cdot T_0$ |
| Heavy traffic: $\lim_{\rho\to1}\mu_{\min}(1-\rho)T(\rho) = h$ | $h$ depends on heterogeneity | $a_1 = h - a_0$ |

The light-traffic value $T_0 = E[\max(X_1, X_2)]$ for independent $\text{Exp}(\mu_i)$ random variables is exact.

---

## 4. The Heavy-Traffic Factor

The key design question is: **what is the heavy-traffic constant $h$?**

For the 2-queue fork-join system, $h$ depends on the heterogeneity ratio $r = \mu_{\max}/\mu_{\min} \geq 1$:

- **Homogeneous** ($r = 1$): Nelson-Tantawi gives $h = 11/8$. Both servers saturate simultaneously, and their correlated sojourn times produce an excess multiplier of $3/8$ above the bottleneck limit.

- **Highly heterogeneous** ($r \to \infty$): As $\lambda \to \mu_{\min}$, the slower server's sojourn time diverges while the faster server's remains finite. The system behaves like a single M/M/1 bottleneck, so $h \to 1$.

### 4.1 Functional Form

We model the transition with a power-law decay:

$$\boxed{h(r) = 1 + \frac{3}{8} \cdot r^{-\beta}, \quad r = \frac{\mu_{\max}}{\mu_{\min}} \geq 1, \quad \beta > 0}$$

**Properties:**
- $h(1) = 11/8$ for all $\beta$ — recovers Nelson-Tantawi exactly
- $h(r) \to 1$ as $r \to \infty$ for all $\beta > 0$ — correct bottleneck limit
- Monotonically decreasing in $r$
- Preserves the rational-function structure of the Reiman-Simon formula
- Single tunable parameter $\beta$ for simulation calibration

At $\beta = 1$: $h(r) = (8r+3)/(8r)$ is a clean rational form.

### 4.2 Calibration

Simulation data reveals a critically important physical fact: **the transition from $h = 11/8$ to $h \approx 1$ is extremely sharp.** Even at $r = 1.5$ (servers with a 3:2 speed ratio), the effective $h$ is already $\approx 1.008$. This means $\beta$ must be large ($\geq 10$) to match observed behavior.

| $r = \mu_{\max}/\mu_{\min}$ | Simulation-implied $h$ | $h(r)$ at $\beta=10$ |
|:-:|:-:|:-:|
| 1.0 | 1.375 (= 11/8) | 1.375 |
| 1.5 | ≈ 1.008 | 1.025 |
| 2.0 | ≈ 1.000 | 1.003 |
| 3.0 | ≈ 1.000 | 1.000 |

**Default: $\beta = 10$**, based on available simulation data. This can be refined through systematic simulation sweeps (see Section 6).

The default of $\beta = 10$ gives errors within ±2.4% across all tested scenarios, and the homogeneous case is exact by construction.

---

## 5. The Formula

Substituting into the Reiman-Simon form:

$$\boxed{T_{\text{LH}}(\lambda) = \frac{a_0 + (h(r) - a_0)\,\rho}{\mu_{\min} - \lambda}}$$

where:
- $\mu_{\min} = \min(\mu_1, \mu_2)$, $\quad r = \mu_{\max}/\mu_{\min}$
- $\rho = \lambda / \mu_{\min}$
- $a_0 = \mu_{\min} \left(\frac{1}{\mu_1} + \frac{1}{\mu_2} - \frac{1}{\mu_1+\mu_2}\right)$
- $h(r) = 1 + \tfrac{3}{8} r^{-\beta}$

**Verification — homogeneous case** ($\mu_1 = \mu_2 = \mu$, $r = 1$, $\beta$ arbitrary):

$$a_0 = \mu \cdot \frac{3}{2\mu} = \frac{3}{2}, \quad h = \frac{11}{8}, \quad a_1 = \frac{11}{8} - \frac{3}{2} = -\frac{1}{8}$$

$$T_{\text{LH}} = \frac{3/2 - \rho/8}{\mu(1-\rho)} = \frac{12 - \rho}{8(\mu - \lambda)} = T_2^{\text{hom}} \quad \checkmark$$

---

## 6. Numerical Results (β = 10)

Results from discrete-event simulation (2M jobs, 100K warmup) versus both approximations:

| $\mu_1$ | $\mu_2$ | $\lambda$ | $T_{\text{sim}}$ | $T_{\text{LH}}$ | Err% | $T_{\text{UL}}$ | Err% |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| 1.0 | 1.0 | 0.3 | 2.088 | 2.089 | +0.06% | 2.089 | +0.06% |
| 1.0 | 1.0 | 0.6 | 3.564 | 3.562 | −0.04% | 3.562 | −0.04% |
| 1.0 | 1.0 | 0.9 | 13.886 | 13.875 | −0.08% | 13.875 | −0.08% |
| 1.0 | 1.5 | 0.3 | 1.705 | 1.698 | −0.41% | 1.716 | +0.67% |
| 1.0 | 1.5 | 0.6 | 2.767 | 2.776 | +0.34% | 2.799 | +1.16% |
| 1.0 | 1.5 | 0.9 | 10.083 | 10.325 | +2.40% | 10.193 | +1.10% |
| 1.0 | 2.0 | 0.3 | 1.582 | 1.595 | +0.85% | 1.590 | +0.54% |
| 1.0 | 2.0 | 0.6 | 2.623 | 2.667 | +1.69% | 2.641 | +0.68% |
| 1.0 | 2.0 | 0.9 | 9.990 | 10.170 | +1.80% | 10.063 | +0.73% |
| 1.0 | 3.0 | 0.6 | 2.546 | 2.583 | +1.47% | 2.554 | +0.30% |
| 1.0 | 5.0 | 0.6 | 2.514 | 2.533 | +0.77% | 2.517 | +0.13% |

$T_{\text{LH}}$ = light-heavy traffic interpolation (this work); $T_{\text{UL}}$ = upper-lower bound interpolation.

Both methods are now exact for the homogeneous case. $T_{\text{UL}}$ is tighter at high heterogeneity; $T_{\text{LH}}$ has stronger theoretical grounding via the Reiman-Simon framework.

---

## 7. Issues Resolved

### 7.1 Prior Code Issue: Ad-Hoc Exponential Step Function

The previous implementation used:
```python
gamma = 1 + 0.375 * math.exp((1 - alpha) * 100)
```
This produces a near-step-function that is $11/8$ only at exact homogeneity and immediately collapses to $1$ for any heterogeneity (e.g., $h \approx 1.003$ at $r = 1.01$). While numerically reasonable, it had no theoretical motivation and was completely inconsistent with the documentation.

### 7.2 Prior Documentation Issue

The documentation described the simple closed-form $T_{\text{LH}} = 1/(\mu_{\min}-\lambda) + 1/\mu_{\max} - 1/(\mu_1+\mu_2)$, which corresponds to $h = 1$ for all $r$. This formula underestimates the homogeneous case by up to 24% at $\rho = 0.9$.

### 7.3 Resolution

The new implementation uses $h(r) = 1 + \tfrac{3}{8} r^{-\beta}$ with $\beta = 10$:
- Grounded in the Reiman-Simon framework
- Code and documentation are consistent
- Exact for the homogeneous case by construction
- Single interpretable parameter $\beta$ with clear physical meaning

---

## 8. Implementation

**File:** `forkjoin/analytical.py`, function `mean_response_time_lh(lam, mu1, mu2, beta=10.0)`.

```python
mu_min = min(mu1, mu2)
mu_max = max(mu1, mu2)
t0 = 1/mu1 + 1/mu2 - 1/(mu1 + mu2)   # light-traffic value
a0 = mu_min * t0
r  = mu_max / mu_min
h  = 1.0 + 0.375 * r**(-beta)          # heavy-traffic factor
a1 = h - a0
rho = lam / mu_min
return (a0 + a1 * rho) / (mu_min - lam)
```

---

## 9. Future Work

### 9.1 Systematic β Calibration

A sweep over a denser simulation grid (varying both $r$ and $\rho$) would provide a more precise calibrated value for $\beta$. The current $\beta = 10$ is based on the limited set of 8 heterogeneous scenarios in the table above.

### 9.2 First-Derivative Matching (Quadratic Numerator)

The Reiman-Simon framework allows a quadratic numerator $(a_0 + a_1\rho + a_2\rho^2)$ by additionally matching the first light-traffic derivative $f'(0)$, computed via:
$$\frac{df}{d\rho}\bigg|_{\rho=0} = \mu \int_{-\infty}^{+\infty} \left[\hat{f}(\{t\}) - \hat{f}(\emptyset)\right] dt$$
This requires computing the conditional expected response time given one prior arrival at time $t$ — non-trivial but feasible. It would also potentially yield an analytic expression for $\beta$ rather than requiring simulation calibration.

### 9.3 Analytic Bound on β

A perturbation analysis of the Flatto-Hahn generating function near $r = 1$ could determine the slope $\partial h/\partial r|_{r=1}$, which fixes $\beta$ analytically. This would fully close the approximation without simulation.

---

## 10. References

- Nelson, R. and Tantawi, A.N. (1988). "Approximate analysis of fork/join synchronization in parallel queues." *IEEE Transactions on Computers*, 37(6), 739–743.
- Flatto, L. and Hahn, S. (1984). "Two parallel queues created by arrivals with two demands I." *SIAM Journal on Applied Mathematics*, 44(5), 1041–1053.
- Reiman, M.I. and Simon, B. (1988). "Light traffic limits of sojourn time distributions in Markovian queueing networks." *Stochastic Models*, 4(2), 191–233.
- Squillante, M. and Tantawi, A.N. (2026). "QC/HPC Resource Allocation/Scheduling." Technical report. Appendix A.
