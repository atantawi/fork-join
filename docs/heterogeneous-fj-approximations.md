# Approximate Methods for Heterogeneous 2-Queue Fork-Join Systems

## 1. Problem Statement

Consider a **fork-join (FJ) queueing system** with two parallel servers having potentially different service rates $\mu_1$ and $\mu_2$. Without loss of generality, assume $\mu_1 \geq \mu_2$, so server 2 is the bottleneck. Jobs arrive as a Poisson process with rate $\lambda$. Each arriving job forks into two tasks that are dispatched simultaneously, one to each server queue; the job departs only after **both** tasks have completed service. Each server operates as an independent M/M/1 FCFS queue.

Define the following notation:

| Symbol | Meaning |
|--------|---------|
| $\lambda$ | Poisson arrival rate |
| $\mu_i$ | Exponential service rate at server $i$ |
| $\mu_{\min} = \mu_2$ | Bottleneck service rate ($\mu_1 \geq \mu_2$) |
| $\mu_{\max} = \mu_1$ | Faster service rate |
| $r = \mu_{\max}/\mu_{\min} \geq 1$ | Heterogeneity ratio |
| $\rho_i = \lambda/\mu_i$ | Utilization of server $i$ |
| $R_i$ | Sojourn time (waiting + service) of a task at server $i$ |
| $T$ | Mean job response time $= E[\max(R_1, R_2)]$ |

**Stability condition:** $\lambda < \mu_{\min}$, equivalently $\rho_{\min} = \lambda/\mu_{\min} < 1$.

**Goal:** Find a closed-form or easily computable expression for $T$ that is accurate across all stable operating points.

### 1.1 Why This is Hard

The fork-join constraint — that a job waits for the **maximum** of two dependent sojourn times — makes exact analysis difficult. The sojourn times $R_1$ and $R_2$ are not independent: they share the same arrival process, inducing positive correlation. As a result:

- Simple independence assumptions (treating the two servers as decoupled M/M/1 queues) **overestimate** $T$, since positive correlation reduces the expected maximum.
- The exact joint distribution of $(R_1, R_2)$ is not available in closed form for the heterogeneous case. Flatto and Hahn [1] derived an exact analysis via elliptic function parametrization, but it yields a complex generating function rather than a tractable expression for $T$.

### 1.2 Known Exact Result: Homogeneous Case

When $\mu_1 = \mu_2 = \mu$ (with $\rho = \lambda/\mu$), Nelson and Tantawi [2] derived the exact result:

$$T_2^{\text{hom}} = \frac{12 - \rho}{8(\mu - \lambda)}$$

This serves as a critical calibration point for any approximation.

---

## 2. Theoretical Bounds

### 2.1 Independent Upper Bound

Baccelli, Makowski, and Shwartz [3] showed that if $R_1$ and $R_2$ were **independent** M/M/1 sojourn times (which they are not — the shared arrival process induces positive correlation), then $E[\max(R_1, R_2)]$ would be:

$$T_{\text{UB}} = \frac{1}{\mu_1 - \lambda} + \frac{1}{\mu_2 - \lambda} - \frac{1}{\mu_1 + \mu_2 - 2\lambda}$$

Since positive correlation reduces the expected maximum, this is an **upper bound**: $T \leq T_{\text{UB}}$.

### 2.2 Bottleneck Lower Bound

By Jensen's inequality, $E[\max(R_1, R_2)] \geq \max(E[R_1], E[R_2])$, giving:

$$T_{\text{bot}} = \max\!\left(\frac{1}{\mu_1 - \lambda},\; \frac{1}{\mu_2 - \lambda}\right) = \frac{1}{\mu_{\min} - \lambda}$$

This is the mean sojourn time of the bottleneck M/M/1 queue alone, and serves as a **lower bound**: $T \geq T_{\text{bot}}$.

### 2.3 Gap Between Bounds

The gap $T_{\text{UB}} - T_{\text{bot}}$ shrinks as heterogeneity increases. In the homogeneous case ($r=1$), the bounds satisfy $T_{\text{bot}} = T_1 = 1/(\mu - \lambda)$ and $T_{\text{UB}} = \frac{3}{2}T_1$, so the bounds differ by a factor of $3/2$. For large $r$, both bounds collapse toward $T_{\text{bot}}$.

---

## 3. Two Approximation Methods

We present two complementary closed-form approximations. Both reduce to the Nelson-Tantawi result in the homogeneous case and respect the correct limiting behavior in light and heavy traffic.

---

## 4. Method 1: Upper-Lower Bound Interpolation ($T_{\text{UL}}$)

### 4.1 Strategy

Interpolate between the upper and lower bounds using the average utilization as the weight:

$$T_{\text{UL}} = (1 - \alpha)\,T_{\text{UB}} + \alpha\,T_{\text{bot}}$$

The weight $\alpha$ must be chosen to:
1. Recover the Nelson-Tantawi result in the homogeneous case
2. Approach $T_{\text{UB}}$ in light traffic ($\lambda \to 0$, nearly no queueing)
3. Approach $T_{\text{bot}}$ in heavy traffic ($\lambda \to \mu_{\min}$, bottleneck dominates)

### 4.2 Choosing the Weight

Define the average utilization $\bar{\rho} = (\rho_1 + \rho_2)/2$ and set $\alpha = \bar{\rho}/4$:

$$\alpha = \frac{\rho_1 + \rho_2}{8}$$

**Verification in the homogeneous case** ($\mu_1 = \mu_2 = \mu$, $\rho = \lambda/\mu$):

- $T_{\text{UB}} = \frac{3}{2} T_1$, $\quad T_{\text{bot}} = T_1$, $\quad \alpha = \rho/4$
- $T_{\text{UL}} = \left(1 - \frac{\rho}{4}\right) \frac{3}{2} T_1 + \frac{\rho}{4} T_1 = \frac{3/2 - 3\rho/8 + \rho/4}{1} T_1 = \frac{12 - \rho}{8} T_1$ ✓

### 4.3 The Formula

$$\boxed{T_{\text{UL}} = \left(1 - \frac{\rho_1+\rho_2}{8}\right) T_{\text{UB}} + \frac{\rho_1+\rho_2}{8}\, T_{\text{bot}}}$$

### 4.4 Properties

| Property | Status |
|----------|--------|
| Exact for homogeneous case | ✓ by construction |
| Light traffic limit ($\lambda \to 0$) | ✓ $\alpha \to 0$, so $T_{\text{UL}} \to T_{\text{UB}}$ |
| Heavy traffic limit ($\lambda \to \mu_{\min}$) | ✓ $T_{\text{UL}} / T_{\text{bot}} \to 1$ |
| Bound compliance $T_{\text{bot}} \leq T_{\text{UL}} \leq T_{\text{UB}}$ | ✓ since $0 \leq \alpha \leq 1/4 < 1$ |
| Closed-form | ✓ |

---

## 5. Method 2: Light-Heavy Traffic Interpolation ($T_{\text{LH}}$)

### 5.1 Strategy: The Reiman-Simon Framework

Following the methodology developed by Reiman and Simon [4] and described in [5, §2.1], we approximate $T(\lambda)$ as a rational function in the bottleneck utilization $\rho = \lambda/\mu_{\min}$:

$$\tilde{T}(\rho) = \frac{a_0 + a_1\,\rho}{\mu_{\min}(1 - \rho)}$$

The denominator $\mu_{\min}(1-\rho) = \mu_{\min} - \lambda$ encodes the correct heavy-traffic singularity. The two coefficients $a_0$ and $a_1$ are determined by matching two conditions: the exact light-traffic value and the heavy-traffic limit.

### 5.2 Light-Traffic Matching

When $\lambda = 0$ (empty system), response time equals the expected maximum of two independent exponential service times:

$$T_0 \equiv T(\rho=0) = \frac{1}{\mu_1} + \frac{1}{\mu_2} - \frac{1}{\mu_1 + \mu_2} = E\!\left[\max\!\left(X_1, X_2\right)\right]$$

where $X_i \sim \text{Exp}(\mu_i)$. The formula $E[\max(X_1, X_2)] = E[X_1] + E[X_2] - E[\min(X_1, X_2)]$ is exact; $\min(X_1, X_2)$ is exponential with rate $\mu_1 + \mu_2$. Setting $\tilde{T}(0) = T_0$ gives:

$$a_0 = \mu_{\min} \cdot T_0$$

### 5.3 Heavy-Traffic Matching

As $\lambda \to \mu_{\min}$ (i.e., $\rho \to 1$), define the heavy-traffic constant:

$$h = \lim_{\rho \to 1} \mu_{\min}(1 - \rho) \cdot T(\rho)$$

Matching this in the approximation gives $a_0 + a_1 = h$, hence $a_1 = h - a_0$.

### 5.4 The Heavy-Traffic Factor $h(r)$

The value of $h$ depends on the heterogeneity ratio $r = \mu_{\max}/\mu_{\min}$:

**Homogeneous case** ($r = 1$): Both servers saturate simultaneously as $\lambda \to \mu$. The correlation between sojourn times remains significant, and the Nelson-Tantawi formula gives $h = 11/8$.

**Heterogeneous case** ($r > 1$): As $\lambda \to \mu_{\min}$, the bottleneck server saturates while the faster server's utilization approaches $\lambda/\mu_{\max} = \mu_{\min}/\mu_{\max} = 1/r < 1$. Its sojourn time remains bounded ($\to 1/(\mu_{\max} - \mu_{\min})$), so the bottleneck dominates and $h \to 1$.

We model the transition with the **power-law factor**:

$$\boxed{h(r) = 1 + \frac{3}{8}\, r^{-\beta}, \quad r = \frac{\mu_{\max}}{\mu_{\min}} \geq 1, \quad \beta > 0}$$

This gives:
- $h(1) = 1 + 3/8 = 11/8$ for all $\beta$ ✓
- $h(r) \to 1$ as $r \to \infty$ for all $\beta > 0$ ✓
- Monotonically decreasing
- Preserves the rational-function structure of $\tilde{T}(\rho)$
- One calibratable parameter $\beta$

**Physical calibration of $\beta$:** Simulation data (Section 7) shows that the transition from $h = 11/8$ to $h \approx 1$ is **extremely sharp**. Even at $r = 1.5$ (a 3:2 speed ratio), the effective $h$ is already approximately $1.008$. This implies $\beta \gg 1$; numerical fitting to the available simulation data gives $\beta \approx 10$ as a reasonable default.

### 5.5 The Formula

Substituting $a_0$ and $a_1 = h(r) - a_0$ into the Reiman-Simon form:

$$\boxed{T_{\text{LH}} = \frac{\mu_{\min} T_0 + \bigl(h(r) - \mu_{\min} T_0\bigr)\,\rho}{\mu_{\min} - \lambda}}$$

where $T_0 = 1/\mu_1 + 1/\mu_2 - 1/(\mu_1+\mu_2)$, $\rho = \lambda/\mu_{\min}$, $r = \mu_{\max}/\mu_{\min}$, and $h(r) = 1 + \tfrac{3}{8}\,r^{-\beta}$.

**Verification — homogeneous case** ($\mu_1 = \mu_2 = \mu$, $r = 1$):

$$a_0 = \mu \cdot \frac{3}{2\mu} = \frac{3}{2}, \qquad h(1) = \frac{11}{8}, \qquad a_1 = \frac{11}{8} - \frac{3}{2} = -\frac{1}{8}$$

$$T_{\text{LH}} = \frac{3/2 - \rho/8}{\mu(1-\rho)} = \frac{12-\rho}{8(\mu-\lambda)} = T_2^{\text{hom}} \quad \checkmark$$

This recovery holds for **all** $\rho$, regardless of $\beta$.

### 5.6 Properties

| Property | Status |
|----------|--------|
| Exact for homogeneous case | ✓ by construction (for any $\beta$) |
| Light traffic limit ($\lambda \to 0$) | ✓ $T_{\text{LH}} \to T_0$ |
| Heavy traffic limit ($\lambda \to \mu_{\min}$) | ✓ $T_{\text{LH}} \sim h(r)/(\mu_{\min} - \lambda)$ |
| Bound compliance | Not guaranteed |
| Theoretical basis | Reiman-Simon rational interpolation [4,5] |
| Tunable parameter | $\beta$ (default 10, calibrated to simulation) |

---

## 6. Comparison of the Two Methods

| Aspect | $T_{\text{UL}}$ (Bound Interpolation) | $T_{\text{LH}}$ (Traffic Interpolation) |
|--------|--------------------------------------|----------------------------------------|
| **Homogeneous case** | Exact (by construction) | Exact (by construction) |
| **Light traffic** | $T \to T_{\text{UB}}$ | $T \to T_0$ (same limit) |
| **Heavy traffic** | $T/T_{\text{bot}} \to 1$ | $T \sim h(r)/(\mu_{\min}-\lambda)$ |
| **Bound compliance** | Always ($T_{\text{bot}} \leq T_{\text{UL}} \leq T_{\text{UB}}$) | Not guaranteed |
| **Theoretical basis** | Convex interpolation of bounds | Reiman-Simon rational approximation |
| **Parameters** | None | $\beta$ (default 10) |
| **Typical errors (heterogeneous)** | $\leq 1.2\%$ | $\leq 2.4\%$ (at $r = 1.5$) |

The methods are **complementary**:
- $T_{\text{UL}}$ provides guaranteed bounds and is more accurate at high heterogeneity.
- $T_{\text{LH}}$ is grounded in the Reiman-Simon framework and is exact at both traffic extremes by design; its $\beta$ parameter can be further refined via simulation.

---

## 7. Numerical Results

Results from discrete-event simulation (2,000,000 jobs per scenario, 100,000 warmup jobs, seed 42). All four quantities — $T_{\text{bot}}$, $T_{\text{UL}}$, $T_{\text{LH}}$ (with $\beta = 10$), and $T_{\text{UB}}$ — are compared against $T_{\text{sim}}$.

| $\mu_1$ | $\mu_2$ | $\lambda$ | $\rho_2$ | $T_{\text{bot}}$ | $T_{\text{sim}}$ | $T_{\text{UL}}$ | Err% | $T_{\text{LH}}$ | Err% | $T_{\text{UB}}$ |
|:------:|:------:|:--------:|:-------:|:---------------:|:---------------:|:---------------:|:----:|:---------------:|:----:|:---------------:|
| 1.0 | 1.0 | 0.3 | 0.30 | 1.429 | 2.088 | 2.089 | +0.06% | 2.089 | +0.06% | 2.143 |
| 1.0 | 1.0 | 0.6 | 0.60 | 2.500 | 3.564 | 3.562 | −0.04% | 3.562 | −0.04% | 3.750 |
| 1.0 | 1.0 | 0.9 | 0.90 | 10.000 | 13.886 | 13.875 | −0.08% | 13.875 | −0.08% | 15.000 |
| 1.0 | 1.5 | 0.3 | 0.20 | 1.429 | 1.705 | 1.716 | +0.67% | 1.698 | −0.41% | 1.736 |
| 1.0 | 1.5 | 0.6 | 0.40 | 2.500 | 2.767 | 2.799 | +1.16% | 2.776 | +0.34% | 2.842 |
| 1.0 | 1.5 | 0.9 | 0.60 | 10.000 | 10.083 | 10.193 | +1.10% | 10.325 | +2.40% | 10.238 |
| 1.0 | 2.0 | 0.3 | 0.15 | 1.429 | 1.582 | 1.590 | +0.54% | 1.595 | +0.85% | 1.600 |
| 1.0 | 2.0 | 0.6 | 0.30 | 2.500 | 2.623 | 2.641 | +0.68% | 2.667 | +1.69% | 2.659 |
| 1.0 | 2.0 | 0.9 | 0.45 | 10.000 | 9.990 | 10.063 | +0.73% | 10.170 | +1.80% | 10.076 |
| 1.0 | 3.0 | 0.6 | 0.20 | 2.500 | 2.546 | 2.554 | +0.30% | 2.583 | +1.47% | 2.560 |
| 1.0 | 5.0 | 0.6 | 0.12 | 2.500 | 2.514 | 2.517 | +0.11% | 2.533 | +0.77% | 2.519 |

### 7.1 Observations

**Homogeneous case** ($\mu_1 = \mu_2 = 1.0$):
Both $T_{\text{UL}}$ and $T_{\text{LH}}$ are essentially exact (errors $\leq 0.1\%$). This is expected: both are designed to recover the Nelson-Tantawi result exactly in this regime.

**Moderate heterogeneity** ($r = 1.5$, i.e., $\mu_2 = 1.5$):
- $T_{\text{UL}}$: errors $+0.67\%$ to $+1.16\%$ (slight systematic overestimate)
- $T_{\text{LH}}$: errors $-0.41\%$ to $+2.40\%$ (mixed sign; worst at heavy load)

The $T_{\text{LH}}$ error at $(r=1.5, \rho=0.9)$ of $+2.40\%$ reflects the inherent limitation of the linear Reiman-Simon numerator at near-homogeneous, heavy-load operating points — a regime where the two effects (heterogeneity and load) interact. The parameter $\beta$ has limited ability to resolve this without a higher-order (quadratic) approximation (see Section 9.1).

**High heterogeneity** ($r \geq 2$):
- $T_{\text{UL}}$: errors $\leq 0.73\%$, improving to $\leq 0.11\%$ at $r = 5$
- $T_{\text{LH}}$: errors $\leq 1.80\%$

As heterogeneity increases, both $T_{\text{UB}}$ and $T_{\text{LH}}$ converge to $T_{\text{bot}}$ from above, and $T_{\text{UL}}$ tracks the simulation closely. In this regime the faster server's contribution to $E[\max(R_1, R_2)]$ is negligible; all that matters is the bottleneck M/M/1.

---

## 8. Summary

The two approximation methods address the lack of a closed-form mean response time for the heterogeneous 2-queue fork-join system:

1. **$T_{\text{UL}}$ (bound interpolation):** A convex combination of the independent upper bound and the bottleneck lower bound, weighted by the average utilization. Exact for the homogeneous case, always within bounds, with errors below $1.2\%$ across all tested scenarios. Recommended as the default.

2. **$T_{\text{LH}}$ (light-heavy traffic interpolation):** A rational Reiman-Simon approximation matching the light-traffic value exactly and the heavy-traffic limit via the factor $h(r) = 1 + \tfrac{3}{8}\,r^{-\beta}$. Exact for the homogeneous case. Maximum errors $\approx 2.4\%$ in tested scenarios. The parameter $\beta$ can be refined via simulation; a quadratic numerator (Section 9.1) could further improve accuracy.

Both methods are closed-form, computationally trivial, and calibrated to be exact in the homogeneous limit.

---

## 9. Future Directions

### 9.1 Quadratic Reiman-Simon Approximation for $T_{\text{LH}}$

The current $T_{\text{LH}}$ uses a linear numerator, matching two conditions. A quadratic numerator $a_0 + a_1\rho + a_2\rho^2$ would match three conditions: the two above plus the first light-traffic derivative $f'(0)$, which Reiman and Simon [4] express as:

$$\frac{df}{d\rho}\bigg|_{\rho=0} = \mu_{\min} \int_{-\infty}^{+\infty} \left[\hat{f}(\{t\}) - \hat{f}(\emptyset)\right] dt$$

where $\hat{f}(\emptyset) = T_0$ is the unconditional light-traffic value and $\hat{f}(\{t\})$ is the conditional expected response time given one prior arrival at time $t$. Computing this requires an analysis of the joint queue state conditioned on one background customer, which is non-trivial but feasible. Matching this derivative would reduce the approximation error in the near-homogeneous, heavy-load regime and may also determine $\beta$ analytically.

### 9.2 Analytic Determination of $\beta$

The sharp transition of $h(r)$ near $r = 1$ is a structural property of the heterogeneous fork-join system. A perturbation analysis of the Flatto-Hahn generating function [1] near $r = 1$ would yield $\partial h/\partial r|_{r=1}$, which directly determines $\beta$ (since $h'(1) = -3\beta/8$). This would fully close the approximation without requiring simulation calibration.

### 9.3 Extension to $K > 2$ Servers

Both methods extend naturally to larger systems:
- $T_{\text{UL}}$ can use the independent upper bound $T_{\text{UB}} = \sum_{i=1}^{K} 1/(\mu_i - \lambda) - \ldots$ (inclusion-exclusion) and the bottleneck lower bound.
- $T_{\text{LH}}$ can employ the same Reiman-Simon framework, with $T_0 = E[\max(X_1, \ldots, X_K)]$ for independent exponentials computed recursively.

---

## References

[1] L. Flatto and S. Hahn, "Two parallel queues created by arrivals with two demands I," *SIAM Journal on Applied Mathematics*, vol. 44, no. 5, pp. 1041–1053, 1984.

[2] R. Nelson and A. N. Tantawi, "Approximate analysis of fork/join synchronization in parallel queues," *IEEE Transactions on Computers*, vol. 37, no. 6, pp. 739–743, 1988.

[3] F. Baccelli, A. M. Makowski, and A. Shwartz, "The fork-join queue and related systems with synchronization constraints: stochastic ordering and computable bounds," *Advances in Applied Probability*, vol. 21, no. 3, pp. 629–660, 1989.

[4] M. I. Reiman and B. Simon, "Light traffic limits of sojourn time distributions in Markovian queueing networks," *Stochastic Models*, vol. 4, no. 2, pp. 191–233, 1988.

[5] M. Squillante and A. N. Tantawi, "QC/HPC Resource Allocation/Scheduling," Technical Report, January 2026. (Section 2.1.)

[6] S. Varma and A. M. Makowski, "Interpolation approximations for symmetric fork-join queues," *Performance Evaluation*, vol. 20, no. 1, pp. 245–265, 1994.
