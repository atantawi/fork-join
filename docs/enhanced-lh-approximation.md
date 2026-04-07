# Enhanced Light-Heavy Traffic Interpolation Approximation ($T_{\text{LH}}^{\text{enh}}$)

## 1. Motivation

The standard $T_{\text{LH}}$ approximation (Section 5 of [`heterogeneous-fj-approximations.md`](heterogeneous-fj-approximations.md)) uses a zero-order Reiman-Simon framework: a linear numerator with two coefficients matched to the zeroth light-traffic derivative $f^{(0)}(0)$ and the heavy-traffic constant $h$. This leaves one degree of freedom unused — the first light-traffic derivative $f^{(1)}(0)$ — which limits accuracy in the near-homogeneous, heavy-load regime where $T_{\text{LH}}$ reaches its worst-case error of $\approx +2.4\%$.

Appendix A of Squillante and Tantawi [1] derives closed-form expressions for both $f^{(0)}(0)$ and $f^{(1)}(0)$ for the heterogeneous two-queue fork-join system (Lemmas 1 and 2). Incorporating the first derivative into the Reiman-Simon framework yields a **first-order** approximation with a quadratic numerator that more faithfully captures the system's behavior across load levels.

---

## 2. The Reiman-Simon First-Order Framework

Following [1, Appendix A], we approximate the mean response time $f(\lambda) = T(\lambda)$ by scaling:

$$g(\lambda) = (\mu^* - \lambda)\, f(\lambda), \qquad \mu^* = \mu_{\min},$$

and interpolating $g(\lambda)$ with a polynomial $\tilde{g}(\lambda)$ matched to $n+1 = 3$ conditions:

| Condition | Meaning |
|-----------|---------|
| $\tilde{g}(0) = g(0)$ | zeroth light-traffic derivative |
| $\tilde{g}'(0) = g'(0)$ | first light-traffic derivative |
| $\tilde{g}(\mu^*) = h$ | heavy-traffic limit |

The three conditions determine three coefficients of a quadratic $\tilde{g}(\lambda) = b_2\lambda^2 + b_1\lambda + b_0$, after which we invert the scaling to obtain $\tilde{f}(\lambda) = \tilde{g}(\lambda)/(\mu^* - \lambda)$.

---

## 3. Light-Traffic Derivatives

### 3.1 Zeroth Derivative — Lemma 1 of [1]

$$f^{(0)}(0) = T_0 = \frac{1}{\mu_1} + \frac{1}{\mu_2} - \frac{1}{\mu_1 + \mu_2} = \frac{\mu_1^2 + \mu_1\mu_2 + \mu_2^2}{\mu_1\mu_2(\mu_1+\mu_2)}$$

This is identical to the light-traffic condition used by $T_{\text{LH}}$.

### 3.2 First Derivative — Lemma 2 of [1]

$$f^{(1)}(0) = \frac{1}{\mu_1^2} + \frac{1}{\mu_2^2} - \frac{2}{(\mu_1+\mu_2)^2} - \frac{2\mu_1\mu_2}{(\mu_1+\mu_2)^4}$$

This was derived in [1] by analyzing the tagged-job sojourn time conditioned on exactly one background arrival at time $t \in (-\infty, \infty)$, then integrating over all $t$.

---

## 4. The Enhanced Formula

Working in terms of $\rho = \lambda/\mu_{\min}$ and letting $\mu_{\min} = \min(\mu_1, \mu_2)$, $\mu_{\max} = \max(\mu_1, \mu_2)$, the three matching conditions give:

$$c_0 = \mu_{\min} \cdot T_0$$

$$c_1 = \mu_{\min}^2 \cdot f^{(1)}(0) - \mu_{\min} \cdot T_0$$

$$c_2 = h - c_1 - c_0$$

where $h = h(r) = 1 + \tfrac{3}{8}\,r^{-\beta}$ is the same heavy-traffic factor used by $T_{\text{LH}}$ (with $r = \mu_{\max}/\mu_{\min}$ and default $\beta = 10$). The enhanced approximation is then:

$$\boxed{T_{\text{LH}}^{\text{enh}} = \frac{c_2\,\rho^2 + c_1\,\rho + c_0}{\mu_{\min}(1 - \rho)}}$$

Explicitly, expanding $c_1$:

$$c_1 = \frac{\mu_{\min}^2}{\mu_{\max}^2} - \frac{\mu_{\min}}{\mu_{\max}} - \frac{2\mu_{\min}^2}{(\mu_{\min}+\mu_{\max})^2} + \frac{\mu_{\min}}{\mu_{\min}+\mu_{\max}} - \frac{2\mu_{\min}^3\mu_{\max}}{(\mu_{\min}+\mu_{\max})^4}$$

---

## 5. Verification: Homogeneous Case

When $\mu_1 = \mu_2 = \mu$ (so $r = 1$, $h = 11/8$, $T_0 = 3/(2\mu)$):

$$f^{(1)}(0) = \frac{2}{\mu^2} - \frac{2}{4\mu^2} - \frac{2\mu^2}{16\mu^4} = \frac{2}{\mu^2} - \frac{1}{2\mu^2} - \frac{1}{8\mu^2} = \frac{11}{8\mu^2}$$

$$c_0 = \frac{3}{2}, \qquad c_1 = \mu^2 \cdot \frac{11}{8\mu^2} - \mu \cdot \frac{3}{2\mu} = \frac{11}{8} - \frac{3}{2} = -\frac{1}{8}$$

$$c_2 = \frac{11}{8} - \left(-\frac{1}{8}\right) - \frac{3}{2} = \frac{12}{8} - \frac{3}{2} = 0$$

With $c_2 = 0$ the formula degenerates to:

$$T_{\text{LH}}^{\text{enh}} = \frac{3/2 - \rho/8}{\mu(1-\rho)} = \frac{12-\rho}{8(\mu-\lambda)} = T_2^{\text{hom}} \quad \checkmark$$

The enhancement vanishes exactly in the homogeneous case; the first derivative provides no additional information there.

---

## 6. Relationship to $T_{\text{LH}}$ and Equation 60 of [1]

Theorem 5 of [1] (equation 60) gives a first-order approximation that uses $f^{(0)}(0)$ and $f^{(1)}(0)$ **but not $h$** — it sets the third coefficient to zero by construction ($c_2 = 0$ in the document's derivation because $a_2 = 0$ is imposed). The result is a linear numerator whose coefficients are determined solely by the two light-traffic derivatives:

$$\tilde{f}_{\text{eq60}}(\lambda) = \frac{\lambda}{\mu_{\min}-\lambda}\,D + T_0, \qquad D = \frac{1}{\mu_{\min}} + \frac{\mu_{\min}}{\mu_{\max}^2} - \frac{2\mu_{\min}}{(\mu_{\min}+\mu_{\max})^2} - \frac{2\mu_{\min}^2\mu_{\max}}{(\mu_{\min}+\mu_{\max})^4}$$

$T_{\text{LH}}^{\text{enh}}$ is the **full** first-order approximation that uses all three conditions simultaneously. It has the same linear numerator form as $T_{\text{LH}}$ but with $a_1$ informed by the first derivative rather than solely by $h$, and with $c_2 \neq 0$ in the heterogeneous case to absorb any mismatch at the heavy-traffic limit.

| Approximation | Light-traffic ($f^{(0)}$) | Light-traffic ($f^{(1)}$) | Heavy-traffic ($h$) | Numerator degree |
|---------------|:---:|:---:|:---:|:---:|
| $T_{\text{LH}}$ | ✓ | — | ✓ | 1 |
| Eq. 60 of [1] | ✓ | ✓ | — | 1 |
| $T_{\text{LH}}^{\text{enh}}$ | ✓ | ✓ | ✓ | **2** |

---

## 7. Numerical Results

Results from discrete-event simulation (2,000,000 jobs, 100,000 warmup, seed 42). Errors are $(T_{\text{approx}} - T_{\text{sim}})/T_{\text{sim}} \times 100\%$.

| $\mu_1$ | $\mu_2$ | $\lambda$ | $\rho_1$ | $T_{\text{sim}}$ | $T_{\text{LH}}$ | Err% | $T_{\text{LH}}^{\text{enh}}$ | Err% |
|:------:|:------:|:--------:|:-------:|:---------------:|:---------------:|:----:|:----------------------------:|:----:|
| 1.0 | 1.0 | 0.3 | 0.30 | 2.088 | 2.089 | +0.07% | 2.089 | +0.07% |
| 1.0 | 1.0 | 0.6 | 0.60 | 3.564 | 3.562 | −0.05% | 3.562 | −0.05% |
| 1.0 | 1.0 | 0.9 | 0.90 | 13.886 | 13.875 | −0.08% | 13.875 | −0.08% |
| 1.0 | 1.5 | 0.3 | 0.30 | 1.705 | 1.698 | −0.40% | 1.710 | +0.32% |
| 1.0 | 1.5 | 0.6 | 0.60 | 2.767 | 2.776 | +0.35% | 2.801 | +1.25% |
| 1.0 | 1.5 | 0.9 | 0.90 | 10.083 | 10.325 | +2.41% | 10.362 | +2.77% |
| 1.0 | 2.0 | 0.3 | 0.30 | 1.582 | 1.595 | +0.87% | 1.589 | +0.45% |
| 1.0 | 2.0 | 0.6 | 0.60 | 2.623 | 2.667 | +1.70% | 2.654 | +1.20% |
| 1.0 | 2.0 | 0.9 | 0.90 | 9.990 | 10.170 | +1.80% | 10.150 | +1.60% |
| 1.0 | 3.0 | 0.6 | 0.60 | 2.546 | 2.583 | +1.45% | 2.561 | +0.57% |
| 1.0 | 5.0 | 0.6 | 0.60 | 2.514 | 2.533 | +0.75% | 2.519 | +0.19% |

---

## 8. Analysis of Results

### 8.1 Homogeneous Case ($\mu_1 = \mu_2 = 1.0$)

Both approximations are identical — the quadratic coefficient $c_2 = 0$ as shown in Section 5 — and errors remain below $0.1\%$ (attributable to finite simulation variance). Both recover the Nelson-Tantawi formula exactly by construction.

### 8.2 Moderate Heterogeneity ($r = 1.5$)

Results are mixed. At light-to-moderate load ($\lambda \in \{0.3, 0.6\}$), the enhanced approximation is in the same ballpark as $T_{\text{LH}}$ but slightly worse. At heavy load ($\lambda = 0.9$), the enhancement makes little difference (+2.41% vs +2.77%). This regime — near-homogeneous and heavily loaded — is precisely where the approximation of $h(r)$ by the power-law $1 + \tfrac{3}{8}r^{-\beta}$ is least reliable, since the true $h$ varies sharply near $r = 1$. A more accurate determination of $h$ (e.g., via the perturbation approach of [2, §9.2]) would likely improve performance here.

### 8.3 High Heterogeneity ($r \geq 2$)

The enhancement provides clear improvement:

| Scenario | $T_{\text{LH}}$ error | $T_{\text{LH}}^{\text{enh}}$ error | Improvement |
|----------|----------------------|------------------------------------|-------------|
| $r=2$, $\lambda=0.3$ | +0.87% | +0.45% | 48% |
| $r=2$, $\lambda=0.6$ | +1.70% | +1.20% | 29% |
| $r=2$, $\lambda=0.9$ | +1.80% | +1.60% | 11% |
| $r=3$, $\lambda=0.6$ | +1.45% | +0.57% | 61% |
| $r=5$, $\lambda=0.6$ | +0.75% | +0.19% | 75% |

The improvement is most pronounced at high heterogeneity and light-to-moderate load — exactly the regime where the additional information in $f^{(1)}(0)$ is most useful (the system is far from saturation and the queue dynamics are governed primarily by service time distributions rather than queueing delays).

### 8.4 Summary of Accuracy

| Regime | $T_{\text{LH}}$ max error | $T_{\text{LH}}^{\text{enh}}$ max error |
|--------|--------------------------|----------------------------------------|
| Homogeneous ($r=1$) | ≤ 0.08% | ≤ 0.08% |
| Moderate heterogeneity ($r=1.5$) | ≤ 2.41% | ≤ 2.77% |
| High heterogeneity ($r \geq 2$) | ≤ 1.80% | ≤ 1.60% |
| All scenarios | ≤ 2.41% | ≤ 2.77% |

The enhanced approximation improves on $T_{\text{LH}}$ for $r \geq 2$ but offers no benefit (and a slight degradation) at $r = 1.5$. The overall worst-case error is comparable.

---

## 9. Properties

| Property | Status |
|----------|--------|
| Exact for homogeneous case ($\mu_1 = \mu_2$) | ✓ (for any $\beta$) |
| Light traffic limit ($\lambda \to 0$): $T \to T_0$ | ✓ |
| First light-traffic derivative matched | ✓ |
| Heavy traffic limit: $T \sim h(r)/(\mu_{\min}-\lambda)$ | ✓ |
| Bound compliance $T_{\text{bot}} \leq T \leq T_{\text{UB}}$ | Not guaranteed |
| Tunable parameter | $\beta$ (default 10) |
| Numerator degree | 2 (quadratic) |

---

## 10. Implementation

The formula is implemented as `mean_response_time_lh_enhanced(lam, mu1, mu2, beta=10.0)` in `forkjoin/analytical.py` and exported from `forkjoin/__init__.py`.

```python
from forkjoin import mean_response_time_lh_enhanced

T = mean_response_time_lh_enhanced(lam=0.6, mu1=1.0, mu2=2.0)
```

---

## References

[1] M. Squillante and A. N. Tantawi, "QC/HPC Resource Allocation/Scheduling," Technical Report, January 2026. (Theorems 4–5, Lemmas 1–2 in Appendix A.)

[2] R. Nelson and A. N. Tantawi, "Approximate analysis of fork/join synchronization in parallel queues," *IEEE Transactions on Computers*, vol. 37, no. 6, pp. 739–743, 1988.

[3] M. I. Reiman and B. Simon, "Light traffic limits of sojourn time distributions in Markovian queueing networks," *Stochastic Models*, vol. 4, no. 2, pp. 191–233, 1988.
