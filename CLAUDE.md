# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Research project on **fork-join (FJ) queues** — a queueing model where arriving jobs split (fork) into n parallel tasks, each joining a separate single-server queue, and the job departs only after all n tasks complete. The primary performance measure is mean job response time T (average time from arrival to departure).

## Domain Context

- **Homogeneous FJ queue (n=2):** Exact result by Nelson & Tantawi (1988): `T = (12 - rho) / (8 * (mu - lambda))` where `rho = lambda/mu`
- **Heterogeneous FJ queue:** Servers have different rates mu_1, mu_2, ..., mu_n. No known closed-form for mean response time; exact analysis by Flatto & Hahn (1984) yields generating functions via elliptic function parametrization but not a simple expression
- **Stability condition:** For n=2 heterogeneous, `lambda < mu_1` and `lambda < mu_2`
- The research goal is to find closed-form or approximate expressions for T in the heterogeneous case

## Key Reference Papers (docs/references/)

| File | Authors | Year | Key Content |
|------|---------|------|-------------|
| Nelson1988.pdf | Nelson, Tantawi | 1988 | Scaling approximation for homogeneous FJ; exact T_2 formula |
| Flatto1984.pdf | Flatto, Hahn | 1984 | Exact analysis of 2-queue heterogeneous FJ via elliptic functions |
| Flatto1985.pdf | Flatto | 1985 | Continuation: limit laws for queue lengths |
| Baccelli1989.pdf | Baccelli, Makowski, Shwartz | 1989 | Bounds on FJ response time (upper: independent queues, lower: D/G/1) |
| Varma1994.pdf | Varma, Makowski | 1994 | Light/heavy traffic interpolation approximations for symmetric FJ |
| Nguyen1994.pdf | Nguyen | 1994 | Heterogeneous customer populations in FJ networks |
| Balsamo2002.pdf | Balsamo, Donatiello, Van Dijk | 1998 | Bound performance models for heterogeneous parallel systems |
| Mohanty2024.pdf | Mohanty, Gautam, Aggarwal, Parag | 2024 | (k,k) FJ on heterogeneous servers; asymptotic independence; upper/lower bounds on mean completion time |

## Key Formulas

- **M/M/1 mean response time:** `T_1 = 1 / (mu - lambda)`
- **Homogeneous FJ (n=2):** `T_2 = (12 - rho) / 8 * T_1`
- **Independent upper bound (n=K):** `T_K <= H_K * T_1` where `H_K` is the K-th harmonic number
- **Scaling approximation (K>=2):** `T_K = [H_K/H_2 + (4/11)(1 - H_K/H_2) * rho] * T_2`
- **T_UL approximation:** `T_UL = (1 - alpha)*T_UB + alpha*T_bot` where `alpha = (rho_1+rho_2)/8`
- **T_LH approximation:** `T_LH = (a0 + a1*rho) / (mu_min - lambda)` where `a0 = mu_min*T0`, `a1 = h(r) - a0`, `h(r) = 1 + 3/8 * r^(-beta)`, `r = mu_max/mu_min`, `beta=10`
- **T_LH_enhanced (first-order):** `T_LHe = (c2*rho^2 + c1*rho + c0) / (mu_min*(1-rho))` where `c0 = mu_min*T0`, `c1 = mu_min^2*f1 - mu_min*T0`, `c2 = h - c1 - c0`, `f1 = 1/mu_min^2 + 1/mu_max^2 - 2/(mu_min+mu_max)^2 - 2*mu_min*mu_max/(mu_min+mu_max)^4`
- The Flatto-Hahn generating function for the heterogeneous 2-queue case uses `P(z,0)` expressed via `sqrt(a_3 - z)` where `a_3` is a root of the discriminant `D_1(z) = [(1+alpha+beta)z - alpha]^2 - 4*beta*z^3` (with alpha=mu_1/lambda, beta=mu_2/lambda)

## Notation Conventions

- `lambda`: Poisson job arrival rate
- `mu`, `mu_i`: exponential service rate(s)
- `mu_min = min(mu_1, mu_2)`, `mu_max = max(mu_1, mu_2)`: bottleneck and faster server rates
- `rho = lambda/mu`: server utilization (homogeneous case)
- `rho_i = lambda/mu_i`: utilization of server i (heterogeneous case)
- `r = mu_max/mu_min`: heterogeneity ratio (r=1 homogeneous, r→∞ highly heterogeneous)
- `T`, `T_K`: mean job response time (K queues)
- `n` or `K`: number of parallel queues/servers
- `H_K`: harmonic number sum(1/i, i=1..K)
- `T0 = 1/mu_1 + 1/mu_2 - 1/(mu_1+mu_2)`: light-traffic limit (E[max(X_1,X_2)] at lambda=0)
- `h(r) = 1 + 3/8 * r^(-beta)`: heavy-traffic factor; h(1)=11/8, h→1 as r→∞
