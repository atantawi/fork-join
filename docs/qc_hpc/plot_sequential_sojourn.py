"""Expected sojourn time for the sequential computational phase, Eq. (2).

Model: open Jackson network with 4 queues.
  Queue 1 (QPU): mu_hat_1 = 1,  arrival rate gamma_1 = p_o * lambda
  Queue 2 (GPU): mu_hat_2 = r,  arrival rate gamma_2 = (1-p_o) * lambda
  Queue 3 (GPU): mu_hat_3 = 1,  arrival rate gamma_3 = p_o * lambda
  Queue 4 (QPU): mu_hat_4 = r,  arrival rate gamma_4 = (1-p_o) * lambda

With p_o = 0.5 all gammas equal 0.5*lambda = 0.5*rho (since mu_hat_1 = 1).
Eq. (2) simplifies to: E[T_PF] = 2/(1 - 0.5*rho) + 2/(r - 0.5*rho)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

RHO_VALS   = [0.4, 0.8, 0.9, 0.95]
R_DISCRETE = np.array([1, 2, 4, 8])
R_CONT     = np.geomspace(1, 8, 400)
P_O        = 0.5


def E_T_PF(rho, r):
    g = P_O * rho  # gamma_i = 0.5 * rho for all i (since p_o = 0.5, mu1 = 1)
    return 2.0 / (1.0 - g) + 2.0 / (r - g)


fig, axes = plt.subplots(2, 2, figsize=(9, 7))
fig.subplots_adjust(hspace=0.38, wspace=0.32)

for idx, rho in enumerate(RHO_VALS):
    ax = axes[idx // 2][idx % 2]

    ax.plot(R_CONT, E_T_PF(rho, R_CONT),
            color="steelblue", linewidth=2.0, zorder=2)
    ax.plot(R_DISCRETE, E_T_PF(rho, R_DISCRETE),
            "o", color="steelblue", markersize=6, zorder=3)

    ax.set_title(rf"$\rho = {rho}$", fontsize=12)
    ax.set_xscale("log", base=2)
    ax.set_xticks(R_DISCRETE)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.tick_params(axis="x", which="minor", bottom=False)
    ax.grid(True, linestyle="--", linewidth=0.6, color="lightgray", zorder=0)

    if idx >= 2:
        ax.set_xlabel(r"$r = \hat{\mu}_2 / \hat{\mu}_1$", fontsize=11)
    ax.set_ylabel(r"$\mathbb{E}[T_{PF}]$", fontsize=11)

fig.suptitle(
    r"Sequential Phase Sojourn Time $\mathbb{E}[T_{PF}]$ (Eq. 2)"
    "\n"
    r"$\hat{\mu}_1 = \hat{\mu}_3 = 1,\;"
    r"\hat{\mu}_2 = \hat{\mu}_4 = r,\;"
    r"p_\circ = 0.5$",
    fontsize=11,
)

out_dir = Path(__file__).parent / "figures"
out_dir.mkdir(exist_ok=True)
fig.savefig(out_dir / "t_pf_sequential_phase.pdf", bbox_inches="tight")
fig.savefig(out_dir / "t_pf_sequential_phase.png", dpi=150, bbox_inches="tight")
print(f"Saved to {out_dir}")
plt.show()
