"""Expected sojourn time for the sequential computational phase, Eq. (2).

Second version: four rho curves on a single plot (linear y-axis).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
import numpy as np

RHO_VALS   = [0.4, 0.8, 0.9, 0.95]
R_DISCRETE = np.array([1, 2, 4, 8])
R_CONT     = np.geomspace(1, 8, 400)
P_O        = 0.5

COLORS  = ["steelblue", "darkorange", "forestgreen", "crimson"]
MARKERS = ["o", "D", "s", "^"]


def E_T_PF(rho, r):
    g = P_O * rho
    return 2.0 / (1.0 - g) + 2.0 / (r - g)


fig, ax = plt.subplots(figsize=(7, 5))

for rho, color, marker in zip(RHO_VALS, COLORS, MARKERS):
    ax.plot(R_CONT, E_T_PF(rho, R_CONT),
            color=color, linewidth=2.0, zorder=2)
    ax.plot(R_DISCRETE, E_T_PF(rho, R_DISCRETE),
            marker, color=color, markersize=6, zorder=3)

ax.set_xscale("log", base=2)
ax.set_xticks(R_DISCRETE)
ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
ax.tick_params(axis="x", which="minor", bottom=False)
ax.set_ylim(2, 8)
ax.grid(True, linestyle="--", linewidth=0.6, color="lightgray", zorder=0)

ax.set_xlabel(r"$r = \hat{\mu}_2 / \hat{\mu}_1$", fontsize=11)
ax.set_ylabel(r"$\mathbb{E}[T_{PF}]$", fontsize=11)

legend_handles = [
    Line2D([0], [0], color=c, linewidth=2.0, marker=m, markersize=6,
           label=rf"$\rho = {rho}$")
    for rho, c, m in zip(RHO_VALS, COLORS, MARKERS)
]
ax.legend(handles=legend_handles, fontsize=10, loc="upper right")

ax.set_title(
    r"Sequential Phase Sojourn Time $\mathbb{E}[T_{PF}]$ (Eq. 2)"
    "\n"
    r"$\hat{\mu}_1 = \hat{\mu}_3 = 1,\;"
    r"\hat{\mu}_2 = \hat{\mu}_4 = r,\;"
    r"p_\circ = 0.5$",
    fontsize=11,
)

out_dir = Path(__file__).parent / "figures"
out_dir.mkdir(exist_ok=True)
fig.savefig(out_dir / "t_pf_sequential_phase_v2.pdf", bbox_inches="tight")
fig.savefig(out_dir / "t_pf_sequential_phase_v2.png", dpi=150, bbox_inches="tight")
print(f"Saved to {out_dir}")
plt.show()
